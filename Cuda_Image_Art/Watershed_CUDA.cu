// Copyright (c) 2023-2026 Steve Young
// Licensed under the MIT License

# include "Watershed_CUDA.cuh"
# include "Utils_CUDA.cuh"

__global__ void update_values(int width, int height, int* pixel_data, bool* ResultsArray)
{
	// In the first part of the watershed processing, new identifier values are written as negative (to disambiguate them as having changed on that step).  Now,
	// we just need to change those negatives to positives.
	// width and height are the dimensions of the pixel_data matrix.
	// pixel_data are the superpixel identifiers.
	// ResultsArray is an array with result information telling us which blocks need cleaning up.
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (ResultsArray[blockIdx.x])
	{
		if (idx < width * height)
		{
			int x = idx % width;
			int y = idx / width;
			int pos = x + y * width;
			int value = pixel_data[pos];
			if (value < 0)
			{
				pixel_data[pos] = -value;
			}
		}
	}
}
__global__ void watershed_step(int level, int width, int height, int* pixel_data, unsigned char* gradient_data, bool* ResultsArray)
{
	// Kernel for computing watershed algorithm over full image.  Each call is one pass, and will need to be repeated until all ResultsArray values are false.
	// In this version, there is no shared memory.  If performance is not good, re-architect to use shared memory.
	// level is the current level for the watershed processing.
	// width and height define the dimensions of the image.
	// pixel_data is a pixel map of the superpixel identifiers (each an integer).  It starts all with all zeroes except the seeds for each superpixel.
	// gradient_data is the gradient information for the image.
	// ResultsArray is array of booleans indicating for each block whether another pass is needed at this level.

	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < width * height)
	{
		int x = idx % width;
		int y = idx / width;
		int pos = x + y * width;
		int value = 0;
		// Is this an open position on the pixel map?
		if (0 == pixel_data[pos])
		{
			// Is the gradient level at or below the target level?
			if (level >= gradient_data[pos])
			{
				int smallest_adjacent = -1; // This will be the smallest positive adjacent superpixel identifier.
				if (y > 0) // Can't look above when at the top of the image
				{
					value = pixel_data[pos - width];
					if ((value > 0) && ((smallest_adjacent < 0) || (value < smallest_adjacent)))
					{
						smallest_adjacent = value;
					}
				}
				if (y < height - 1) // Make sure not at bottom of image.
				{
					value = pixel_data[pos + width];
					if ((value > 0) && ((smallest_adjacent < 0) || (value < smallest_adjacent)))
					{
						smallest_adjacent = value;
					}
				}
				if (x > 0) // Make sure not at left side.
				{
					value = pixel_data[pos - 1];
					if ((value > 0) && ((smallest_adjacent < 0) || (value < smallest_adjacent)))
					{
						smallest_adjacent = value;
					}
				}
				if (x < width - 1) // Make sure not on the right side.
				{
					value = pixel_data[pos + 1];
					if ((value > 0) && ((smallest_adjacent < 0) || (value < smallest_adjacent)))
					{
						smallest_adjacent = value;
					}
				}
				if (smallest_adjacent > 0) // This will be true if there is an adjacent value larger than zero, and it will always defer to the smaller identifier.
				{
					pixel_data[pos] = -smallest_adjacent; // Set to negative, so other threads don't read this as a superpixel until after this iteration is complete.
					ResultsArray[blockIdx.x] = true; // Having expanded a superpixel, another pass is needed at this level.
				}
			}
		}
	}
}

__global__ void calc_bbox(int width, int height, int* pixel_data, int* d_bbox_array, int sp_count)
{
	// This is the kernel for calculating the bounding box and size for each superpixel.
	// Each thread corresponds to one superpixel identifier.
	// width and height are the dimensions of the pixel_data matrix.
	// d_bbox_array is an array that will be filled out with the bounding box and size for each
	// superpixel.  The order is: identifier, x0, y0, x1, y1, size.
	// sp_count is the number of superpixels.

	extern __shared__ int row_data[];
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	// Set local variables for this superpixel, set to values outside of the actual matrix.
	int x0 = width;
	int y0 = height;
	int x1 = -1;
	int y1 = -1;
	int size = 0;

	// Get identifier for this superpixel.
	int identifier;
	if (idx < sp_count)
	{
		identifier = d_bbox_array[idx * 6];
	}

	// Loop through all rows.
	for (int y = 0; y < height; ++y)
	{
		// First task is to load in the row data.
		for (int x = threadIdx.x; x < width; x += threads_per_block) // Stride to gather all elements of the row.
		{
			row_data[x] = pixel_data[x + y * width];
		}
		__syncthreads(); // Need to ensure that all row data has been copied, so wait for all threads to reach this point.

		if (idx < sp_count)
		{
			bool identifier_seen = false;
			for (int x = 0; x < width; ++x)
			{
				if (identifier == row_data[x])
				{
					identifier_seen = true;
					size++;
					// Update horizontal values for bounding box.
					if (x < x0)
					{
						x0 = x;
					}
					if (x > x1)
					{
						x1 = x;
					}
				}
			}
			// Update vertical values for bounding box.
			if (identifier_seen)
			{
				if (y < y0)
				{
					y0 = y;
				}
				if (y > y1)
				{
					y1 = y;
				}
			}
		}
		__syncthreads(); // Need to sync again, so all threads are done with this row before loading the next row.
	}

	// Now, wrap up by inserting calculated values into array.
	if (idx < sp_count)
	{
		int pos = idx * 6;
		d_bbox_array[pos + 1] = x0;
		d_bbox_array[pos + 2] = y0;
		d_bbox_array[pos + 3] = x1;
		d_bbox_array[pos + 4] = y1;
		d_bbox_array[pos + 5] = size;
	}
}

__global__ void find_initial_cell_regions(int width, int height, int xdiv, int ydiv, int buffer, int* pixdata, unsigned char* edge_data)
{
	// CUDA kernel for finding the initial regions within a region of the image corresponding to the min value.
	// width and height are dimensions for the overall image, and map to the size of pixdata nd edge_data.
	// xdiv and ydiv are the size of idividual cells.
	// buffer is the number of pixels on each edge that are not eligible to be seeds.
	// pixdata is the superpixel identifier data for the image, where an integer identifier is assigned to each pixel.
	// edge_data is the gradient data for the image, where edges get higher numbers.

	// This kernel will be called with each block corresponding to one cell.
	int n_cell = blockIdx.x;

	// Shared data across the threads of this cell.
	__shared__ unsigned char min; // Minimum value of edge_data in this region.
	__shared__ unsigned char min_array[threads_per_block]; // One minimum calculated for each thread.
	__shared__ int region_pos; // Position within the region.
	__shared__ int seed_identifier; // Identifier to record in each contiguous area of min.
	__shared__ bool need_flood_fill; // Indicates to threads whether a flood-fill is needed.
	__shared__ bool flood_fill_done; // Flag for completion of flood fill.

	// Calculate the number of cells in each direction.
	int n_y_cells = (height + ydiv / 2) / ydiv;
	int n_x_cells = (width + xdiv / 2) / xdiv;
	//int num_cells = n_x_cells * n_y_cells;

	//	Region calculations:
	int i = n_cell % n_x_cells;
	int j = n_cell / n_x_cells;
	int x_offset = buffer + i * xdiv;
	int y_offset = buffer + j * ydiv;
	int region_width = xdiv - 2 * buffer;
	int region_height = ydiv - 2 * buffer;

	// First task is to calculate the minimum value of edge_data in the region.

	if ((x_offset < width) && (y_offset < height) && (x_offset >= 0) && (y_offset >= 0)) // Do basic dimension testing.
	{
		int N = width * height;
		// Adjust portion dimensions if necessary.
		if (x_offset + region_width > width)
		{
			region_width = width - x_offset;
		}
		if (y_offset + region_height > height)
		{
			region_height = height - y_offset;
		}
		if ((region_width > 0) && (region_height > 0))
		{
			int N_region_elements = region_width * region_height; // Number of elements in the target portion of the array.

			int idx = threadIdx.x;
			min_array[idx] = 255; // Maximum value for unsigned char.
			for (int i = idx; i < N_region_elements; i += threads_per_block) // Striding through portion of the array.
			{
				int portion_x = i % region_width; // x value within window of the target region.
				int portion_y = i / region_width; // y value within window of the target region.
				int x = x_offset + portion_x; // Overall x value.
				int y = y_offset + portion_y; // Overall y value.
				int pos = x + y * width; // Position of the element within the full array.
				if (edge_data[pos] < min_array[idx])
				{
					min_array[idx] = edge_data[pos];
				}
			}
			__syncthreads(); // Wait for all threads to finish.

			// Reduce the thread_min array to one element.  Don't assume that threads_per_block is a power of two (although it almost always is).
			unsigned int working_set = threads_per_block;
			while (working_set > 1)
			{
				int half_ceiling = (working_set + 1) / 2; // Handle case where working_set is odd.
				if ((idx < half_ceiling) && (idx + half_ceiling < working_set))
				{
					if (min_array[idx + half_ceiling] < min_array[idx])
					{
						min_array[idx] = min_array[idx + half_ceiling];
					}
				}
				working_set = half_ceiling;
				__syncthreads(); // Each pass through the reduction process needs to wait for all threads to complete.
			}

			// Put final answer in min.
			if (0 == idx)
			{
				min = min_array[0];
			}
			__syncthreads(); // Now, we have the minimum value of the region in edge_data in min.

			// Now, we need to go through the pixels of the region sequentially.
			if (0 == idx)
			{
				region_pos = 0;
				seed_identifier = 1;
			}
			__syncthreads();
			while (region_pos < N_region_elements)
			{
				int region_x = region_pos % region_width;
				int region_y = region_pos / region_width;
				int full_pos = region_x + x_offset + (region_y + y_offset) * width;
				if (0 == idx)
				{
					if ((0 == pixdata[full_pos]) && (min == edge_data[full_pos]))
					{
						pixdata[full_pos] = seed_identifier;
						need_flood_fill = true;
						flood_fill_done = false;
					}
					else {
						need_flood_fill = false;
					}
				}
				__syncthreads();

				// All threads check to see whether flood fill is needed.
				if (need_flood_fill)
				{
					// Do flood-fill based on min and seed_identifier.

					while (false == flood_fill_done)
					{
						__syncthreads();
						if (0 == threadIdx.x)
						{
							flood_fill_done = true; // Start off each round with an assumption that we are done.
						}
						__syncthreads();

						// In the non-CUDA code, the growth of each contiguous region is not limited by the buffer zones.  Use flood-specific
						// variables here.
						int flood_region_width = xdiv;
						int flood_region_height = ydiv;
						int flood_x_offset = i * xdiv;
						int flood_y_offset = j * ydiv;
						if (flood_x_offset + flood_region_width > width)
						{
							flood_region_width = width - flood_x_offset;
						}
						if (flood_y_offset + flood_region_height > height)
						{
							flood_region_height = height - flood_y_offset;
						}
						int flood_N_region_elements = flood_region_width * flood_region_height;

						for (int flood_region_pos = threadIdx.x; flood_region_pos < flood_N_region_elements; flood_region_pos += threads_per_block) // flood_region_pos is the position within the region.  Stride across the region.
						{
							int flood_thread_x = flood_region_pos % flood_region_width; // Region x value for this thread.
							int flood_thread_y = flood_region_pos / flood_region_width; // Region y value for this thread.
							if ((flood_thread_x + flood_x_offset < width) && (flood_thread_y + flood_y_offset < height)) // Ensure we are within the overall image.
							{
								int flood_full_pos = flood_thread_x + flood_x_offset + (flood_thread_y + flood_y_offset) * width;
								if ((min == edge_data[flood_full_pos]) && (seed_identifier != pixdata[flood_full_pos])) // We are only able to expand if we are on a pixel with value min.
								{
									// Check adjacent pixels.

									// Top
									if (flood_thread_y > 0) // Don't look outside of defined region.
									{
										if (seed_identifier == pixdata[flood_full_pos - width])
										{
											pixdata[flood_full_pos] = seed_identifier;
											if (flood_fill_done)
											{
												flood_fill_done = false; // Overrides earlier assumption that we are done.
											}
										}
									}

									// Bottom
									if (flood_thread_y < flood_region_height - 1) // Don't look outside of defined region.
									{
										if (seed_identifier == pixdata[flood_full_pos + width])
										{
											pixdata[flood_full_pos] = seed_identifier;
											if (flood_fill_done)
											{
												flood_fill_done = false; // Overrides earlier assumption that we are done.
											}
										}
									}

									// Left
									if (flood_thread_x > 0) // Don't look outside of defined region.
									{
										if (seed_identifier == pixdata[flood_full_pos - 1])
										{
											pixdata[flood_full_pos] = seed_identifier;
											if (flood_fill_done)
											{
												flood_fill_done = false; // Overrides earlier assumption that we are done.
											}
										}
									}

									// Right
									if (flood_thread_x < flood_region_width - 1) // Don't look outside of defined region.
									{
										if (seed_identifier == pixdata[flood_full_pos + 1])
										{
											pixdata[flood_full_pos] = seed_identifier;
											if (flood_fill_done)
											{
												flood_fill_done = false; // Overrides earlier assumption that we are done.
											}
										}
									}

								}
							}
						}
						__syncthreads(); // Wait for all threads to complete for this round.
					}
				}

				if (0 == idx)
				{
					if (need_flood_fill)
					{
						seed_identifier = seed_identifier + 1;
					}
					region_pos = region_pos + 1;
				}
				__syncthreads();
			}

		}
	}
}

__global__ void cell_grow(int width, int height, int xdiv, int ydiv, unsigned char* edge_data, int* temp_spdata, int* delta_spdata, bool* collision_data, int* stage_sizes)
{
	// CUDA kernel for growing regions in a cell to find the one that overtakes any others.
	// width and height are dimensions for the overall image, and map to the size of temp_spdata and edge_data.
	// xdiv and ydiv are the size of idividual cells.
	// edge_data is the gradient data for the image, where edges get higher numbers.
	// temp_spdata is the copy of spdata that is used for growing the superpixels.
	// delta_spdata is used for holding a temporary value before it is written to temp_spdata.
	// collision_data is used for tracking superpixel intersections when they happen.
	// stage_sizes holds the size of each superpixel.

	// This kernel will be called with each block corresponding to one cell.
	int n_cell = blockIdx.x;

	// Shared data across the threads of this cell.
	__shared__ bool continue_step;
	__shared__ int numbers_seen[threads_per_block];
	__shared__ unsigned char level; // Current level
	__shared__ bool all_done; // Indicates that there are no longer multiple superpixels in this cell.
	__shared__ int iter_pos; // Position value while iterating through the region.
	__shared__ int contig_identifier[4]; // The four potential identifiers that are contiguous to a location.
	__shared__ int contig_size[4]; // Holds the size of each superpixel included in contig_identifiers.
	__shared__ int target_identifier; // Identifier for absorption step.

	// Calculate the number of cells in each direction.
	int n_x_cells = (width + xdiv / 2) / xdiv;
	if (0 == n_x_cells)
	{
		n_x_cells = 1;
	}

	//	Region calculations (no buffer is used here):
	int cell_x = n_cell % n_x_cells;
	int cell_y = n_cell / n_x_cells;
	int x_offset = cell_x * xdiv;
	int y_offset = cell_y * ydiv;
	int region_width = xdiv;
	int region_height = ydiv;

	if ((x_offset < width) && (y_offset < height) && (x_offset >= 0) && (y_offset >= 0)) // Do basic dimension testing.
	{
		// Adjust portion dimensions if necessary.
		if (x_offset + region_width > width)
		{
			region_width = width - x_offset;
		}
		if (y_offset + region_height > height)
		{
			region_height = height - y_offset;
		}
		if ((region_width > 0) && (region_height > 0))
		{
			int N_region_elements = region_width * region_height; // Number of elements in the target portion of the array.

			int idx = threadIdx.x;
			if (0 == idx)
			{
				continue_step = true;
				level = 0;
				all_done = false;
			}
			__syncthreads();

			while (false == all_done) // Loops through levels from 0 to 255, unless stopped sooner.
			{
				while (continue_step && (false == all_done)) // Loops through growth steps at the current level.
				{
					// First, we clear out and then calculate stage_sizes.
					for (int i = idx; i < N_region_elements; i += threads_per_block) // Striding through portion of the array.
					{
						int portion_x = i % region_width; // x value within window of the target region.
						int portion_y = i / region_width; // y value within window of the target region.
						int x = x_offset + portion_x; // Overall x value.
						int y = y_offset + portion_y; // Overall y value.
						int pos = x + y * width; // Position of the element within the full array.
						stage_sizes[pos] = 0;
					}
					__syncthreads();
					for (int i = idx; i < N_region_elements; i += threads_per_block) // Striding through portion of the array.
					{
						int portion_x = i % region_width; // x value within window of the target region.
						int portion_y = i / region_width; // y value within window of the target region.
						int x = x_offset + portion_x; // Overall x value.
						int y = y_offset + portion_y; // Overall y value.
						int pos = x + y * width; // Position of the element within the full array.
						int value = temp_spdata[pos];
						if (value > 0)
						{
							// value becomes the region-specific location to update.
							portion_x = value % region_width; // x value within window of the target region.
							portion_y = value / region_width; // y value within window of the target region.
							x = x_offset + portion_x; // Overall x value.
							y = y_offset + portion_y; // Overall y value.
							pos = x + y * width; // Position of the element within the full array.
							atomicAdd(&stage_sizes[pos], 1);
						}
					}
					__syncthreads();

					// We need to look for cases where an open pixel would, on the next step, lead to a collision between
					// two or more regions.  These are set to true in collision_data.
					for (int i = idx; i < N_region_elements; i += threads_per_block) // Striding through portion of the array.
					{
						int portion_x = i % region_width; // x value within window of the target region.
						int portion_y = i / region_width; // y value within window of the target region.
						int x = x_offset + portion_x; // Overall x value.
						int y = y_offset + portion_y; // Overall y value.
						int pos = x + y * width; // Position of the element within the full array.
						collision_data[pos] = false;
						if ((0 == temp_spdata[pos]) && (level >= edge_data[pos])) // Open position at the right level.
						{
							int value = 0;
							// Check adjacent pixels.

							// Top
							if (portion_y > 0) // Don't look outside of defined region.
							{
								value = temp_spdata[pos - width]; // No need to test, since this is first.
							}

							// Bottom
							if (portion_y < region_height - 1) // Don't look outside of defined region.
							{
								if (value > 0) // At least one has already been seen.
								{
									if ((temp_spdata[pos + width] > 0) && (temp_spdata[pos + width] != value))
									{
										collision_data[pos] = true;
									}
								}
								else {
									value = temp_spdata[pos + width];
								}
							}

							// Left
							if (portion_x > 0) // Don't look outside of defined region.
							{
								if (value > 0) // At least one has already been seen.
								{
									if ((temp_spdata[pos - 1] > 0) && (temp_spdata[pos - 1] != value))
									{
										collision_data[pos] = true;
									}
								}
								else {
									value = temp_spdata[pos - 1];
								}
							}

							// Right
							if (portion_x < region_width - 1) // Don't look outside of defined region.
							{
								if (value > 0) // At least one has already been seen.
								{
									if ((temp_spdata[pos + 1] > 0) && (temp_spdata[pos + 1] != value))
									{
										collision_data[pos] = true;
									}
								} // No need for the else here, since there is at most one adjacent superpixel value.
							}
						}
					}

					// We also need to check for cases where two superpixels have already grown into each other.
					for (int i = idx; i < N_region_elements; i += threads_per_block) // Striding through portion of the array.
					{
						int portion_x = i % region_width; // x value within window of the target region.
						int portion_y = i / region_width; // y value within window of the target region.
						int x = x_offset + portion_x; // Overall x value.
						int y = y_offset + portion_y; // Overall y value.
						int pos = x + y * width; // Position of the element within the full array.
						if ((false == collision_data[pos]) && (temp_spdata[pos] > 0))
						{
							int value = 0;
							// Check adjacent pixels.

							// Top
							if (portion_y > 0) // Don't look outside of defined region.
							{
								value = temp_spdata[pos - width];
								if ((value > 0) && (value != temp_spdata[pos]))
								{
									collision_data[pos] = true;
								}
							}

							// Bottom
							if (portion_y < region_height - 1) // Don't look outside of defined region.
							{
								value = temp_spdata[pos + width];
								if ((value > 0) && (value != temp_spdata[pos]))
								{
									collision_data[pos] = true;
								}
							}

							// Left
							if (portion_x > 0) // Don't look outside of defined region.
							{
								value = temp_spdata[pos - 1];
								if ((value > 0) && (value != temp_spdata[pos]))
								{
									collision_data[pos] = true;
								}
							}

							// Right
							if (portion_x < region_width - 1) // Don't look outside of defined region.
							{
								value = temp_spdata[pos + 1];
								if ((value > 0) && (value != temp_spdata[pos]))
								{
									collision_data[pos] = true;
								}
							}

						}
					}


					// Now, the collision_data matrix defines all upcoming intersections between superpixels
					// that need to be addressed.

					// We need to go through and address each upcoming intersection by switching the values of
					// each superpixel to that of the larger of the ones that are part of the intersection (with
					// a tie going to the one with a smaller identifier).

					// We'll use a single thread to track the progress.
					if (0 == idx)
					{
						iter_pos = 0;
					}
					__syncthreads();

					while (iter_pos < N_region_elements) // Step through all elements in the region.
					{

						// Find the identifiers contiguous to the current location defined by iter_pos.
						if (0 == idx)
						{
							int portion_x = iter_pos % region_width; // x value within window of the target region.
							int portion_y = iter_pos / region_width; // y value within window of the target region.
							int x = x_offset + portion_x; // Overall x value.
							int y = y_offset + portion_y; // Overall y value.
							int pos = x + y * width; // Position of the element within the full array.
							int value;

							contig_identifier[0] = 0;
							contig_identifier[1] = 0;
							contig_identifier[2] = 0;
							contig_identifier[3] = 0;
							contig_size[0] = 0;
							contig_size[1] = 0;
							contig_size[2] = 0;
							contig_size[3] = 0;
							if (collision_data[pos]) // Find identifiers for contiguous pixels.
							{
								// Top
								if (portion_y > 0) // Don't look outside of defined region.
								{
									value = temp_spdata[pos - width];
									if (value > 0)
									{
										contig_identifier[0] = value; // No need to test, since this is first.
										int size_pos = x_offset + value % region_width + width * (y_offset + value / region_width);
										contig_size[0] = stage_sizes[size_pos];
									}
								}

								// Bottom
								if (portion_y < (region_height - 1))
								{
									value = temp_spdata[pos + width];
									if (value > 0)
									{
										int index = 0;
										while ((index < 4) && (contig_identifier[index] > 0)) // Position of index is already filled.
										{
											if (value == contig_identifier[index]) // This is a repeat.
											{
												index = 4; // Ends while loop.
											}
											index++;
										}
										if (index < 4)
										{
											contig_identifier[index] = value;
											int size_pos = x_offset + value % region_width + width * (y_offset + value / region_width);
											contig_size[index] = stage_sizes[size_pos];
										}
									}
								}

								// Left
								if (portion_x > 0) // Don't look outside of defined region.
								{
									value = temp_spdata[pos - 1];
									if (value > 0)
									{
										int index = 0;
										while ((index < 4) && (contig_identifier[index] > 0)) // Position of index is already filled.
										{
											if (value == contig_identifier[index]) // This is a repeat.
											{
												index = 4; // Ends while loop.
											}
											index++;
										}
										if (index < 4)
										{
											contig_identifier[index] = value;
											int size_pos = x_offset + value % region_width + width * (y_offset + value / region_width);
											contig_size[index] = stage_sizes[size_pos];
										}
									}
								}

								// Right
								if (portion_x < region_width - 1) // Don't look outside of defined region.
								{
									value = temp_spdata[pos + 1];
									if (value > 0)
									{
										int index = 0;
										while ((index < 4) && (contig_identifier[index] > 0)) // Position of index is already filled.
										{
											if (value == contig_identifier[index]) // This is a repeat.
											{
												index = 4; // Ends while loop.
											}
											index++;
										}
										if (index < 4)
										{
											contig_identifier[index] = value;
											int size_pos = x_offset + value % region_width + width * (y_offset + value / region_width);
											contig_size[index] = stage_sizes[size_pos];
										}
									}
								}
							}
						}
						__syncthreads();

						// Now, there should be two or more identifiers in contig_identifier if collision_data for this location is true.
						// It may not be, though, because of previous absorb steps.

						if (contig_identifier[1] > 0)
						{
							int index;
							// Now, contig_identifier contains multiple identifiers that need to be combined.
							// contig_size contains the size of each superpixel corresponding to those identifiers.

							// Absorb step
							// First, find the identifier of the largest superpixel.
							if (0 == idx)
							{
								int target_index = 0;
								index = 1;
								while ((index < 4) && (contig_size[index] > 0))
								{
									if (contig_size[target_index] < contig_size[index])
									{
										target_index = index;
									}
									else if (contig_size[target_index] == contig_size[index]) // Special handling for ties.
									{
										if (contig_identifier[index] < contig_identifier[target_index])
										{
											target_index = index;
										}
									}
									index++;
								}
								target_identifier = contig_identifier[target_index];
							}
							__syncthreads();

							// The value of target_identifier is the value that all identifiers that show up in contig_identifier should be changed to.
							for (int i = idx; i < N_region_elements; i += threads_per_block) // Striding through portion of the array.
							{
								// Find the position for this element.
								int local_portion_x = i % region_width; // x value within window of the target region.
								int local_portion_y = i / region_width; // y value within window of the target region.
								int local_x = x_offset + local_portion_x; // Overall x value.
								int local_y = y_offset + local_portion_y; // Overall y value.
								int local_pos = local_x + local_y * width; // Position of the element within the full array.
								for (index = 0; index < 4; ++index)
								{
									int local_identifier = contig_identifier[index];
									if ((local_identifier > 0) && (local_identifier != target_identifier))
									{
										if (temp_spdata[local_pos] == local_identifier)
										{
											temp_spdata[local_pos] = target_identifier; // Because only this thread can write to this location, no conflict.
										}
									}
								}
							}
							__syncthreads();

							// Update the sizes
							if (0 == idx)
							{
								int number_shifted = 0;
								index = 0;
								while (index < 4)
								{
									if ((contig_identifier[index] > 0) && (contig_identifier[index] != target_identifier))
									{
										number_shifted += contig_size[index];
										int size_pos = x_offset + contig_identifier[index] % region_width + width * (y_offset + contig_identifier[index] / region_width);
										stage_sizes[size_pos] = 0;
									}
									index++;
								}
								int size_pos = x_offset + target_identifier % region_width + width * (y_offset + target_identifier / region_width);
								stage_sizes[size_pos] += number_shifted;
							}
							__syncthreads();
						}


						if (0 == idx)
						{
							iter_pos++;
						}
						__syncthreads();
					}

					// All collision points have been addressed.
					// Now, grow each superpixel by one pixel.

					if (0 == idx)
					{
						continue_step = false; // Only sets to true if there is growth.
					}
					__syncthreads();

					// First, zero out the delta_spdata matrix.
					for (int i = idx; i < N_region_elements; i += threads_per_block) // Striding through portion of the array.
					{
						int portion_x = i % region_width; // x value within window of the target region.
						int portion_y = i / region_width; // y value within window of the target region.
						int x = x_offset + portion_x; // Overall x value.
						int y = y_offset + portion_y; // Overall y value.
						int pos = x + y * width; // Position of the element within the full array.
						delta_spdata[pos] = 0;
					}
					__syncthreads();

					for (int i = idx; i < N_region_elements; i += threads_per_block) // Striding through portion of the array.
					{
						int portion_x = i % region_width; // x value within window of the target region.
						int portion_y = i / region_width; // y value within window of the target region.
						int x = x_offset + portion_x; // Overall x value.
						int y = y_offset + portion_y; // Overall y value.
						int pos = x + y * width; // Position of the element within the full array.
						int value = 0;

						if ((0 == temp_spdata[pos]) && (edge_data[pos] <= level))
						{
							// Top
							if (portion_y > 0) // Don't look outside of defined region.
							{
								value = temp_spdata[pos - width];
								if (value > 0)
								{
									delta_spdata[pos] = value; // Write to delta matrix, for later consolidation with spdata.
								}
							}

							// Bottom
							if ((value == 0) && (portion_y < region_height - 1)) // Don't look outside of defined region.
							{
								value = temp_spdata[pos + width];
								if (value > 0)
								{
									delta_spdata[pos] = value;
								}
							}

							// Left
							if ((value == 0) && (portion_x > 0)) // Don't look outside of defined region.
							{
								value = temp_spdata[pos - 1];
								if (value > 0)
								{
									delta_spdata[pos] = value;
								}
							}

							// Right
							if ((value == 0) && (portion_x < region_width - 1)) // Don't look outside of defined region.
							{
								value = temp_spdata[pos + 1];
								if (value > 0)
								{
									delta_spdata[pos] = value;
								}
							}
						}
					}
					__syncthreads();

					// Now that the growth step is done, update with the delta matrix.
					for (int i = idx; i < N_region_elements; i += threads_per_block) // Striding through portion of the array.
					{
						int portion_x = i % region_width; // x value within window of the target region.
						int portion_y = i / region_width; // y value within window of the target region.
						int x = x_offset + portion_x; // Overall x value.
						int y = y_offset + portion_y; // Overall y value.
						int pos = x + y * width; // Position of the element within the full array.
						int value = delta_spdata[pos];
						if (value > 0) // This value was updated in the growth step above.
						{
							temp_spdata[pos] = value;
							if (false == continue_step)
							{
								continue_step = true;
							}
						}
					}
					__syncthreads();

					// Last part of main loop:
					// Count unique region numbers.  If more than one, need to continue.
					if (0 == idx)
					{
						all_done = true;
					}
					__syncthreads();

					numbers_seen[idx] = 0; // Initialize value.
					for (int i = idx; i < N_region_elements; i += threads_per_block) // Striding through portion of the array.
					{
						int portion_x = i % region_width; // x value within window of the target region.
						int portion_y = i / region_width; // y value within window of the target region.
						int x = x_offset + portion_x; // Overall x value.
						int y = y_offset + portion_y; // Overall y value.
						int pos = x + y * width; // Position of the element within the full array.
						int value = temp_spdata[pos]; // Value of this element.
						if (value > 0)
						{
							if (0 == numbers_seen[idx])
							{
								numbers_seen[idx] = value;
							}
							else {
								if (value != numbers_seen[idx])
								{
									all_done = false; // There is more than one numbered region.  This thread happened to hit two.
								}
							}
						}
					}
					__syncthreads(); // Wait for all threads to finish.

					// Use single thread to check on numbers seen.
					if (all_done && (0 == idx)) // all_done being true here may not be final, but if it is false, then no need for this loop.
					{
						int value = 0;
						for (int i = 0; all_done && (i < threads_per_block); ++i)
						{
							if (numbers_seen[i] > 0)
							{
								if (0 == value)
								{
									value = numbers_seen[i];
								}
								else {
									if (value != numbers_seen[i]) // Since value has a different number, there are at least two superpixels.
									{
										all_done = false;
									}
								}
							}
						}
					}
					__syncthreads(); // Now, all_done tells us whether to keep going.
				} // continue_step tells us whether to do another growth step at this level.
				// Increase level.
				__syncthreads();
				if (0 == idx)
				{
					if (level < 255)
					{
						level += 1;
						continue_step = true;
					}
					else {
						all_done = true; // Should not get here with all_done set to false, but we're done anyway.
					}
				}
				__syncthreads();
			}
		}
	}
}

__global__ void resolve_seeds(int width, int height, int xdiv, int ydiv, int buffer, int* spdata, int* temp_spdata, int* d_seeds)
{
	// CUDA kernel for determining the seed for each cell.
	// width and height are dimensions for the overall image, and map to the size of spdata and temp_spdata.
	// xdiv and ydiv are the size of idividual cells.
	// buffer is the space on the outside of the cell where a seed cannot be located.
	// temp_spdata is the copy of spdata that was used for growing the superpixels.
	// d_seeds is the array into which seed results are placed.  They are ordered by cell number, with x followed by y.

	// This kernel will be called with each block corresponding to one cell.
	int n_cell = blockIdx.x;

	// Calculate the number of cells in each direction.
	int n_x_cells = (width + xdiv / 2) / xdiv;
	if (0 == n_x_cells)
	{
		n_x_cells = 1;
	}

	__shared__ int pos_array[threads_per_block]; // One position calculated for each thread.
	__shared__ int identifier;

	//	Region calculations (no buffer is used here):
	int cell_x = n_cell % n_x_cells;
	int cell_y = n_cell / n_x_cells;
	int x_offset = cell_x * xdiv;
	int y_offset = cell_y * ydiv;
	int region_width = xdiv;
	int region_height = ydiv;

	int idx = threadIdx.x;
	if (0 == idx)
	{
		d_seeds[2 * n_cell] = -1; // Default values in case this kernel is called with bad parameters.
		d_seeds[2 * n_cell + 1] = -1;
	}
	__syncthreads();

	if ((x_offset < width) && (y_offset < height) && (x_offset >= 0) && (y_offset >= 0)) // Do basic dimension testing.
	{
		// Adjust portion dimensions if necessary.
		if (x_offset + region_width > width)
		{
			region_width = width - x_offset;
		}
		if (y_offset + region_height > height)
		{
			region_height = height - y_offset;
		}
		if ((region_width > 0) && (region_height > 0))
		{
			int N_region_elements = region_width * region_height; // Number of elements in the target portion of the array.

			if (0 == idx)
			{
				// Find the identifier for the remaining superpixel in this cell.
				int value = 0;
				for (int i = 0; (i < N_region_elements) && (value == 0); ++i) // Striding through portion of the array, stop when non-zero found.
				{
					int portion_x = i % region_width; // x value within window of the target region.
					int portion_y = i / region_width; // y value within window of the target region.
					int x = x_offset + portion_x; // Overall x value.
					int y = y_offset + portion_y; // Overall y value.
					int pos = x + y * width; // Position of the element within the full array.
					value = temp_spdata[pos];
				}
				if (value > 0)
				{
					identifier = value;
				}
				else {
					identifier = -1; // No identifier found.
				}
			}
			__syncthreads();

			if (identifier > 0)
			{
				// Find region position with lowest number that maps to an identifier in spdata.
				int first_pos_found = -1; // Initialize to invalid value.
				for (int local_pos = idx; (local_pos < N_region_elements) && (first_pos_found < 0); local_pos += threads_per_block) // Striding through portion of the array.
				{
					int portion_x = local_pos % region_width; // x value within window of the target region.
					int portion_y = local_pos / region_width; // y value within window of the target region.
					int x = x_offset + portion_x; // Overall x value.
					int y = y_offset + portion_y; // Overall y value.
					int pos = x + y * width; // Position of the element within the full array.
					if ((portion_x >= buffer) && (portion_y >= buffer) && (portion_x < region_width - buffer) && (portion_y < region_height - buffer))
					{
						int value = spdata[pos];
						if (identifier == value)
						{
							first_pos_found = pos; // Assign value of position in overall image.
						}
					}
				}
				pos_array[idx] = first_pos_found;
				__syncthreads(); // Wait for all threads to finish.			

				// Now, reduce the answer array to one element.  Don't assume that threads_per_block is a power of two (although it almost always is).
				unsigned int working_set = threads_per_block;
				while (working_set > 1)
				{
					int half_ceiling = (working_set + 1) / 2; // Handle case where working_set is odd.
					if ((idx < half_ceiling) && (idx + half_ceiling < working_set))
					{
						if (pos_array[idx + half_ceiling] >= 0)
						{
							if ((pos_array[idx] < 0) || (pos_array[idx + half_ceiling] < pos_array[idx]))
							{
								pos_array[idx] = pos_array[idx + half_ceiling];
							}
						}
					}
					working_set = half_ceiling;
					__syncthreads(); // Each pass through the reduction process needs to wait for all threads to complete.
				}

				if (0 == idx)
				{
					if (pos_array[0] >= 0) // Found a valid seed.
					{
						// Seeds need to be returned in reference to the overall image.
						int x = pos_array[0] % width;
						int y = pos_array[0] / width; 
						d_seeds[2 * n_cell] = x;
						d_seeds[2 * n_cell + 1] = y;
					}
				}
				__syncthreads();
			}
		}
	}
}

bool c_watershed_grow(int x, int y, SuperPixel* list_head)
{
	// CUDA version of watershed growth algorithm.
	// x and y define the size of the image.
	// list_head is a pointer to the first SuperPixel.
	// 
	// device_pixel_data is the on-device version of the superpixel data, where each location is a SuperPixel identifier.
	// device_gradient_data is the on-device version of the gradient image, where each location is between 0 (no gradient) to 255 (maximum gradient).

	bool ret = true;
	bool* c_result = BoolArray(1, 1, false);
	cudaError_t cudaStatus;
	int blocksPerGrid = (x * y + threads_per_block - 1) / threads_per_block;

	// First, we need to set all identifiers on the pixel_data to zero.
	SPixelData* host_data = list_head->GetPixelData();
	ret = host_data->Reset();
	if (false == ret)
	{
		throw std::runtime_error("Failed to reset SPixelData in c_watershed_grow.\n");
		return ret;
	}

	// Now, we need to set each SuperPixel seed to its identifier value.
	SuperPixel* current = list_head->GetHead();
	while (ret && (NULL != current))
	{
		PointPair seed_location = current->GetSeed();
		int identifier = current->GetIdentifier();
		ret = host_data->SetPixel(seed_location.x, seed_location.y, identifier);
		current = current->GetNext();
	}

	if (ret)
	{
		// Place the host superpixel data on the CUDA device.
		ret = CopyFromHost(host_data->GetData(), x * y, host_data->GetDeviceData());
		if (false == ret)
		{
			throw std::runtime_error("Failed to copy superpixel data to CUDA device.\n");
			return ret;
		}

		bool* ResultsArray = BoolArray(blocksPerGrid, 1, false); // Set up the results array, so each block can flag whether more processing is needed.
		bool continue_with_level = true;
		for (int level = 0; level < 256; level++) // The levels for the watershed algorithm.
		{
			std::cout << ".";
			continue_with_level = true;
			int count = 0;
			while (continue_with_level)
			{
				count++;
				ResetBoolArray(ResultsArray, blocksPerGrid, 1, false);
				watershed_step << < blocksPerGrid, threads_per_block >> > (level, x, y, host_data->GetDeviceData(), list_head->GetGradient()->GetCData(), ResultsArray);
				cudaDeviceSynchronize();
				cudaStatus = cudaGetLastError();
				if (cudaStatus != cudaSuccess)
				{
					throw std::runtime_error("Failed in watershed_step.\n");
					return false;
				}
				// Now, change the new identifier values from negative to positive.  They are initially stored as negatives, to distinguish them from the pre-state.
				update_values << < blocksPerGrid, threads_per_block >> > (x, y, host_data->GetDeviceData(), ResultsArray);
				cudaDeviceSynchronize();
				cudaStatus = cudaGetLastError();
				if (cudaStatus != cudaSuccess)
				{
					throw std::runtime_error("Failed in update_values in c_watershed_grow.\n");
					return false;
				}
				continue_with_level = TestResult(blocksPerGrid, 1, ResultsArray, c_result);  // If any expansion of superpixels in the last run, need to do another at this level.
			}
		}

		// Copy the CUDA device superpixel data to the host.
		ret = CopyToHost(host_data->GetDeviceData(), x * y, host_data->GetData());
		if (false == ret)
		{
			throw std::runtime_error("Failed to copy superpixel data to host.\n");
			return ret;
		}

		// Need to recompute the bounding box, size, and last_level for each SuperPixel.		
		ret = c_find_bbox(x, y, list_head);
		if (false == ret)
		{
			throw std::runtime_error("Failed to calculate bounding boxes in c_watershed_grow.\n");
			return ret;
		}

		// Release allocated memory.
		ret = FreeBoolArray(c_result);
		if (false == ret)
		{
			throw std::runtime_error("Error freeing placeholder boolean value on CUDA device.\n");
			return ret;
		}
		ret = FreeBoolArray(ResultsArray);
		if (false == ret)
		{
			throw std::runtime_error("Error freeing ResultsArray on CUDA device.\n");
			return ret;
		}
	}
	return ret;
}

bool c_find_bbox(int x, int y, SuperPixel* list_head)
{
	// CUDA function to populate bounding boxes in superpixels based on pixelmap.
	// x and y define size of image.
	// list_head is a pointer to the first superpixel.

	bool ret = true;
	cudaError_t cudaStatus;

	// Get pixeldata from the SuperPixels.
	SPixelData* sp_data = list_head->GetPixelData();

	// We need the count for the number of superpixels.
	SuperPixel* current = list_head->GetHead();
	int sp_count = 0; // The number of superpixels.
	while (NULL != current)
	{
		sp_count++;
		current = current->GetNext();
	}

	// Figure out how many grids we will need, so each thread can be mapped to a single superpixel.
	int blocksPerGrid = (sp_count + threads_per_block - 1) / threads_per_block;

	// Allocate memory for on-host and on-device arrays to hold the bounding box and size information for
	// each superpixel.
	// Since each bounding box takes four integers, and the size takes 1, there would be five integers per
	// superpixel.  However, in some cases, the sequence of superpixel identifiers may not increase monotonically,
	// and there could be gaps.  So, we need to also include the actual identifiers in the arrays.
	// The layout is going to be: identifier, x0, y0, x1, y1, size.

	int* h_bbox_array = (int*)malloc(6 * sizeof(int) * sp_count);
	if (NULL == h_bbox_array)
	{
		throw std::runtime_error("Failed to allocate memory for h_bbox_array in c_find_bbox.\n");
		return false;
	}

	int* d_bbox_array = IntArray(6 * sp_count, 1, false);
	if (NULL == d_bbox_array)
	{
		throw std::runtime_error("Failed to allocate memory for d_bbox_array in c_find_bbox.\n");
		return false;
	}

	// Now, set values for the identifiers.
	current = list_head->GetHead();
	sp_count = 0;
	while (NULL != current)
	{
		h_bbox_array[sp_count * 6] = current->GetIdentifier();
		sp_count++;
		current = current->GetNext();
	}

	// Copy the array to the device.
	cudaStatus = cudaMemcpy(d_bbox_array, h_bbox_array, 6 * sizeof(int) * sp_count, cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed to copy memory to device in c_find_bbox.\n");
		return false;
	}

	// Begin calculations.
	// *** Test to confirm that x*sizeof(int) is not too large for shared memory.  If it is, then break this into strips.

	calc_bbox << < blocksPerGrid, threads_per_block, x * sizeof(int) >> > (x, y, sp_data->GetDeviceData(), d_bbox_array, sp_count);
	cudaDeviceSynchronize();
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed in calc_bbox.\n");
		return false;
	}

	// Copy the array back to the host.
	cudaStatus = cudaMemcpy(h_bbox_array, d_bbox_array, 6 * sizeof(int) * sp_count, cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed to copy memory to host in c_find_bbox.\n");
		return false;
	}

	// Put calculated data into superpixels.
	current = list_head->GetHead();
	sp_count = 0;
	while (NULL != current)
	{
		int pos = sp_count * 6;
		RectQuad bbox;
		bbox.x0 = h_bbox_array[pos + 1];
		bbox.y0 = h_bbox_array[pos + 2];
		bbox.x1 = h_bbox_array[pos + 3];
		bbox.y1 = h_bbox_array[pos + 4];
		current->SetWindow(bbox);
		current->SetSize(h_bbox_array[pos + 5]);
		current->FindEdgePixels();
		current->SetLevelComplete(255);
		sp_count++;
		current = current->GetNext();
	}

	FreeIntArray(d_bbox_array);
	free(h_bbox_array);
	return ret;
}

std::vector <PointPair> c_find_seeds(int width, int height, GradData* edge, SPixelData* pixeldata, SuperPixel* list_head, int xdiv, int ydiv, int buffer)
{
	std::vector <PointPair> ret;
	ret.clear();

	cudaError_t cudaStatus;
	int count = 0;
	pixeldata->Reset();
	unsigned char* edge_data = edge->GetCData();
	int* spdata = pixeldata->GetDeviceData();
	// Clear spdata, since we will use this for growing seeds in each cell to determine the best seed.
	ResetIntArray(spdata, width, height, 0);

	if ((buffer >= xdiv / 2) || (buffer >= ydiv / 2))
	{
		buffer = 0;
	}

	// Calculate the number of cells.  This will be the number of blocks used.
	int n_y_cells = (height + ydiv / 2) / ydiv;
	int n_x_cells = (width + xdiv / 2) / xdiv;
	int num_cells = n_x_cells * n_y_cells;

	find_initial_cell_regions << <num_cells, threads_per_block >> > (width, height, xdiv, ydiv, buffer, spdata, edge_data);
	cudaDeviceSynchronize();
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed in find_initial_cell_regions.\n");
		return ret;
	}

	// Set up buffer arrays for the process of growing the superpixels, the sizes of the superpixels, and to find the intersection locations.
	int* temp_spdata = IntArray(width, height, false);
	if (NULL == temp_spdata)
	{
		throw std::runtime_error("Unable to allocate temp_spdata on device in c_find_seeds.\n");
		return ret;
	}
	int* delta_spdata = IntArray(width, height, false);
	if (NULL == delta_spdata)
	{
		throw std::runtime_error("Unable to allocate delta_spdata on device in c_find_seeds.\n");
		return ret;
	}
	int* stage_sizes = IntArray(width, height, false);
	if (NULL == stage_sizes)
	{
		throw std::runtime_error("Unable to allocate stage_sizes on device in c_find_seeds.\n");
		return ret;
	}
	bool* collision_data = BoolArray(width, height, false);
	if (NULL == collision_data)
	{
		throw std::runtime_error("Unable to allocate collision_data on device in c_find_seeds.\n");
		return ret;
	}
	cudaStatus = cudaMemcpy(temp_spdata, spdata, width * height * sizeof(int), cudaMemcpyDeviceToDevice);
	if (cudaStatus != cudaSuccess)
	{
		std::cout << "Error copying memory in c_find_seeds.\n";
		return ret;
	}

	cell_grow << <num_cells, threads_per_block >> > (width, height, xdiv, ydiv, edge_data, temp_spdata, delta_spdata, collision_data, stage_sizes);
	cudaDeviceSynchronize();
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed in cell_grow.\n");
		return ret;
	}

	FreeIntArray(stage_sizes);

	// Now we need to find the identifier from temp_spdata, and use it to get the earliest match of that identifier from
	// spdata, using a row-first process.  That is then the seed location for the cell.

	int* d_seeds = IntArray(2 * num_cells, 1, false); // Array to hold seed values on CUDA device.
	if (NULL == d_seeds)
	{
		throw std::runtime_error("Unable to allocate d_seeds on device in c_find_seeds.\n");
		return ret;
	}
	int* h_seeds = (int*)malloc(2 * num_cells * sizeof(int));// Array to hold seed values on host.
	if (NULL == h_seeds)
	{
		throw std::runtime_error("Unable to allocate h_seeds on host in c_find_seeds.\n");
		return ret;
	}

	resolve_seeds << <num_cells, threads_per_block >> > (width, height, xdiv, ydiv, buffer, spdata, temp_spdata, d_seeds);
	cudaDeviceSynchronize();
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed in resolve_seeds.\n");
		return ret;
	}
	// Copy the seed data back to the host.
	cudaStatus = cudaMemcpy(h_seeds, d_seeds, sizeof(int) * 2 * num_cells, cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed to copy memory to host in c_find_bbox.\n");
		return ret;
	}

	PointPair seed;
	for (int i = 0; i < num_cells; ++i)
	{
		if ((h_seeds[i * 2] >= 0) && (h_seeds[i * 2 + 1] >= 0))
		{
			seed.x = h_seeds[i * 2];
			seed.y = h_seeds[i * 2 + 1];
			ret.push_back(seed);
		}
		else {
			seed.x = -1;
			seed.y = -1;
		}
	}

	free(h_seeds);
	FreeIntArray(d_seeds);
	FreeIntArray(delta_spdata);
	FreeIntArray(temp_spdata);
	FreeBoolArray(collision_data);
	return ret;
}