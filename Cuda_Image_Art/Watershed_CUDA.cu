// Copyright (c) 2023-2026 Steve Young
// Licensed under the MIT License

# include "Watershed_CUDA.cuh"
# include "Utils_CUDA.cuh"

__global__ void watershed_step_2(int width, int height, int* pixel_data, bool* ResultsArray)
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

__global__ void watershed_step_1(int level, int width, int height, int* pixel_data, unsigned char* gradient_data, bool* ResultsArray)
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
					if (!ResultsArray[blockIdx.x])
					{
						ResultsArray[blockIdx.x] = true; // Having expanded a superpixel, another pass is needed at this level.
					}
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

__device__ void region_flood_fill(int width, int height, int x_offset, int y_offset, int region_width, int region_height, unsigned char min, int seed_identifier, int* pixdata, unsigned char* edge_data, int* mask, int mask_identifier)
{
	// This is a kernel to do flood-fill in a defined region of pixdata based on min, edge_data, and a seed_identifier.

	__shared__ bool fill_done; // Flag for completion of flood fill.
	int idx = threadIdx.x;
	bool ignore_mask = (NULL == mask) || (mask_identifier < 1);

	if (0 == idx)
	{
		fill_done = false;
	}
	__syncthreads();

	while (false == fill_done)
	{
		if (0 == idx)
		{
			fill_done = true; // Start off each round with an assumption that we are done.
		}
		__syncthreads();

		int N_region_elements = region_width * region_height;

		for (int region_pos = idx; region_pos < N_region_elements; region_pos += threads_per_block) // region_pos is the position within the region.  Stride across the region.
		{
			int thread_x = region_pos % region_width; // Region x value for this thread.
			int thread_y = region_pos / region_width; // Region y value for this thread.
			if ((thread_x + x_offset < width) && (thread_y + y_offset < height)) // Ensure we are within the overall image.
			{
				int full_pos = thread_x + x_offset + (thread_y + y_offset) * width;
				if ((min == edge_data[full_pos]) && (0 == pixdata[full_pos]) && (ignore_mask || (mask_identifier == mask[full_pos]))) // We are only able to expand if we are on a pixel with gradient at min and no identifier.
				{
					// Check adjacent pixels.

					// Top
					if (thread_y > 0) // Don't look outside of defined region.
					{
						if ((seed_identifier == pixdata[full_pos - width]) && (ignore_mask || (mask_identifier == mask[full_pos - width])))
						{
							pixdata[full_pos] = seed_identifier;
							if (fill_done)
							{
								fill_done = false; // Overrides earlier assumption that we are done.
							}
						}
					}

					// Bottom
					if (thread_y < region_height - 1) // Don't look outside of defined region.
					{
						if ((seed_identifier == pixdata[full_pos + width]) && (ignore_mask || (mask_identifier == mask[full_pos + width])))
						{
							pixdata[full_pos] = seed_identifier;
							if (fill_done)
							{
								fill_done = false; // Overrides earlier assumption that we are done.
							}
						}
					}

					// Left
					if (thread_x > 0) // Don't look outside of defined region.
					{
						if ((seed_identifier == pixdata[full_pos - 1]) && (ignore_mask || (mask_identifier == mask[full_pos - 1])))
						{
							pixdata[full_pos] = seed_identifier;
							if (fill_done)
							{
								fill_done = false; // Overrides earlier assumption that we are done.
							}
						}
					}

					// Right
					if (thread_x < region_width - 1) // Don't look outside of defined region.
					{
						if ((seed_identifier == pixdata[full_pos + 1]) && (ignore_mask || (mask_identifier == mask[full_pos + 1])))
						{
							pixdata[full_pos] = seed_identifier;
							if (fill_done)
							{
								fill_done = false; // Overrides earlier assumption that we are done.
							}
						}
					}

				}
			}
		}
		__syncthreads(); // Wait for all threads to complete for this round.
	}
}

__global__ void find_initial_cell_regions(int width, int height, int xdiv, int ydiv, int buffer, int* pixdata, unsigned char* edge_data)
{
	// CUDA kernel for finding the initial regions within a cell of the image corresponding to the min value.
	// width and height are dimensions for the overall image, and map to the size of pixdata nd edge_data.
	// xdiv and ydiv are the size of idividual cells.
	// buffer is the number of pixels on each edge that are not eligible to be seeds.
	// pixdata is the superpixel identifier data for the image, where an integer identifier is assigned to each pixel.
	// edge_data is the gradient data for the image, where edges get higher numbers.
	// mask is a pointer to a matrix of superpixel identifiers.
	// mask_identifier is the identifier in the mask that defines where this kernel operates.

	// This kernel will be called with each block corresponding to one cell.
	int n_cell = blockIdx.x; // Indicates which cell we are working in.
	int idx = threadIdx.x;

	// Shared data across the threads of this cell.
	__shared__ unsigned char min; // Minimum value of edge_data in this region.
	__shared__ unsigned char min_array[threads_per_block]; // One minimum calculated for each thread.
	__shared__ int cell_pos; // Position within the cell.
	__shared__ int seed_identifier; // Identifier to record in each contiguous area of min.
	__shared__ bool need_flood_fill; // Indicates to threads whether a flood-fill is needed.

	// Calculate the number of cells in each direction.
	int n_y_cells = (height + ydiv / 2) / ydiv;
	int n_x_cells = (width + xdiv / 2) / xdiv;

	//	Region calculations:
	int cell_i = n_cell % n_x_cells; // Column number for cell.
	int cell_j = n_cell / n_x_cells; // Row number for cell.
	int x_offset = buffer + cell_i * xdiv;
	int y_offset = buffer + cell_j * ydiv;
	int cell_width = xdiv - 2 * buffer;
	int cell_height = ydiv - 2 * buffer;

	if ((x_offset < width) && (y_offset < height) && (x_offset >= 0) && (y_offset >= 0)) // Do basic dimension testing.
	{
		int N = width * height;
		// Adjust portion dimensions if necessary.
		if (x_offset + cell_width > width)
		{
			cell_width = width - x_offset;
		}
		if (y_offset + cell_height > height)
		{
			cell_height = height - y_offset;
		}
		if ((cell_width > 0) && (cell_height > 0))
		{
			int N_cell_elements = cell_width * cell_height; // Number of elements in the cell.

			// First task is to calculate the minimum value of edge_data in the cell.
			min_uChar_Array_portion_main(edge_data, width, x_offset, y_offset, cell_width, cell_height, min_array, NULL, 0);

			// Put final answer in min.
			if (0 == idx)
			{
				min = min_array[0];
			}
			__syncthreads(); // Now, we have the minimum value of the region in edge_data in min.

			// Now, we need to go through the pixels of the region sequentially.
			if (0 == idx)
			{
				cell_pos = 0;
				seed_identifier = 1;
			}
			__syncthreads();
			while (cell_pos < N_cell_elements)
			{
				int cell_x = cell_pos % cell_width;
				int cell_y = cell_pos / cell_width;
				int full_pos = cell_x + x_offset + (cell_y + y_offset) * width;
				if (0 == idx)
				{
					if ((0 == pixdata[full_pos]) && (min == edge_data[full_pos]))
					{
						pixdata[full_pos] = seed_identifier;
						need_flood_fill = true;
					}
					else {
						need_flood_fill = false;
					}
				}
				__syncthreads();

				// All threads check to see whether flood fill is needed.
				if (need_flood_fill)
				{
					// In the non-CUDA code, the growth of each contiguous region is not limited by the buffer zones.  Use flood-specific
					// variables here.
					int flood_region_width = xdiv;
					int flood_region_height = ydiv;
					int flood_x_offset = cell_i * xdiv;
					int flood_y_offset = cell_j * ydiv;
					if (flood_x_offset + flood_region_width > width)
					{
						flood_region_width = width - flood_x_offset;
					}
					if (flood_y_offset + flood_region_height > height)
					{
						flood_region_height = height - flood_y_offset;
					}
					region_flood_fill(width, height, flood_x_offset, flood_y_offset, flood_region_width, flood_region_height, min, seed_identifier, pixdata, edge_data, NULL, 0);
				}

				if (0 == idx)
				{
					if (need_flood_fill)
					{
						seed_identifier = seed_identifier + 1;
					}
					cell_pos = cell_pos + 1;
				}
				__syncthreads();
			}
		}
	}
}

__global__ void find_seeds_cell_grow(int width, int height, int xdiv, int ydiv, unsigned char* edge_data, int* temp_spdata, int* delta_spdata, bool* collision_data, int* stage_sizes)
{
	// CUDA kernel for growing regions in a cell to find the one that overtakes any others.
	// This version uses __device__ functions.
	// 
	// width and height are dimensions for the overall image, and map to the size of temp_spdata and edge_data.
	// xdiv and ydiv are the size of idividual cells.
	// edge_data is the gradient data for the image, where edges get higher numbers.
	// temp_spdata is the copy of spdata that is used for growing the superpixels.
	// delta_spdata is used for holding a temporary value before it is written to temp_spdata.
	// collision_data is used for tracking superpixel intersections when they happen.
	// stage_sizes holds the size of each superpixel.

	// This kernel will be called with each block corresponding to one cell.
	int n_cell = blockIdx.x;

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
	int cell_width = xdiv;
	int cell_height = ydiv;

	if ((x_offset < width) && (y_offset < height) && (x_offset >= 0) && (y_offset >= 0)) // Do basic dimension testing.  Applies for block, so __syncthread()-safe.
	{
		// Adjust portion dimensions if necessary.
		if (x_offset + cell_width > width)
		{
			cell_width = width - x_offset;
		}
		if (y_offset + cell_height > height)
		{
			cell_height = height - y_offset;
		}
		if ((cell_width > 0) && (cell_height > 0))
		{
			find_seeds_cell_grow_main_old(width, height, x_offset, y_offset, cell_width, cell_height, edge_data, temp_spdata, delta_spdata, collision_data, stage_sizes, NULL, 0);
		}
	}
}

__device__ void find_seeds_cell_grow_main_old(int width, int height, int x_offset, int y_offset, int cell_width, int cell_height, unsigned char* edge_data, int* temp_spdata, int* delta_spdata, bool* collision_data, int* stage_sizes, int* mask, int mask_identifier)
{
	// Main function for growing in cell.
	// Inherits same arguments as find_seeds_cell_grow, plus:
	// cell_width and cell_height are the dimensions of this cell (which may be different from xdiv and ydiv).
	// If present, mask and mask_identifier limit where growth may occur.

	// Shared data across the threads of this cell.
	__shared__ bool continue_step;
	__shared__ int numbers_seen[threads_per_block];
	__shared__ unsigned char level; // Current level
	__shared__ bool all_done; // Indicates that there are no longer multiple superpixels in this cell.

	int idx = threadIdx.x;
	int N_cell_elements = cell_width * cell_height; // Number of elements in the cell.

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
			// First, we calculate stage_sizes.
			calc_stage_sizes_old(width, N_cell_elements, cell_width, x_offset, y_offset, temp_spdata, stage_sizes);
			__syncthreads();

			// We need to look for cases where an open pixel would, on the next step, lead to a collision between
			// two or more regions.  These are set to true in collision_data.
			detect_collision_old(width, cell_width, cell_height, x_offset, y_offset, level, edge_data, temp_spdata, collision_data);
			__syncthreads();

			// We also need to check for cases where two superpixels have already grown into each other.
			detect_touching_regions(width, cell_width, cell_height, x_offset, y_offset, temp_spdata, collision_data, NULL, 0);

			// Now, the collision_data matrix defines all upcoming intersections between superpixels
			// that need to be addressed.

			// We need to go through and address each upcoming intersection by switching the values of
			// each superpixel to that of the larger of the ones that are part of the intersection (with
			// a tie going to the one with a smaller identifier).

			handle_collisions_old(width, cell_width, cell_height, x_offset, y_offset, temp_spdata, collision_data, stage_sizes);

			// All collision points have been addressed.
			// Now, grow each superpixel by one pixel.
			if (0 == idx)
			{
				continue_step = false; // Only sets to true if there is growth.
			}
			__syncthreads();

			// First, zero out the delta_spdata matrix.
			zero_out_cell_old(width, cell_width, cell_height, x_offset, y_offset, delta_spdata);
			__syncthreads();

			// Perform the growth step.
			growth_step_1(width, cell_width, cell_height, x_offset, y_offset, level, edge_data, temp_spdata, delta_spdata, mask, mask_identifier);
			__syncthreads();

			// Now that the growth step is done, update with the delta matrix.
			growth_step_2(width, cell_width, cell_height, x_offset, y_offset, temp_spdata, delta_spdata, &continue_step, NULL, 0);
			__syncthreads();

			// Last part of main loop:
			// Count unique region numbers.  If more than one, need to continue.
			if (0 == idx)
			{
				all_done = true;
			}
			__syncthreads();

			detect_multiple_regions(width, cell_width, cell_height, x_offset, y_offset, numbers_seen, temp_spdata, &all_done, NULL, 0);
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

__device__ void calc_stage_sizes_old(int width, int N_cell_elements, int cell_width, int x_offset, int y_offset, int* temp_spdata, int* stage_sizes)
{
	// Function to clear and then fill out the size of each numbered region in temp_dpdata.
	// width is the width of the overall image.
	// n_cell_elements is the number of elements in this cell.
	// cell_width is the width of this cell.
	// x_offset and y_offset are the x and y locations of the beginning of this cell within the image.
	// temp_spdata is the matrix holding superpixel identifier for each image pixel.
	// stage_sizes is a matrix holding the size of each superpixel, indexed by identifier.

	int idx = threadIdx.x;

	// Clear out stage_sizes
	for (int i = idx; i < N_cell_elements; i += threads_per_block) // Striding through array.
	{
		int local_x = i % cell_width; // x value within cell.
		int local_y = i / cell_width; // y value within cell.
		int x = x_offset + local_x; // Overall x value.
		int y = y_offset + local_y; // Overall y value.
		int pos = x + y * width; // Position of the element within the full array.
		stage_sizes[pos] = 0;
	}
	__syncthreads();

	// Calculate current sizes.
	for (int i = idx; i < N_cell_elements; i += threads_per_block) // Striding through portion of the array.
	{
		int local_x = i % cell_width; // x value within cell.
		int local_y = i / cell_width; // y value within cell.
		int x = x_offset + local_x; // Overall x value.
		int y = y_offset + local_y; // Overall y value.
		int pos = x + y * width; // Position of the element within the full array.
		int value = temp_spdata[pos];
		if (value > 0)
		{
			// value is now the cell-specific location to update.  This works because the number of superpixels in a cell
			// must be equal or smaller than the number of pixels in the cell.  In practice, this results in much more
			// memory being allocated in stage_sizes than is used.
			// By offsetting the position based on the cell offset, we ensure that we do not
			// write over stage size data being used by other blocks (cells).
			local_x = value % cell_width; // x value within cell.
			local_y = value / cell_width; // y value within cell.
			x = x_offset + local_x; // Overall x value.
			y = y_offset + local_y; // Overall y value.
			pos = x + y * width; // Position of the element within the full array.
			atomicAdd(&stage_sizes[pos], 1);
		}
	}
}

__device__ void detect_collision_old(int width, int cell_width, int cell_height, int x_offset, int y_offset, unsigned char level, unsigned char* edge_data, int* temp_spdata, bool* collision_data)
{
	// Function to detect that growth in the next iteration will cause two regions to touch.
	// Arguments have their values from find_seeds_cell_grow_main function.

	int idx = threadIdx.x;
	int N_cell_elements = cell_width * cell_height; // Number of elements in the cell.

	for (int i = idx; i < N_cell_elements; i += threads_per_block) // Striding through portion of the array.
	{
		int local_x = i % cell_width; // x value within cell.
		int local_y = i / cell_width; // y value within cell.
		int x = x_offset + local_x; // Overall x value.
		int y = y_offset + local_y; // Overall y value.
		int pos = x + y * width; // Position of the element within the full array.
		collision_data[pos] = false;
		if ((0 == temp_spdata[pos]) && (level >= edge_data[pos])) // Open position at the right level.
		{
			int value = 0; // Used to store value of contiguous superpixels detected.
			// Check adjacent pixels.

			// Top
			if (local_y > 0) // Don't look outside of cell.
			{
				value = temp_spdata[pos - width]; // No need to test, since this is first.
			}

			// Bottom
			if (local_y < cell_height - 1) // Don't look outside of cell.
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
			if (local_x > 0) // Don't look outside of cell.
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
			if (local_x < cell_width - 1) // Don't look outside of cell.
			{
				if (value > 0) // At least one has already been seen.
				{
					if ((temp_spdata[pos + 1] > 0) && (temp_spdata[pos + 1] != value))
					{
						collision_data[pos] = true;
					}
				} // No need for the else here, since there is at most one adjacent superpixel value, so no collision.
			}
		}
	}
}

__device__ void detect_touching_regions(int width, int cell_width, int cell_height, int x_offset, int y_offset, int* temp_spdata, bool* collision_data, int* mask, int mask_identifier)
{
	int idx = threadIdx.x;
	int N_cell_elements = cell_width * cell_height; // Number of elements in the cell.

	for (int i = idx; i < N_cell_elements; i += threads_per_block) // Striding through cell.
	{
		int portion_x = i % cell_width; // x value within cell.
		int portion_y = i / cell_width; // y value within cell.
		int x = x_offset + portion_x; // Overall x value.
		int y = y_offset + portion_y; // Overall y value.
		int pos = x + y * width; // Position of the element within the full array.
		if ((NULL == mask) || (mask_identifier == mask[pos]))
		{
			if ((false == collision_data[pos]) && (temp_spdata[pos] > 0)) // Not a current collision location, and is part of a superpixel region.
			{
				int value = 0;
				// Check adjacent pixels.

				// Top
				if (portion_y > 0) // Don't look outside of defined region.
				{
					if ((NULL == mask) || (mask_identifier == mask[pos - width]))
					{
						value = temp_spdata[pos - width];
						if ((value > 0) && (value != temp_spdata[pos]))
						{
							collision_data[pos] = true;
						}
					}
				}

				// Bottom
				if (portion_y < cell_height - 1) // Don't look outside of defined region.
				{
					if ((NULL == mask) || (mask_identifier == mask[pos + width]))
					{
						value = temp_spdata[pos + width];
						if ((value > 0) && (value != temp_spdata[pos]))
						{
							collision_data[pos] = true;
						}
					}
				}

				// Left
				if (portion_x > 0) // Don't look outside of defined region.
				{
					if ((NULL == mask) || (mask_identifier == mask[pos - 1]))
					{
						value = temp_spdata[pos - 1];
						if ((value > 0) && (value != temp_spdata[pos]))
						{
							collision_data[pos] = true;
						}
					}
				}

				// Right
				if (portion_x < cell_width - 1) // Don't look outside of defined region.
				{
					if ((NULL == mask) || (mask_identifier == mask[pos + 1]))
					{
						value = temp_spdata[pos + 1];
						if ((value > 0) && (value != temp_spdata[pos]))
						{
							collision_data[pos] = true;
						}
					}
				}
			}
		}
	}
}

__device__ void handle_collisions_old(int width, int cell_width, int cell_height, int x_offset, int y_offset, int* temp_spdata, bool* collision_data, int* stage_sizes)
{
	// Resolve each instance of a collision between superpixels, by having one superpixel absorbed
	// by the other, which involves replacing the identifiers of one with the other.  The larger
	// superpixel prevails, and in the case of a tie, the one with a lower identifier value wins.

	__shared__ int target_identifier; // Identifier for absorption step.
	__shared__ int contig_identifier[4]; // The four potential identifiers that are contiguous to a location.
	__shared__ int contig_size[4]; // Holds the size of each superpixel included in contig_identifiers.
	__shared__ int iter_pos; // Position value while iterating through the region.

	int idx = threadIdx.x;
	int N_cell_elements = cell_width * cell_height; // Number of elements in the cell.

	// We'll use a single thread to track the progress.
	if (0 == idx)
	{
		iter_pos = 0;
	}
	__syncthreads();

	while (iter_pos < N_cell_elements) // Step through all elements in the region.
	{

		// Find the identifiers that are contiguous to the current location defined by iter_pos.
		if (0 == idx)
		{
			target_identifier = find_target_identifier_old(width, cell_width, cell_height, x_offset, y_offset, iter_pos, contig_identifier, contig_size, temp_spdata, collision_data, stage_sizes);
		}
		__syncthreads();

		// Absorb the contiguous superpixels into the target_identifier superpixel.
		if (target_identifier > 0) // True if there is a target identifier that all contiguous superpixels should be changed to.
		{
			// The value of target_identifier is the value that all identifiers that show up in contig_identifier should be changed to.
			for (int i = idx; i < N_cell_elements; i += threads_per_block) // Striding through cell.
			{
				// Find the position for this element.
				int local_portion_x = i % cell_width; // x value within cell.
				int local_portion_y = i / cell_width; // y value within cell.
				int local_x = x_offset + local_portion_x; // Overall x value.
				int local_y = y_offset + local_portion_y; // Overall y value.
				int local_pos = local_x + local_y * width; // Position of the element within the full array.
				for (int index = 0; index < 4; ++index)
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
				int index = 0;
				while (index < 4)
				{
					if ((contig_identifier[index] > 0) && (contig_identifier[index] != target_identifier))
					{
						number_shifted += contig_size[index];
						int size_pos = x_offset + contig_identifier[index] % cell_width + width * (y_offset + contig_identifier[index] / cell_width);
						stage_sizes[size_pos] = 0;
					}
					index++;
				}
				int size_pos = x_offset + target_identifier % cell_width + width * (y_offset + target_identifier / cell_width);
				stage_sizes[size_pos] += number_shifted;
			}
		}
		__syncthreads();

		if (0 == idx)
		{
			iter_pos++;
		}
		__syncthreads();
	}
}

__device__ int find_target_identifier_old(int width, int cell_width, int cell_height, int x_offset, int y_offset, int iter_pos, int* contig_identifier, int* contig_size, int* temp_spdata, bool* collision_data, int* stage_sizes)
{
	// Find the target_identifier which corresponds to the identifier that all contiguous superpixels should be changed to.
	// int iter_pos is the position within the cell being processed.
	// This function is single-threaded and should be called only by the first thread in the block.

	int local_x = iter_pos % cell_width; // x value within cell.
	int local_y = iter_pos / cell_width; // y value within cell.
	int x = x_offset + local_x; // Overall x value.
	int y = y_offset + local_y; // Overall y value.
	int pos = x + y * width; // Position of the element within the full array.
	int value;
	int target_identifier = -1; // Return value.  Negative means targets does not exist.

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
		if (local_y > 0) // Don't look outside of cell.
		{
			value = temp_spdata[pos - width];
			if (value > 0)
			{
				contig_identifier[0] = value; // No need to test, since this is first.
				int size_pos = x_offset + value % cell_width + width * (y_offset + value / cell_width);
				contig_size[0] = stage_sizes[size_pos]; // Get the superpixel size from stage_sizes.
			}
		}

		// Bottom
		if (local_y < (cell_height - 1)) // Look only inside of cell.
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
					int size_pos = x_offset + value % cell_width + width * (y_offset + value / cell_width);
					contig_size[index] = stage_sizes[size_pos]; // Get the superpixel size from stage_sizes.
				}
			}
		}

		// Left
		if (local_x > 0) // Don't look outside of cell.
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
					int size_pos = x_offset + value % cell_width + width * (y_offset + value / cell_width);
					contig_size[index] = stage_sizes[size_pos]; // Get the superpixel size from stage_sizes.
				}
			}
		}

		// Right
		if (local_x < cell_width - 1) // Don't look outside of cell.
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
					int size_pos = x_offset + value % cell_width + width * (y_offset + value / cell_width);
					contig_size[index] = stage_sizes[size_pos]; // Get the superpixel size from stage_sizes.
				}
			}
		}

		// Now, there should be two or more identifiers in contig_identifier if collision_data for this location is true.
		// It may not be, though, because of previous absorb steps.

		if (contig_identifier[1] > 0)
		{
			int index;
			// We know that contig_identifier contains multiple identifiers that need to be combined.
			// contig_size contains the size of each superpixel corresponding to those identifiers.

			// Find the identifier of the largest superpixel.
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
	}
	return target_identifier;
}

__device__ void zero_out_cell_old(int width, int cell_width, int cell_height, int x_offset, int y_offset, int* matrix)
{
	// Zero out values in a defined cell portion of a superpixel materix.
	int idx = threadIdx.x;
	int N_cell_elements = cell_width * cell_height; // Number of elements in the cell.

	for (int i = idx; i < N_cell_elements; i += threads_per_block) // Striding through portion of the array.
	{
		int portion_x = i % cell_width; // x value within window of the target region.
		int portion_y = i / cell_width; // y value within window of the target region.
		int x = x_offset + portion_x; // Overall x value.
		int y = y_offset + portion_y; // Overall y value.
		int pos = x + y * width; // Position of the element within the full array.
		matrix[pos] = 0;
	}
}

__device__ void growth_step_1(int width, int cell_width, int cell_height, int x_offset, int y_offset, unsigned char level, unsigned char* edge_data, int* temp_spdata, int* delta_spdata, int* mask, int mask_identifier)
{
	// Peform single growth step in a cell, with level variable defining the value of comparison in edge_data.

	int idx = threadIdx.x;
	bool ignore_mask = (NULL == mask) || (mask_identifier < 1);
	int N_cell_elements = cell_width * cell_height; // Number of elements in the cell.
	for (int i = idx; i < N_cell_elements; i += threads_per_block) // Striding through portion of the array.
	{
		int portion_x = i % cell_width; // x value within window of the target region.
		int portion_y = i / cell_width; // y value within window of the target region.
		int x = x_offset + portion_x; // Overall x value.
		int y = y_offset + portion_y; // Overall y value.
		int pos = x + y * width; // Position of the element within the full array.
		int value = 0;

		if ((0 == temp_spdata[pos]) && (edge_data[pos] <= level) && (ignore_mask || (mask_identifier == mask[pos])))
		{
			// Rely on the fact that there should not be growth to a given location from different superpixels, due to collision steps.

			// Top
			if (portion_y > 0) // Don't look outside of cell.
			{
				if ((NULL == mask) || (mask_identifier == mask[pos - width]))
				{
					value = temp_spdata[pos - width];
					if (value > 0)
					{
						delta_spdata[pos] = value; // Write to delta matrix, for later consolidation with temp_spdata.
					}
				}
			}

			// Bottom
			if ((value == 0) && (portion_y < cell_height - 1)) // Don't look outside of cell.
			{
				if ((NULL == mask) || (mask_identifier == mask[pos + width]))
				{
					value = temp_spdata[pos + width];
					if (value > 0)
					{
						delta_spdata[pos] = value;
					}
				}
			}

			// Left
			if ((value == 0) && (portion_x > 0)) // Don't look outside of cell.
			{
				if ((NULL == mask) || (mask_identifier == mask[pos - 1]))
				{
					value = temp_spdata[pos - 1];
					if (value > 0)
					{
						delta_spdata[pos] = value;
					}
				}
			}

			// Right
			if ((value == 0) && (portion_x < cell_width - 1)) // Don't look outside of cell.
			{
				if ((NULL == mask) || (mask_identifier == mask[pos + 1]))
				{
					value = temp_spdata[pos + 1];
					if (value > 0)
					{
						delta_spdata[pos] = value;
					}
				}
			}
		}
	}
}

__device__ void growth_step_2(int width, int cell_width, int cell_height, int x_offset, int y_offset, int* temp_spdata, int* delta_spdata, bool* continue_step_ptr, int* mask, int mask_identifier)
{
	int idx = threadIdx.x;
	int N_cell_elements = cell_width * cell_height; // Number of elements in the cell.

	for (int i = idx; i < N_cell_elements; i += threads_per_block) // Striding through portion of the array.
	{
		int portion_x = i % cell_width; // x value within window of the target region.
		int portion_y = i / cell_width; // y value within window of the target region.
		int x = x_offset + portion_x; // Overall x value.
		int y = y_offset + portion_y; // Overall y value.
		int pos = x + y * width; // Position of the element within the full array.
		if ((NULL == mask) || (mask_identifier == mask[pos]))
		{
			int value = delta_spdata[pos];
			if (value > 0) // This value was updated in the first growth step.
			{
				temp_spdata[pos] = value;
				if (false == *continue_step_ptr)
				{
					*continue_step_ptr = true;
				}
			}
		}
	}
}

__device__ void detect_multiple_regions(int width, int cell_width, int cell_height, int x_offset, int y_offset, int* numbers_seen, int* temp_spdata, bool* all_done_ptr, int* mask, int mask_identifier)
{
	// Set all_done (through all_done_ptr) to false if there is more than one region identifier in the cell.
	// all_done (shared variable) should be set to true before calling this function.

	int idx = threadIdx.x;
	int N_cell_elements = cell_width * cell_height; // Number of elements in the cell.

	numbers_seen[idx] = 0; // Initialize value.
	for (int i = idx; *all_done_ptr && (i < N_cell_elements); i += threads_per_block) // Striding through portion of the array.
	{
		int portion_x = i % cell_width; // x value within window of the target region.
		int portion_y = i / cell_width; // y value within window of the target region.
		int x = x_offset + portion_x; // Overall x value.
		int y = y_offset + portion_y; // Overall y value.
		int pos = x + y * width; // Position of the element within the full array.
		if ((NULL == mask) || (mask_identifier == mask[pos]))
		{
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
						*all_done_ptr = false; // There is more than one numbered region.  This thread happened to hit two.
					}
				}
			}
		}
	}
	__syncthreads(); // Wait for all threads to finish.

	// Use single thread to check on numbers seen.
	if (*all_done_ptr && (0 == idx)) // all_done being true here may not be final, but if it is false, then no need for this loop.
	{
		int value = 0;
		for (int i = 0; *all_done_ptr && (i < threads_per_block); ++i)
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
						*all_done_ptr = false;
					}
				}
			}
		}
	}
}

__global__ void resolve_seeds_initial(int width, int height, int xdiv, int ydiv, int buffer, int* region_spdata, int* temp_spdata, int* d_seeds)
{
	// CUDA kernel for determining the initial seed for each cell.  Assumes that there is only one superpixel region in the cell.
	// width and height are dimensions for the overall image, and map to the size of spdata and temp_spdata.
	// xdiv and ydiv are the size of idividual cells.
	// buffer is the space on the outside of the cell where a seed cannot be located.
	// region_spdata is the original numbered regions at min level.
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

	// Call main function with zero offset, because this covers the full image.
	resolve_seeds_main_old(width, height, 0, 0, xdiv, ydiv, n_x_cells, buffer, region_spdata, temp_spdata, d_seeds);
}

__device__ void resolve_seeds_main_old(int width, int height, int x_region_offset, int y_region_offset, int dx, int dy, int n_x_cells, int buffer, int* region_spdata, int* temp_spdata, int* d_seeds)
{
	// CUDA kernel for determining the seed for each cell.  Assumes that there is only one superpixel region in the cell.
	// width and height are dimensions for the overall image, and map to the size of local_pixdata nd edge_data.
	// x_region_offset and y_region_offset are the beginning of the region being split (offsets within the overall image matrix).
	// dx and dy are the size of idividual cells within the region to be split.
	// n_x_cells is the number of columns in the cell.  Number of rows is not a parameter, because the number of blocks implies it.
	// buffer is the number of pixels on each edge that are not eligible to be seeds.
	// region_spdata is the original numbered regions at min level.
	// temp_spdata is the copy of region_spdata that was used for growing the superpixels, so it should have only one numbered region.
	// d_seeds is the array into which seed results are placed.  They are ordered by cell number, with x followed by y.

	// This kernel will be called with each block corresponding to one cell.
	int n_cell = blockIdx.x;

	__shared__ int pos_array[threads_per_block]; // One position calculated for each thread.
	__shared__ int identifier;

	//	Region calculations (no buffer is used here):
	int cell_x = n_cell % n_x_cells;
	int cell_y = n_cell / n_x_cells;
	int x_cell_offset = cell_x * dx + x_region_offset;
	int y_cell_offset = cell_y * dy + y_region_offset;

	int idx = threadIdx.x;
	if (0 == idx)
	{
		d_seeds[2 * n_cell] = -1; // Default values in case this kernel is called with bad parameters.
		d_seeds[2 * n_cell + 1] = -1;
	}
	__syncthreads();

	if ((x_cell_offset < width) && (y_cell_offset < height) && (x_cell_offset >= 0) && (y_cell_offset >= 0)) // Do basic dimension testing (same for block, so not __syncthreads problem).
	{
		// Adjust portion dimensions if necessary.
		if (x_cell_offset + dx > width)
		{
			dx = width - x_cell_offset;
		}
		if (y_cell_offset + dy > height)
		{
			dy = height - y_cell_offset;
		}
		if ((dx > 0) && (dy > 0))
		{
			int N_cell_elements = dx * dy; // Number of elements in the cell.

			if (0 == idx)
			{
				// Find the identifier for the remaining superpixel in this cell.
				int value = 0;
				for (int i = 0; (i < N_cell_elements) && (value == 0); ++i) // Going through cell elements, stop when non-zero found.
				{
					int local_x = i % dx; // x value within cell.
					int local_y = i / dx; // y value within cell.
					int x = x_cell_offset + local_x; // Overall x value.
					int y = y_cell_offset + local_y; // Overall y value.
					int pos = x + y * width; // Position of this element within the full array.
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
				// Find the position with lowest number that maps to an identifier in spdata.
				int first_pos_found = -1; // Variable is local to this thread.  Initialize to invalid value.
				for (int local_pos = idx; (local_pos < N_cell_elements) && (first_pos_found < 0); local_pos += threads_per_block) // Striding through cell.
				{
					int local_x = local_pos % dx; // x value within the cell.
					int local_y = local_pos / dx; // y value within the cell.
					int x = x_cell_offset + local_x; // Overall x value.
					int y = y_cell_offset + local_y; // Overall y value.
					int pos = x + y * width; // Position of the element within the full array.
					if ((local_x >= buffer) && (local_y >= buffer) && (local_x < dx - buffer) && (local_y < dy - buffer))
					{
						int value = region_spdata[pos];
						if (identifier == value)
						{
							first_pos_found = pos; // Assign value of position in overall image.
						}
					}
				}
				pos_array[idx] = first_pos_found; // Write to the location that is specific to this thread.
				__syncthreads(); // Wait for all threads to finish.			

				// Now, reduce the answer array to one element.  Don't assume that threads_per_block is a power of two (although it almost always is).
				unsigned int working_set = threads_per_block;
				while (working_set > 1)
				{
					int half_ceiling = (working_set + 1) / 2; // Handle case where working_set is odd.
					if ((idx < half_ceiling) && (idx + half_ceiling < working_set))
					{
						if (pos_array[idx + half_ceiling] >= 0) // Screen out negative values for threads that failed to find a value.
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

bool c_watershed_grow(int width, int height, SuperPixel* list_head)
{
	// CUDA version of watershed growth algorithm.
	// width and height define the size of the image.
	// list_head is a pointer to the first SuperPixel.

	bool ret = true;
	bool* c_result = BoolArray(1, 1, false);
	cudaError_t cudaStatus;
	int blocksPerGrid = (width * height + threads_per_block - 1) / threads_per_block;

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
		ret = CopyFromHost(host_data->GetData(), width * height, host_data->GetDeviceData());
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
			while (continue_with_level)
			{
				ResetBoolArray(ResultsArray, blocksPerGrid, 1, false);
				watershed_step_1 << < blocksPerGrid, threads_per_block >> > (level, width, height, host_data->GetDeviceData(), list_head->GetGradient()->GetCData(), ResultsArray);
				cudaDeviceSynchronize();
				cudaStatus = cudaGetLastError();
				if (cudaStatus != cudaSuccess)
				{
					throw std::runtime_error("Failed in watershed_step.\n");
					return false;
				}
				// Now, change the new identifier values from negative to positive.  They are initially stored as negatives, to distinguish them from the pre-state.
				watershed_step_2 << < blocksPerGrid, threads_per_block >> > (width, height, host_data->GetDeviceData(), ResultsArray);
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
		ret = CopyToHost(host_data->GetDeviceData(), width * height, host_data->GetData());
		if (false == ret)
		{
			throw std::runtime_error("Failed to copy superpixel data to host.\n");
			return ret;
		}

		// Need to recompute the bounding box, size, and last_level for each SuperPixel.		
		ret = c_find_bbox(width, height, list_head);
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

std::vector <PointPair> c_find_seeds(int width, int height, GradData* edge, SPixelData* pixeldata, int xdiv, int ydiv, int buffer)
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

	find_seeds_cell_grow << <num_cells, threads_per_block >> > (width, height, xdiv, ydiv, edge_data, temp_spdata, delta_spdata, collision_data, stage_sizes);
	cudaDeviceSynchronize();
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed in find_seeds_cell_grow.\n");
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

	resolve_seeds_initial << <num_cells, threads_per_block >> > (width, height, xdiv, ydiv, buffer, spdata, temp_spdata, d_seeds);
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

// Below is the new code for calculating seeds for superpixels that need to be split.

__global__ void resolve_seeds_split(int width, int height, int* d_geometry, int buffer, int* region_spdata, int* temp_spdata, int* d_seeds, int max_nx, int max_ny, int* mask, int* d_identifiers)
{
	// CUDA kernel for determining the split candidate seed for each cell.  Assumes that there is only one superpixel region in the cell.
	// width and height are dimensions for the overall image, and map to the size of spdata and temp_spdata.
	// x_region_offset and y_region_offset are the beginning of the region being split (offsets within the overall image matrix). 
	// xdiv and ydiv are the size of idividual cells.
	// n_x_cells is the number of columns of cells in the region to be evaluated.  Number of rows is not a parameter, because the number of blocks implies it.
	// buffer is the space on the outside of the cell where a seed cannot be located.
	// region_spdata is the original numbered regions at min level.
	// temp_spdata is the copy of spdata that was used for growing the superpixels.
	// d_seeds is the array into which seed results are placed.  They are ordered by cell number, with x followed by y.

	// This kernel will be called with each block corresponding to one cell.
	int n_cell = blockIdx.x;
	int n_superpixel = blockIdx.y;
	int idx = threadIdx.x;

	// Get mask_identifier.
	int mask_identifier = d_identifiers[n_superpixel];

	// Make local_stage_sizes a shared memory (shared within the block).
	extern __shared__ int stage_sizes[];  // Array sized for the largest number of regions encountered in any cell.

	// Get geometry from d_geometry.
	int offset_x = d_geometry[10 * n_superpixel];
	int offset_y = d_geometry[10 * n_superpixel + 1];
	int region_width = d_geometry[10 * n_superpixel + 2];
	int region_height = d_geometry[10 * n_superpixel + 3];
	int nx = d_geometry[10 * n_superpixel + 4];
	int ny = d_geometry[10 * n_superpixel + 5];
	int dx = d_geometry[10 * n_superpixel + 6];
	int dy = d_geometry[10 * n_superpixel + 7];

	// Find row for d_seeds:
	int* local_seeds = &d_seeds[2 * n_superpixel * max_nx * max_ny];

	//	Region calculations:
	int cell_i = n_cell % nx;
	int cell_j = n_cell / nx;
	int cell_x_offset = buffer + cell_i * dx + offset_x; // The offset within the overall image to the left side of this cell (within the border).
	int cell_y_offset = buffer + cell_j * dy + offset_y; // The offset within the overall image to the top of this cell (within the border).
	int cell_width = dx - 2 * buffer; // The width of the cell, not counting the border.
	int cell_height = dy - 2 * buffer; // The height of the cell, not counting the border.

	resolve_seeds_main(width, height, offset_x, offset_y, dx, dy, nx, ny, buffer, region_spdata, temp_spdata, local_seeds, mask, mask_identifier);
}

__device__ void resolve_seeds_main(int width, int height, int x_region_offset, int y_region_offset, int dx, int dy, int nx, int ny, int buffer, int* region_spdata, int* temp_spdata, int* d_seeds, int* mask, int mask_identifier)
{
	// CUDA kernel for determining the seed for each cell.  Assumes that there is only one superpixel region in the cell.
	// width and height are dimensions for the overall image, and map to the size of local_pixdata nd edge_data.
	// x_region_offset and y_region_offset are the beginning of the region being split (offsets within the overall image matrix).
	// dx and dy are the size of individual cells within the region to be split.
	// nx is the number of cell columns in the region.  ny is the number of rows of cells in the region.
	// buffer is the number of pixels on each edge that are not eligible to be seeds.
	// region_spdata is the original numbered regions at min level.
	// temp_spdata is the copy of region_spdata that was used for growing the superpixels, so it should have only one numbered region.
	// d_seeds is the array into which seed results are placed.  They are ordered by cell number, with x followed by y.

	// This kernel will be called with each block corresponding to one cell.
	int n_cell = blockIdx.x;
	int idx = threadIdx.x;

	__shared__ int pos_array[threads_per_block]; // One position calculated for each thread.
	__shared__ int identifier;

	//	Region calculations (no buffer is used here):
	int cell_x = n_cell % nx;
	int cell_y = n_cell / nx;
	int x_cell_offset = cell_x * dx + x_region_offset;
	int y_cell_offset = cell_y * dy + y_region_offset;

	if (0 == idx)
	{
		d_seeds[2 * n_cell] = -1; // Default values.
		d_seeds[2 * n_cell + 1] = -1;
	}
	__syncthreads();

	if (n_cell < nx * ny) // Confirm this cell is within the set defined by nx and ny.
	{
		if ((x_cell_offset < width) && (y_cell_offset < height) && (x_cell_offset >= 0) && (y_cell_offset >= 0)) // Do basic dimension testing (same for block, so not __syncthreads problem).
		{
			// Adjust portion dimensions if necessary.
			if (x_cell_offset + dx > width)
			{
				dx = width - x_cell_offset;
			}
			if (y_cell_offset + dy > height)
			{
				dy = height - y_cell_offset;
			}
			if ((dx > 0) && (dy > 0))
			{
				int N_cell_elements = dx * dy; // Number of elements in the cell.

				if (0 == idx)
				{
					// Find the identifier for the remaining superpixel in this cell.
					int value = 0;
					for (int i = 0; (i < N_cell_elements) && (value == 0); ++i) // Going through cell elements, stop when non-zero found.
					{
						int local_x = i % dx; // x value within cell.
						int local_y = i / dx; // y value within cell.
						int x = x_cell_offset + local_x; // Overall x value.
						int y = y_cell_offset + local_y; // Overall y value.
						int pos = x + y * width; // Position of this element within the full array.
						if ((NULL == mask) || (mask_identifier == mask[pos]))
						{
							value = temp_spdata[pos];
						}
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
					// Find the position with lowest number that maps to an identifier in region_spdata.
					int first_pos_found = -1; // Variable is local to this thread.  Initialize to invalid value.
					for (int local_pos = idx; (local_pos < N_cell_elements) && (first_pos_found < 0); local_pos += threads_per_block) // Striding through cell.
					{
						int local_x = local_pos % dx; // x value within the cell.
						int local_y = local_pos / dx; // y value within the cell.
						int x = x_cell_offset + local_x; // Overall x value.
						int y = y_cell_offset + local_y; // Overall y value.
						int pos = x + y * width; // Position of the element within the full array.
						if ((NULL == mask) || (mask_identifier == mask[pos]))
						{
							if ((local_x >= buffer) && (local_y >= buffer) && (local_x < dx - buffer) && (local_y < dy - buffer))
							{
								int value = region_spdata[pos];
								if (identifier == value)
								{
									first_pos_found = pos; // Assign value of position in overall image.
								}
							}
						}
					}
					pos_array[idx] = first_pos_found; // Write to the location that is specific to this thread.
					__syncthreads(); // Wait for all threads to finish.			

					// Now, reduce the answer array to one element.  Don't assume that threads_per_block is a power of two (although it almost always is).
					unsigned int working_set = threads_per_block;
					while (working_set > 1)
					{
						int half_ceiling = (working_set + 1) / 2; // Handle case where working_set is odd.
						if ((idx < half_ceiling) && (idx + half_ceiling < working_set))
						{
							if (pos_array[idx + half_ceiling] >= 0) // Screen out negative values for threads that failed to find a value.
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
}

__device__ void zero_out_cell(int width, int cell_width, int cell_height, int x_offset, int y_offset, int* matrix, int* mask, int mask_identifier)
{
	// Zero out values in a defined cell portion of a superpixel materix.
	int idx = threadIdx.x;
	int N_cell_elements = cell_width * cell_height; // Number of elements in the cell.

	for (int i = idx; i < N_cell_elements; i += threads_per_block) // Striding through portion of the array.
	{
		int portion_x = i % cell_width; // x value within window of the target region.
		int portion_y = i / cell_width; // y value within window of the target region.
		int x = x_offset + portion_x; // Overall x value.
		int y = y_offset + portion_y; // Overall y value.
		int pos = x + y * width; // Position of the element within the full array.
		if ((NULL == mask) || (mask_identifier == mask[pos]))
		{
			matrix[pos] = 0;
		}
	}
}

__device__ void detect_collision(int width, int cell_width, int cell_height, int x_offset, int y_offset, unsigned char level, unsigned char* edge_data, int* temp_spdata, bool* collision_data, int* mask, int mask_identifier)
{
	// Function to detect that growth in the next iteration will cause two regions to touch.
	// Arguments have their values from find_seeds_cell_grow_main function.

	int idx = threadIdx.x;
	int N_cell_elements = cell_width * cell_height; // Number of elements in the cell.

	for (int i = idx; i < N_cell_elements; i += threads_per_block) // Striding through portion of the array.
	{
		int local_x = i % cell_width; // x value within cell.
		int local_y = i / cell_width; // y value within cell.
		int x = x_offset + local_x; // Overall x value.
		int y = y_offset + local_y; // Overall y value.
		int pos = x + y * width; // Position of the element within the full array.
		if ((NULL == mask) || (mask_identifier == mask[pos]))
		{
			collision_data[pos] = false;
			if ((0 == temp_spdata[pos]) && (level >= edge_data[pos])) // Open position at the right level.
			{
				int value = 0; // Used to store value of contiguous superpixels detected.
				// Check adjacent pixels.

				// Top
				if (local_y > 0) // Don't look outside of cell.
				{
					if ((NULL == mask) || (mask_identifier == mask[pos - width]))
					{
						value = temp_spdata[pos - width]; // No need to test, since this is first.
					}
				}

				// Bottom
				if (local_y < cell_height - 1) // Don't look outside of cell.
				{
					if ((NULL == mask) || (mask_identifier == mask[pos + width]))
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
				}

				// Left
				if (local_x > 0) // Don't look outside of cell.
				{
					if ((NULL == mask) || (mask_identifier == mask[pos - 1]))
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
				}

				// Right
				if (local_x < cell_width - 1) // Don't look outside of cell.
				{
					if ((NULL == mask) || (mask_identifier == mask[pos + 1]))
					{
						if (value > 0) // At least one has already been seen.
						{
							if ((temp_spdata[pos + 1] > 0) && (temp_spdata[pos + 1] != value))
							{
								collision_data[pos] = true;
							}
						} // No need for the else here, since there is at most one adjacent superpixel value, so no collision.
					}
				}
			}
		}
	}
}

__device__ int find_target_identifier(int width, int cell_width, int cell_height, int x_offset, int y_offset, int iter_pos, int* contig_identifier, int* contig_size, int* temp_spdata, bool* collision_data, int* stage_sizes, int* mask, int mask_identifier)
{
	// Find the target_identifier which corresponds to the identifier that all contiguous superpixels should be changed to.
	// int iter_pos is the position within the cell being processed.
	// This function is single-threaded and should be called only by the first thread in the block.

	int local_x = iter_pos % cell_width; // x value within cell.
	int local_y = iter_pos / cell_width; // y value within cell.
	int x = x_offset + local_x; // Overall x value.
	int y = y_offset + local_y; // Overall y value.
	int pos = x + y * width; // Position of the element within the full array.
	int value;
	int target_identifier = -1; // Return value.  Negative means targets does not exist.

	contig_identifier[0] = 0;
	contig_identifier[1] = 0;
	contig_identifier[2] = 0;
	contig_identifier[3] = 0;
	contig_size[0] = 0;
	contig_size[1] = 0;
	contig_size[2] = 0;
	contig_size[3] = 0;

	if ((NULL == mask) || (mask_identifier == mask[pos]))
	{
		if (collision_data[pos]) // Find identifiers for contiguous pixels.
		{
			// Top
			if (local_y > 0) // Don't look outside of cell.
			{
				if ((NULL == mask) || (mask_identifier == mask[pos - width]))
				{
					value = temp_spdata[pos - width];
					if (value > 0)
					{
						contig_identifier[0] = value; // No need to test, since this is first.
						contig_size[0] = stage_sizes[value - 1]; // Get the superpixel size from stage_sizes.
					}
				}
			}

			// Bottom
			if (local_y < (cell_height - 1)) // Look only inside of cell.
			{
				if ((NULL == mask) || (mask_identifier == mask[pos + width]))
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
							contig_size[index] = stage_sizes[value - 1]; // Get the superpixel size from stage_sizes.
						}
					}
				}
			}

			// Left
			if (local_x > 0) // Don't look outside of cell.
			{
				if ((NULL == mask) || (mask_identifier == mask[pos - 1]))
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
							contig_size[index] = stage_sizes[value - 1]; // Get the superpixel size from stage_sizes.
						}
					}
				}
			}

			// Right
			if (local_x < cell_width - 1) // Don't look outside of cell.
			{
				if ((NULL == mask) || (mask_identifier == mask[pos + 1]))
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
							contig_size[index] = stage_sizes[value - 1]; // Get the superpixel size from stage_sizes.
						}
					}
				}
			}

			// Now, there should be two or more identifiers in contig_identifier if collision_data for this location is true.
			// It may not be, though, because of previous absorb steps.

			if (contig_identifier[1] > 0)
			{
				int index;
				// We know that contig_identifier contains multiple identifiers that need to be combined.
				// contig_size contains the size of each superpixel corresponding to those identifiers.

				// Find the identifier of the largest superpixel.
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
		}
	}
	return target_identifier;
}

__device__ void handle_collisions(int width, int cell_width, int cell_height, int x_offset, int y_offset, int* temp_spdata, bool* collision_data, int* stage_sizes, int* mask, int mask_identifier)
{
	// Resolve each instance of a collision between superpixels, by having one superpixel absorbed
	// by the other, which involves replacing the identifiers of one with the other.  The larger
	// superpixel prevails, and in the case of a tie, the one with a lower identifier value wins.

	__shared__ int target_identifier; // Identifier for absorption step.
	__shared__ int contig_identifier[4]; // The four potential identifiers that are contiguous to a location.
	__shared__ int contig_size[4]; // Holds the size of each superpixel included in contig_identifiers.
	__shared__ int iter_pos; // Position value while iterating through the region.

	int idx = threadIdx.x;
	int N_cell_elements = cell_width * cell_height; // Number of elements in the cell.

	// We'll use a single thread to track the progress.
	if (0 == idx)
	{
		iter_pos = 0;
	}
	__syncthreads();

	while (iter_pos < N_cell_elements) // Step through all elements in the region.
	{
		// Find the identifiers that are contiguous to the current location defined by iter_pos.
		if (0 == idx)
		{
			target_identifier = find_target_identifier(width, cell_width, cell_height, x_offset, y_offset, iter_pos, contig_identifier, contig_size, temp_spdata, collision_data, stage_sizes, mask, mask_identifier);
		}
		__syncthreads();

		// Absorb the contiguous superpixels into the target_identifier superpixel.
		if (target_identifier > 0) // True if there is a target identifier that all contiguous superpixels should be changed to.
		{
			// The value of target_identifier is the value that all identifiers that show up in contig_identifier should be changed to.
			for (int i = idx; i < N_cell_elements; i += threads_per_block) // Striding through cell.
			{
				// Find the position for this element.
				int local_portion_x = i % cell_width; // x value within cell.
				int local_portion_y = i / cell_width; // y value within cell.
				int local_x = x_offset + local_portion_x; // Overall x value.
				int local_y = y_offset + local_portion_y; // Overall y value.
				int local_pos = local_x + local_y * width; // Position of the element within the full array.
				if ((NULL == mask) || (mask_identifier == mask[local_pos]))
				{
					for (int index = 0; index < 4; ++index)
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
			}
			__syncthreads();

			// Update the sizes
			if (0 == idx)
			{
				int number_shifted = 0;
				int index = 0;
				while (index < 4)
				{
					if ((contig_identifier[index] > 0) && (contig_identifier[index] != target_identifier))
					{
						stage_sizes[contig_identifier[index] - 1] = 0;
					}
					index++;
				}
				stage_sizes[target_identifier - 1] += number_shifted;
			}
		}
		__syncthreads();

		if (0 == idx)
		{
			iter_pos++;
		}
		__syncthreads();
	}
}

__device__ void calc_stage_sizes(int width, int N_cell_elements, int cell_width, int x_offset, int y_offset, int* temp_spdata, int* stage_sizes, int max_identifiers, int* mask, int mask_identifier)
{
	// Function to clear and then fill out the size of each numbered region in temp_dpdata.
	// width is the width of the overall image.
	// n_cell_elements is the number of elements in this cell.
	// cell_width is the width of this cell.
	// x_offset and y_offset are the x and y locations of the beginning of this cell within the image.
	// temp_spdata is the matrix holding superpixel identifier for each image pixel.
	// stage_sizes is a matrix holding the size of each superpixel, indexed by identifier.
	// max_identifiers is the number of integers allocated for stage_sizes.

	int idx = threadIdx.x;

	// Clear out stage_sizes
	for (int i = idx; i < max_identifiers; i += threads_per_block) // Striding through array.
	{
		stage_sizes[i] = 0;
	}
	__syncthreads();

	// Calculate current sizes.
	for (int i = idx; i < N_cell_elements; i += threads_per_block) // Striding through portion of the array.
	{
		int local_x = i % cell_width; // x value within cell.
		int local_y = i / cell_width; // y value within cell.
		int x = x_offset + local_x; // Overall x value.
		int y = y_offset + local_y; // Overall y value.
		int pos = x + y * width; // Position of the element within the full array.
		if ((NULL == mask) || (mask_identifier == mask[pos]))
		{
			int value = temp_spdata[pos];
			if (value > 0)
			{
				atomicAdd(&stage_sizes[value - 1], 1);
			}
		}
	}
}

__device__ void find_seeds_cell_grow_main(int width, int height, int x_offset, int y_offset, int cell_width, int cell_height, unsigned char* edge_data, int* temp_spdata, int* delta_spdata, bool* collision_data, int* stage_sizes, int max_identifiers, int* mask, int mask_identifier)
{
	// Main function for growing in cell.
	// Inherits same arguments as find_seeds_cell_grow, plus:
	// cell_width and cell_height are the dimensions of this cell (which may be different from xdiv and ydiv).
	// If present, mask and mask_identifier limit where growth may occur.

	// Shared data across the threads of this cell.
	__shared__ bool continue_step;
	__shared__ int numbers_seen[threads_per_block];
	__shared__ unsigned char level; // Current level
	__shared__ bool all_done; // Indicates that there are no longer multiple superpixels in this cell.

	int idx = threadIdx.x;
	int N_cell_elements = cell_width * cell_height; // Number of elements in the cell.

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
			// First, we calculate stage_sizes.
			calc_stage_sizes(width, N_cell_elements, cell_width, x_offset, y_offset, temp_spdata, stage_sizes, max_identifiers, mask, mask_identifier);
			__syncthreads();

			// We need to look for cases where an open pixel would, on the next step, lead to a collision between
			// two or more regions.  These are set to true in collision_data.
			detect_collision(width, cell_width, cell_height, x_offset, y_offset, level, edge_data, temp_spdata, collision_data, mask, mask_identifier);
			__syncthreads();

			// We also need to check for cases where two superpixels have already grown into each other.
			detect_touching_regions(width, cell_width, cell_height, x_offset, y_offset, temp_spdata, collision_data, mask, mask_identifier);

			// Now, the collision_data matrix defines all upcoming intersections between superpixels
			// that need to be addressed.

			// We need to go through and address each upcoming intersection by switching the values of
			// each superpixel to that of the larger of the ones that are part of the intersection (with
			// a tie going to the one with a smaller identifier).

			handle_collisions(width, cell_width, cell_height, x_offset, y_offset, temp_spdata, collision_data, stage_sizes, mask, mask_identifier);

			// All collision points have been addressed.
			// Now, grow each superpixel by one pixel.
			if (0 == idx)
			{
				continue_step = false; // Only sets to true if there is growth.
			}
			__syncthreads();

			// First, zero out the delta_spdata matrix.
			zero_out_cell(width, cell_width, cell_height, x_offset, y_offset, delta_spdata, mask, mask_identifier);
			__syncthreads();

			// Perform the growth step.
			growth_step_1(width, cell_width, cell_height, x_offset, y_offset, level, edge_data, temp_spdata, delta_spdata, mask, mask_identifier);
			__syncthreads();

			// Now that the growth step is done, update with the delta matrix.
			growth_step_2(width, cell_width, cell_height, x_offset, y_offset, temp_spdata, delta_spdata, &continue_step, mask, mask_identifier);
			__syncthreads();

			// Last part of main loop:
			// Count unique region numbers.  If more than one, need to continue.
			if (0 == idx)
			{
				all_done = true;
			}
			__syncthreads();

			detect_multiple_regions(width, cell_width, cell_height, x_offset, y_offset, numbers_seen, temp_spdata, &all_done, mask, mask_identifier);
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

__global__ void split_cell_grow(int width, int height, int* d_geometry, int buffer, unsigned char* edge_data, int* temp_spdata, int* delta_spdata, bool* collision_data, int max_identifiers, int* mask, int* d_identifiers)
{

	// CUDA kernel for growing regions in a cell, as part of the splitting function, to find the one that overtakes any others.
	// width and height are dimensions for the overall image, and map to the size of temp_spdata and edge_data.
	// offset_x and offset_y are the distance into image where the cells for interrogation are located.
	// nx and ny are the number of cells in the x and y directions to be evaluated.
	// dx and dy are the dimensions of each cell to be evaluated.
	// edge_data is the gradient data for the image, where edges get higher numbers.
	// temp_spdata is the copy of spdata that is used for growing the superpixels.
	// delta_spdata is used for holding a temporary value before it is written to temp_spdata.
	// collision_data is used for tracking superpixel intersections when they happen.
	// max_identifiers is the most identifiers in any cell.  This number of integers is the size of the dynamically-allocated stage_sizes.
	// mask is a matrix with the current superpixel data.
	// mask_identifier is the identifier of the current superpixel.

	// This kernel will be called with each block corresponding to one cell.
	int n_cell = blockIdx.x;
	int n_superpixel = blockIdx.y;
	int idx = threadIdx.x;

	// Get mask_identifier.
	int mask_identifier = d_identifiers[n_superpixel];

	// Make local_stage_sizes a shared memory (shared within the block).
	extern __shared__ int stage_sizes[];  // Array sized for the largest number of regions encountered in any cell.

	// Get geometry from d_geometry.
	int offset_x = d_geometry[10 * n_superpixel];
	int offset_y = d_geometry[10 * n_superpixel + 1];
	int region_width = d_geometry[10 * n_superpixel + 2];
	int region_height = d_geometry[10 * n_superpixel + 3];
	int nx = d_geometry[10 * n_superpixel + 4];
	int ny = d_geometry[10 * n_superpixel + 5];
	int dx = d_geometry[10 * n_superpixel + 6];
	int dy = d_geometry[10 * n_superpixel + 7];

	//	Region calculations:
	int cell_i = n_cell % nx;
	int cell_j = n_cell / nx;
	int cell_x_offset = buffer + cell_i * dx + offset_x; // The offset within the overall image to the left side of this cell (within the border).
	int cell_y_offset = buffer + cell_j * dy + offset_y; // The offset within the overall image to the top of this cell (within the border).
	int cell_width = dx - 2 * buffer; // The width of the cell, not counting the border.
	int cell_height = dy - 2 * buffer; // The height of the cell, not counting the border.

	if ((offset_x < width) && (offset_y < height) && (offset_x >= 0) && (offset_y >= 0)) // Do basic dimension testing.  Applies for block, so __syncthread()-safe.
	{
		// Adjust portion dimensions if necessary.
		if (cell_x_offset + cell_width > width)
		{
			cell_width = width - cell_x_offset;
		}
		if (cell_y_offset + cell_height > height)
		{
			cell_height = height - cell_y_offset;
		}
		if ((cell_width > 0) && (cell_height > 0))
		{
			find_seeds_cell_grow_main(width, height, cell_x_offset, cell_y_offset, cell_width, cell_height, edge_data, temp_spdata, delta_spdata, collision_data, stage_sizes, max_identifiers, mask, mask_identifier);
		}
	}
}

__global__ void find_split_cell_regions(int width, int height, int* d_geometry, int buffer, int* mask, int* d_identifiers, int* local_pixdata, unsigned char* edge_data)
{
	// CUDA kernel for finding the potential split cell regions within regions of the image corresponding to the min value.
	// d_geometry defines offset_x, offset_y, region_width, region_height, nx, ny, dx, and dy.
	// Each region is a subset of the overall image, offset by x_region_offset and y_region_offset.  Within the region are cells
	// that are of dx by dy dimensions, but with borders of size buffer on each side.  There are nx cells in the x direction
	// in the region, and ny cells in the y direction in the region.
	// 
	// width and height are dimensions for the overall image, and map to the size of local_pixdata nd edge_data.
	// x_region_offset and y_region_offset are the beginning of the region being split (offsets within the overall image matrix).
	// dx and dy are the size of individual cells within the region to be split.
	// nx and ny are the number of cells in each direction that define the cells where we look for new seeds.
	// buffer is the number of pixels on each edge that are not eligible to be seeds.
	// mask is the superpixel identifier data for the image, where an integer identifier is assigned to each pixel.
	// d_identifiers contains the values in mask that act as a mask for the discovery of seed locations.
	// local_pixdata is a blank copy of mask, used locally for examining growth of superpixels at different seed locations.
	// edge_data is the gradient data for the image, where edges get higher numbers.

	// This kernel will be called with each block corresponding to one cell.
	// The x direction for blocks defines the cell number, within the nx by ny framework.
	// The y direction for blocks defines the index number for the superpixel being calculated.

	int n_cell = blockIdx.x;
	int n_superpixel = blockIdx.y;
	int idx = threadIdx.x;

	// Shared data across the threads of this cell.
	__shared__ unsigned char min; // Minimum value of edge_data in this region.
	__shared__ unsigned char min_array[threads_per_block]; // One minimum calculated for each thread.
	__shared__ int cell_pos; // Position within the cell.
	__shared__ int seed_identifier; // Identifier to record in each contiguous area of min.
	__shared__ bool need_flood_fill; // Indicates to threads whether a flood-fill is needed.

	// Get geometry from d_geometry.
	int offset_x = d_geometry[10 * n_superpixel];
	int offset_y = d_geometry[10 * n_superpixel + 1];
	int region_width = d_geometry[10 * n_superpixel + 2];
	int region_height = d_geometry[10 * n_superpixel + 3];
	int nx = d_geometry[10 * n_superpixel + 4];
	int ny = d_geometry[10 * n_superpixel + 5];
	int dx = d_geometry[10 * n_superpixel + 6];
	int dy = d_geometry[10 * n_superpixel + 7];

	// Get mask_identifier.
	int mask_identifier = d_identifiers[n_superpixel];

	if ((nx > 0) && (ny > 0) && (n_cell < nx * ny)) // Test for empty entry, or cell outside of range.
	{
		//	Region calculations:
		int cell_i = n_cell % nx;
		int cell_j = n_cell / nx;
		int cell_x_offset = buffer + cell_i * dx + offset_x; // The offset within the overall image to the left side of this cell (within the border).
		int cell_y_offset = buffer + cell_j * dy + offset_y; // The offset within the overall image to the top of this cell (within the border).
		int cell_width = dx - 2 * buffer; // The width of the cell, within the border.
		int cell_height = dy - 2 * buffer; // The height of the cell, within the border.

		// First task is to calculate the minimum value of edge_data in the cell, using the mask.
		if ((offset_x < width) && (offset_y < height) && (offset_x >= 0) && (offset_y >= 0)) // Do basic dimension testing.
		{
			// Adjust cell dimensions if necessary.
			if (cell_x_offset + cell_width > width)
			{
				cell_width = width - cell_x_offset;
			}
			if (cell_y_offset + cell_height > height)
			{
				cell_height = height - cell_y_offset;
			}
			if ((cell_width > 0) && (cell_height > 0))
			{
				int N_cell_elements = cell_width * cell_height; // Number of elements in the cell.

				// First task is to calculate the minimum value of edge_data in the cell.
				min_uChar_Array_portion_main(edge_data, width, cell_x_offset, cell_y_offset, cell_width, cell_height, min_array, mask, mask_identifier);

				// Put final answer in min.
				if (0 == idx)
				{
					min = min_array[0];
				}
				__syncthreads(); // Now, we have the minimum value of this cell in edge_data in min.

				// Now, we need to go through the pixels of the region sequentially.
				if (0 == idx)
				{
					cell_pos = 0;
					seed_identifier = 1;
				}
				__syncthreads();
				while (cell_pos < N_cell_elements)
				{
					int cell_x = cell_pos % cell_width;
					int cell_y = cell_pos / cell_width;
					int full_pos = cell_x + cell_x_offset + (cell_y + cell_y_offset) * width;
					if (0 == idx)
					{
						if ((mask_identifier == mask[full_pos]) && (0 == local_pixdata[full_pos]) && (min == edge_data[full_pos]))
						{
							local_pixdata[full_pos] = seed_identifier;
							need_flood_fill = true;
						}
						else {
							need_flood_fill = false;
						}
					}
					__syncthreads();

					// All threads check to see whether flood fill is needed.
					if (need_flood_fill)
					{
						// In the non-CUDA code, the growth of each contiguous region is not limited by the buffer zones.  Use flood-specific
						// variables here.
						int flood_region_width = dx;
						int flood_region_height = dy;
						int flood_x_offset = cell_i * dx + offset_x;
						int flood_y_offset = cell_j * dy + offset_y;
						if (flood_x_offset + flood_region_width > width)
						{
							flood_region_width = width - flood_x_offset;
						}
						if (flood_y_offset + flood_region_height > height)
						{
							flood_region_height = height - flood_y_offset;
						}
						region_flood_fill(width, height, flood_x_offset, flood_y_offset, flood_region_width, flood_region_height, min, seed_identifier, local_pixdata, edge_data, mask, mask_identifier);
					}

					if (0 == idx)
					{
						if (need_flood_fill)
						{
							seed_identifier = seed_identifier + 1;
						}
						cell_pos += 1;
					}
					__syncthreads();
				}
			}
		}
	}
}


__global__ void split_grow(int width, int height, int* d_geometry, unsigned char* image_data, int color_channels, int* mask, unsigned char* edge, int* temp_pixdata, int* d_seeds, int max_nx, int max_ny, int* d_identifiers)
{
	// This kernel goes through each cell of a superpixel (one superpixel per block) and tests each candidate seed
	//   against the original seed for the superpixel.  The best seed, if any, will be written out to the first cell
	//   location of d_seeds (overwriting the input seed, which is no longer needed).
	// 
	// width and height are dimensions for the overall image and the edge, mask and temp_pixdata matrices.
	// d_geometry is an array with all geometric information for the cells, ordered by superpixel.
	// color_channels is the number of color channels in the image.  This is generally 3, but in the future it could
	//   be more (but only the first three are used).
	// mask is a matrix of superpixel identifiers.
	// d_seeds is an array of candidate seed values for each cell of each superpixel.
	// max_nx and max_nx are the maximum number of cells in each direction that may be investigated for each superpixel.
	// d_identifiers is an array of the identifiers for each superpixel, which are used along with the mask matrix to
	//   determine which locations belong to each superpixel.
	// 
	// This kernel will be called with each block corresponding to one superpixel.
	// The x direction for blocks defines the cell number, within the nx by ny framework.
	// The y direction for blocks defines the index number for the superpixel being calculated.

	int n_superpixel = blockIdx.x;
	int idx = threadIdx.x;

	// Memory shared across the block (superpixel).
	__shared__ int level; // Current level
	__shared__ bool continue_at_level;
	__shared__ double Color_R_1[threads_per_block]; // Variables with _1 correspond to the initial seed.
	__shared__ double Color_G_1[threads_per_block];
	__shared__ double Color_B_1[threads_per_block];
	__shared__ double Color_R_2[threads_per_block]; // Variables with _2 correspond to the candidate seed corresponding to a cell.
	__shared__ double Color_G_2[threads_per_block];
	__shared__ double Color_B_2[threads_per_block];
	__shared__ double pixel_count[threads_per_block];
	__shared__ double pixel_count_1[threads_per_block];
	__shared__ double pixel_count_2[threads_per_block];
	__shared__ float Color_Diff;
	__shared__ int out_seed_x;
	__shared__ int out_seed_y;

	if (0 == idx)
	{
		Color_Diff = 0;
		out_seed_x = -1;
		out_seed_y = -1;
	}

	// Get geometry from d_geometry.
	int offset_x = d_geometry[10 * n_superpixel];
	int offset_y = d_geometry[10 * n_superpixel + 1];
	int region_width = d_geometry[10 * n_superpixel + 2];
	int region_height = d_geometry[10 * n_superpixel + 3];
	int nx = d_geometry[10 * n_superpixel + 4];
	int ny = d_geometry[10 * n_superpixel + 5];
	int dx = d_geometry[10 * n_superpixel + 6];
	int dy = d_geometry[10 * n_superpixel + 7];
	int initial_seed_x = d_geometry[10 * n_superpixel + 8]; // This is the initial seed for input, but we write the output seed here as well.
	int initial_seed_y = d_geometry[10 * n_superpixel + 9];

	PointPair resolved_seed;
	resolved_seed.x = -1;
	resolved_seed.y = -1;

	// Get mask_identifier.
	int mask_identifier = d_identifiers[n_superpixel];

	// Find row for d_seeds:
	int* local_seeds = &d_seeds[2 * n_superpixel * max_nx * max_ny];

	int n_region_elements = region_width * region_height;
	int initial_seed_pos = initial_seed_x + initial_seed_y * width; // Location within temp_pixdata of initial seed.

	if ((initial_seed_x < offset_x) || (initial_seed_x >= offset_x + region_width) || (initial_seed_y < offset_y) || (initial_seed_y >= offset_y + region_height) || (mask_identifier != mask[initial_seed_pos]) || (color_channels < 3)) // This should not happen.
	{
		if (0 == idx)
		{
			local_seeds[0] = -1;
			local_seeds[1] = -1;
			d_geometry[10 * n_superpixel + 8] = -1;
			d_geometry[10 * n_superpixel + 9] = -1;
		}
		return;
	}

	for (int n_cell = 0; n_cell < nx * ny; ++n_cell) // Loop through each cell to compare the corresponding seed to the initial seed.
	{
		int candidate_seed_x = local_seeds[n_cell * 2];
		int candidate_seed_y = local_seeds[n_cell * 2 + 1];
		int candidate_seed_pos = candidate_seed_x + candidate_seed_y * width;// Location within temp_pixdata of candidate seed.

		if ((candidate_seed_x >= offset_x) && (candidate_seed_y >= offset_y) &&
			(candidate_seed_x < offset_x + region_width) && (candidate_seed_y < offset_y + region_height) &&
			(mask_identifier == mask[candidate_seed_pos])) // See if candidate seed is on mask.  Not all cells are used, those that aren't will have -1, -1 seed location.
		{
			// First step is to reset the portion of temp_pixdata corresponding to the mask.
			for (int region_pos = idx; region_pos < n_region_elements; region_pos += threads_per_block) // Stride through region.
			{
				int region_x = region_pos % region_width;
				int region_y = region_pos / region_width;
				int x = offset_x + region_x;
				int y = offset_y + region_y;
				int pos = x + y * width;
				if ((NULL == mask) || (mask_identifier == mask[pos]))
				{
					if (initial_seed_pos == pos)
					{
						temp_pixdata[pos] = 1; // Used for initial seed.
					}
					else if (candidate_seed_pos == pos)
					{
						temp_pixdata[pos] = 2; // Used for candidate seed.
					}
					else
					{
						temp_pixdata[pos] = 0;
					}
				}
			}
			__syncthreads();

			// Now, grow the superpixels on the mask.
			if (0 == idx)
			{
				level = 0;
				continue_at_level = true;
			}
			__syncthreads();

			while (level <= 255)
			{
				while (continue_at_level)
				{
					// Set continue_at_level to false.
					if (0 == idx)
					{
						continue_at_level = false;
					}
					__syncthreads();

					// Process step 1.
					for (int region_pos = idx; region_pos < n_region_elements; region_pos += threads_per_block) // Stride through region.
					{
						int region_x = region_pos % region_width;
						int region_y = region_pos / region_width;
						int x = offset_x + region_x;
						int y = offset_y + region_y;
						int pos = x + y * width;
						int value = 0;
						if ((NULL == mask) || (mask_identifier == mask[pos]))
						{
							// Is this an open position on the pixel map?
							if (0 == temp_pixdata[pos])
							{
								// Is the gradient level at or below the target level?
								if (level >= edge[pos])
								{
									int smallest_adjacent = -1; // This will be the smallest positive adjacent superpixel identifier.						
									if ((y > 0) && ((NULL == mask) || (mask_identifier == mask[pos - width]))) // Can't look above when at the top of the image
									{
										value = temp_pixdata[pos - width];
										if ((value > 0) && ((smallest_adjacent < 0) || (value < smallest_adjacent)))
										{
											smallest_adjacent = value;
										}
									}
									if ((y < height - 1) && ((NULL == mask) || (mask_identifier == mask[pos + width]))) // Make sure not at bottom of image.
									{
										value = temp_pixdata[pos + width];
										if ((value > 0) && ((smallest_adjacent < 0) || (value < smallest_adjacent)))
										{
											smallest_adjacent = value;
										}
									}
									if ((x > 0) && ((NULL == mask) || (mask_identifier == mask[pos - 1]))) // Make sure not at left side.
									{
										value = temp_pixdata[pos - 1];
										if ((value > 0) && ((smallest_adjacent < 0) || (value < smallest_adjacent)))
										{
											smallest_adjacent = value;
										}
									}
									if ((x < width - 1) && ((NULL == mask) || (mask_identifier == mask[pos + 1]))) // Make sure not on the right side.
									{
										value = temp_pixdata[pos + 1];
										if ((value > 0) && ((smallest_adjacent < 0) || (value < smallest_adjacent)))
										{
											smallest_adjacent = value;
										}
									}
									if (smallest_adjacent > 0) // This will be true if there is an adjacent value larger than zero, and it will always defer to the smaller identifier.
									{
										temp_pixdata[pos] = -smallest_adjacent; // Set to negative, so other threads don't read this as a superpixel until after this iteration is complete.
										if (!continue_at_level)
										{
											continue_at_level = true; // Having expanded a superpixel, another pass is needed at this level.
										}
									}
								}
							}
						}
					}
					__syncthreads();

					// Process step 2.
					if (continue_at_level) // No need for step 2 if no elements were changed in step 1.
					{
						for (int region_pos = idx; region_pos < n_region_elements; region_pos += threads_per_block) // Stride through region.
						{
							int region_x = region_pos % region_width;
							int region_y = region_pos / region_width;
							int x = offset_x + region_x;
							int y = offset_y + region_y;
							int pos = x + y * width;
							if ((NULL == mask) || (mask_identifier == mask[pos]))
							{
								int value = temp_pixdata[pos];
								if (value < 0)
								{
									temp_pixdata[pos] = -value;
								}
							}
						}
					}
					__syncthreads();
					// If any elements were changed, continue_step will be true, and this level will be repeated.
				}
				__syncthreads();
				if (0 == idx)
				{
					level += 1;
					continue_at_level = true;
				}
				__syncthreads();
			}

			// Growth is done for this candidate seed.  Evaluate sizes and color variation.
			Color_R_1[idx] = 0;
			Color_G_1[idx] = 0;
			Color_B_1[idx] = 0;
			pixel_count_1[idx] = 0;
			Color_R_2[idx] = 0;
			Color_G_2[idx] = 0;
			Color_B_2[idx] = 0;
			pixel_count_2[idx] = 0;
			pixel_count[idx] = 0;

			for (int region_pos = idx; region_pos < n_region_elements; region_pos += threads_per_block) // Stride through region.
			{
				int region_x = region_pos % region_width;
				int region_y = region_pos / region_width;
				int x = offset_x + region_x;
				int y = offset_y + region_y;
				int pos = x + y * width;
				int image_pos = color_channels * pos;
				if ((NULL == mask) || (mask_identifier == mask[pos]))
				{
					pixel_count[idx] += 1;
					if (1 == temp_pixdata[pos]) // Initial seed.
					{
						// Record color components squared.
						Color_R_1[idx] += image_data[image_pos] * image_data[image_pos];
						Color_G_1[idx] += image_data[image_pos + 1] * image_data[image_pos + 1];
						Color_B_1[idx] += image_data[image_pos + 2] * image_data[image_pos + 2];
						pixel_count_1[idx] += 1;
					}
					else if (2 == temp_pixdata[pos])
					{
						Color_R_2[idx] += image_data[image_pos] * image_data[image_pos];
						Color_G_2[idx] += image_data[image_pos + 1] * image_data[image_pos + 1];
						Color_B_2[idx] += image_data[image_pos + 2] * image_data[image_pos + 2];
						pixel_count_2[idx] += 1;
					}

				}
			}
			__syncthreads();
			// Now, reduce vectors to a single value.
			unsigned int working_set = threads_per_block;
			while (working_set > 1)
			{
				int half_ceiling = (working_set + 1) / 2; // Handle case where working_set is odd.
				if ((idx < half_ceiling) && (idx + half_ceiling < working_set))
				{
					pixel_count[idx] += pixel_count[idx + half_ceiling]; // Total pixels in mask
					Color_R_1[idx] += Color_R_1[idx + half_ceiling]; // Red value for seed 1 superpixel
					Color_G_1[idx] += Color_G_1[idx + half_ceiling];
					Color_B_1[idx] += Color_B_1[idx + half_ceiling];
					pixel_count_1[idx] += pixel_count_1[idx + half_ceiling]; // Pixels in seed 1 superpixel
					Color_R_2[idx] += Color_R_2[idx + half_ceiling]; // Red value for seed 2 superpixel
					Color_G_2[idx] += Color_G_2[idx + half_ceiling];
					Color_B_2[idx] += Color_B_2[idx + half_ceiling];
					pixel_count_2[idx] += pixel_count_2[idx + half_ceiling]; // Pixels in seed 2 superpixel
				}
				working_set = half_ceiling;
				__syncthreads(); // Each pass through the reduction process needs to wait for all threads to complete.
			}
			// Use sizes to determine whether this is a valid candidate for further evaluation.
			if ((0 == idx) && (pixel_count_1[0] > 0.1 * pixel_count[0]) && (pixel_count_2[0] > 0.1 * pixel_count[0]))
			{
				// Take square root of summed squares, put result in value for idx == 1.
				Color_R_1[0] = sqrt(Color_R_1[0] / pixel_count_1[0]);
				Color_G_1[0] = sqrt(Color_G_1[0] / pixel_count_1[0]);
				Color_B_1[0] = sqrt(Color_B_1[0] / pixel_count_1[0]);
				Color_R_2[0] = sqrt(Color_R_2[0] / pixel_count_2[0]);
				Color_G_2[0] = sqrt(Color_G_2[0] / pixel_count_2[0]);
				Color_B_2[0] = sqrt(Color_B_2[0] / pixel_count_2[0]);

				// Calculate color difference.

				double R_diff = Color_R_1[0] - Color_R_2[0];
				double G_diff = Color_G_1[0] - Color_G_2[0];
				double B_diff = Color_B_1[0] - Color_B_2[0];
				float current_color_diff = sqrt(R_diff * R_diff + G_diff * G_diff + B_diff * B_diff);
				if (current_color_diff > Color_Diff)
				{
					out_seed_x = candidate_seed_x;
					out_seed_y = candidate_seed_y;
					Color_Diff = current_color_diff;
				}
			}
			__syncthreads();
		}
	}
	// Done with analysis of all cells.  Send back output seed (in place of where the initial seed was stored in d_geometry).
	// If we did not find a suitable replacement, then the output seed will be -1, -1.  The caller can use this to determine
	// whether a suitable new seed was found, so it won't need to compare to the input seed to see if it changed.
	if (0 == idx)
	{
		local_seeds[0] = out_seed_x;
		local_seeds[1] = out_seed_y;
		local_seeds[2] = 50;
		local_seeds[3] = 50;
		d_geometry[10 * n_superpixel + 8] = out_seed_x;
		d_geometry[10 * n_superpixel + 9] = out_seed_y;
	}
	__syncthreads();
}

std::vector<PointPair> c_split_seeds(SuperPixel* head, std::vector<int>split_identifiers, ImageData* image)
{
	// Function to find new seeds for superpixels that need to be split.
	// head is the first superpixel in the current set.
	// split_identifiers is the set of identifiers for the superpixels that need to be split.
	// image is the ImageData structure for the underlying image, which includes color data.

	cudaError_t cudaStatus;

	int width = head->GetPixelData()->GetWidth();
	int height = head->GetPixelData()->GetHeight();
	head->GetPixelData()->SyncToDevice(); // Ensure that pixel data is on the device.
	int* mask = head->GetPixelData()->GetDeviceData();
	unsigned char* edge_data = head->GetGradient()->GetCData();

	int* regions_pixdata = IntArray(width, height, true, 0); // Working space to be modified in calculating seeds.

	int max_nx = 7;
	int max_ny = 7;
	int buffer = 2;

	int color_channels = image->GetColorChannels();
	unsigned char* image_data = image->GetCData();

	std::vector<PointPair> add_seeds;
	add_seeds.clear();
	int num_identifiers = split_identifiers.size();
	if (num_identifiers > 0)
	{
		// h_geometry holds these int values for each superpixel:
		// offset_x, offset_y, region_width, region_height, nx, ny, dx, dy, seed_x, seed_y
		int* h_geometry = (int*)malloc(10 * sizeof(int) * num_identifiers);
		if (NULL == h_geometry)
		{
			throw std::runtime_error("Failed to allocate memory for h_geometry in c_split_seeds.\n");
			return add_seeds;
		}

		// h_identifiers holds the identifier values of each superpixel.
		int* h_identifiers = (int*)malloc(sizeof(int) * num_identifiers);
		if (NULL == h_identifiers)
		{
			throw std::runtime_error("Failed to allocate memory for h_identifiers in c_split_seeds.\n");
			return add_seeds;
		}

		std::vector<int>::iterator it;
		int superpixel_count = 0;
		for (it = split_identifiers.begin(); it != split_identifiers.end(); ++it)
		{
			int local_identifier = *it;
			SuperPixel* current = head->GetByIdentifier(local_identifier);
			h_identifiers[superpixel_count] = local_identifier;

			RectQuad local_window = current->GetWindow();
			int nx, ny, dx, dy;
			if ((local_window.x1 - local_window.x0) >= max_nx * 10)
			{
				nx = 7;
				dx = (local_window.x1 - local_window.x0) / max_nx;
			}
			else {
				nx = (local_window.x1 - local_window.x0) / 10;
				if (nx < 1)
				{
					nx = 1;
				}
				dx = (local_window.x1 - local_window.x0) / nx;
			}
			if ((local_window.y1 - local_window.y0) >= max_ny * 10)
			{
				ny = 7;
				dy = (local_window.y1 - local_window.y0) / max_ny;
			}
			else {
				ny = (local_window.y1 - local_window.y0) / 10;
				if (ny < 1)
				{
					ny = 1;
				}
				dy = (local_window.y1 - local_window.y0) / ny;
			}

			if ((nx < 2) && (ny < 2))
			{
				nx = 0;
				ny = 0;
				dx = 0;
				dy = 0;
			}

			h_geometry[10 * superpixel_count] = local_window.x0;
			h_geometry[10 * superpixel_count + 1] = local_window.y0;
			h_geometry[10 * superpixel_count + 2] = local_window.x1 - local_window.x0 + 1;
			h_geometry[10 * superpixel_count + 3] = local_window.y1 - local_window.y0 + 1;
			h_geometry[10 * superpixel_count + 4] = nx;
			h_geometry[10 * superpixel_count + 5] = ny;
			h_geometry[10 * superpixel_count + 6] = dx;
			h_geometry[10 * superpixel_count + 7] = dy;
			h_geometry[10 * superpixel_count + 8] = current->GetSeed().x;
			h_geometry[10 * superpixel_count + 9] = current->GetSeed().y;
			superpixel_count++;
		}

		// d_geometry and d_identifiers are on-device counterparts of h_geometry and h_identifiers.
		// Allocate memory on device, and then copy memory over.
		int* d_geometry = NULL;
		cudaStatus = cudaMalloc(&d_geometry, 10 * sizeof(int) * num_identifiers);
		if (cudaStatus != cudaSuccess)
		{
			throw std::runtime_error("Failed to allocate memory for d_geometry in c_split_seeds.\n");
			return add_seeds;
		}
		cudaStatus = cudaMemcpy(d_geometry, h_geometry, 10 * sizeof(int) * num_identifiers, cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess)
		{
			throw std::runtime_error("Failed to copy h_geometry to device in c_split_seeds.\n");
			return add_seeds;
		}

		int* d_identifiers = NULL;
		cudaStatus = cudaMalloc(&d_identifiers, sizeof(int) * num_identifiers);
		if (cudaStatus != cudaSuccess)
		{
			throw std::runtime_error("Failed to allocate memory for d_identifiers in c_split_seeds.\n");
			return add_seeds;
		}
		cudaStatus = cudaMemcpy(d_identifiers, h_identifiers, sizeof(int) * num_identifiers, cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess)
		{
			throw std::runtime_error("Failed to copy h_identifiers to device in c_split_seeds.\n");
			return add_seeds;
		}

		dim3 numBlocks(max_nx * max_ny, num_identifiers);
		find_split_cell_regions << <numBlocks, threads_per_block >> > (width, height, d_geometry, buffer, mask, d_identifiers, regions_pixdata, edge_data);
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess)
		{
			throw std::runtime_error("Failed in find_split_cell_regions.\n");
			return add_seeds;
		}

		//std::string filename = "D:\\temp\\regions.png";
		//WriteOutIntArray(regions_pixdata, width, height, filename, 0, 100);


		// Need to find the region with the largest identifier, to be able to allocate memory for stage_sizes.
		int* max_array = IntArray(threads_per_block, 1, true, 0);  // *** This should probably be shared memory within the block.
		int max_identifiers = 0;
		max_Int_Array << <1, threads_per_block >> > (regions_pixdata, width, height, max_array);
		cudaStatus = cudaMemcpy(&max_identifiers, max_array, sizeof(int), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess)
		{
			throw std::runtime_error("Failed to copy max value from device in c_split_seeds.\n");
			return add_seeds;
		}

		// Now, max_identifiers should hold the largest number of distinct regions there are in any cell in the image.

		// Set up buffer arrays for the process of growing the superpixels, the sizes of the superpixels, and to find the intersection locations.


		int* temp_pixdata = IntArray(width, height, false);
		if (NULL == temp_pixdata)
		{
			throw std::runtime_error("Unable to allocate temp_pixdata on device in c_find_split_seeds.\n");
			return add_seeds;
		}
		int* delta_spdata = IntArray(width, height, false);
		if (NULL == delta_spdata)
		{
			throw std::runtime_error("Unable to allocate delta_spdata on device in c_find_split_seeds.\n");
			return add_seeds;
		}
		bool* collision_data = BoolArray(width, height, false);
		if (NULL == collision_data)
		{
			throw std::runtime_error("Unable to allocate collision_data on device in c_find_split_seeds.\n");
			return add_seeds;
		}
		// Copy cell region data to temp_pixdata for resolving which region will determine the seed for the cell.
		cudaStatus = cudaMemcpy(temp_pixdata, regions_pixdata, width * height * sizeof(int), cudaMemcpyDeviceToDevice);
		if (cudaStatus != cudaSuccess)
		{
			std::cout << "Error copying memory in c_find_split_seeds.\n";
			return add_seeds;
		}

		// Protect against (highly unlikely) case of so many max_identifiers in a cell that too much shared memory would be allocated.
		if (max_identifiers > 10000)
		{
			throw std::runtime_error("Too many regions in a cell for splitting.\n");
			return add_seeds;
		}

		//std::string filename = "D:\\temp\\regions_1.png";
		//WriteOutIntArray(temp_pixdata, width, height, filename, 0, 100);

		split_cell_grow << <numBlocks, threads_per_block, max_identifiers * sizeof(int) >> > (width, height, d_geometry, 0, edge_data, temp_pixdata, delta_spdata, collision_data, max_identifiers, mask, d_identifiers);
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess)
		{
			throw std::runtime_error("Failed in split_cell_grow_old.\n");
			return add_seeds;
		}

		//filename = "D:\\temp\\regions_2.png";
		//WriteOutIntArray(temp_pixdata, width, height, filename, 0, 100);

		// Now we need to find the identifier from temp_spdata, and use it to get the earliest match of that identifier from
		// spdata, using a row-first process.  That is then the seed location for the cell.

		int* d_seeds = IntArray(2 * max_nx * max_ny * num_identifiers, 1, false); // Array to hold seed values on CUDA device.
		if (NULL == d_seeds)
		{
			throw std::runtime_error("Unable to allocate d_seeds on device in c_find_seeds.\n");
			return add_seeds;
		}
		int* h_seeds = (int*)malloc(2 * max_nx * max_ny * num_identifiers * sizeof(int));// Array to hold seed values on host.
		if (NULL == h_seeds)
		{
			throw std::runtime_error("Unable to allocate h_seeds on host in c_find_seeds.\n");
			return add_seeds;
		}

		// Determine the candidate seed for each cell of each superpixel.
		resolve_seeds_split << <numBlocks, threads_per_block >> > (width, height, d_geometry, buffer, regions_pixdata, temp_pixdata, d_seeds, max_nx, max_ny, mask, d_identifiers);
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess)
		{
			throw std::runtime_error("Failed in resolve_seeds_split.\n");
			return add_seeds;
		}

		//filename = "D:\\temp\\candidate_seeds.png";
		//WriteOutIntArray(d_seeds, 2 * max_nx * max_ny, num_identifiers, filename, 0, 255);

		// Grow each candidate seed, along with the initial seed, to determine which is the best candidate for each superpixel.
		split_grow << <num_identifiers, threads_per_block >> > (width, height, d_geometry, image_data, color_channels, mask, edge_data, temp_pixdata, d_seeds, max_nx, max_ny, d_identifiers);
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess)
		{
			throw std::runtime_error("Failed in split_grow.\n");
			return add_seeds;
		}

		//filename = "D:\\temp\\d_geometry.png";
		//WriteOutIntArray(d_geometry, 10, num_identifiers, filename, 0, 255);

		//filename = "D:\\temp\\d_seeds.png";
		//WriteOutIntArray(d_seeds, 2 * max_nx * max_ny, num_identifiers, filename, 0, 255);

		// Copy d_geometry to h_geometry and then loop through each superpixel in h_geometry and extract the seed values at offsets 8 and 9.
		// For any seed values that are not -1, -1, add those to add_seeds.
		cudaStatus = cudaMemcpy(h_geometry, d_geometry, 10 * sizeof(int) * num_identifiers, cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess)
		{
			throw std::runtime_error("Failed to copy h_geometry to host in c_split_seeds.\n");
			return add_seeds;
		}

		for (int i = 0; i < num_identifiers; ++i)
		{
			PointPair candidate_seed;
			candidate_seed.x = h_geometry[i * 10 + 8];
			candidate_seed.y = h_geometry[i * 10 + 9];
			if ((candidate_seed.x >= 0) && (candidate_seed.y >= 0))
			{
				add_seeds.push_back(candidate_seed);
			}
		}

		FreeIntArray(d_seeds);
		free(h_seeds);
		FreeIntArray(max_array);
		FreeIntArray(temp_pixdata);
		FreeIntArray(delta_spdata);
		FreeBoolArray(collision_data);
		FreeIntArray(regions_pixdata);
		cudaFree(d_identifiers);
		cudaFree(d_geometry);
		free(h_identifiers);
		free(h_geometry);
	}
	return add_seeds;
}
