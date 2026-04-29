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
				if (level == 0)
				{
					std::string results_filename = "ResultsMap_" + std::to_string(level) + "_" + std::to_string(count) + ".png";
					WriteOutBoolArray(ResultsArray, blocksPerGrid, 1, results_filename);
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