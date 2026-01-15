// Copyright (c) 2023-2025 Steve Young
// Licensed under the MIT License

# include "Image_CUDA.cuh"

__global__ void process_dilate_erode_rect_1(int x, int y, unsigned char* input, unsigned char* output, bool isdilate, int struct_size)
{
	// This function is only for a rectangular structuring element.
	// This kernel takes in:
	// x, y - image width and height
	// input - image data in row-major format, one unsigned char per pixel.
	// output - image data for the transformed image, after undergoing dilation or erosion.
	// isdilate - true for dilation, false for erosion.
	// struct_size - horizontal and vertical size of the structuring element.
	// 
	// When calling this kernel, a third calling parameter is needed, which is the amount of shared
	// memory needed.  It requires struct_size*(blockDim.x + struct_size - 1) unsigned chars.
	// 
	// This version of the kernel was intended to be faster, by using more shared memory in order
	// to reduce the number of reads from global memory.  Tests show it to be only slightly
	// faster than the _2 version, and the simpler process of calling the _2 version (which does
	// not need careful checking for memory limitations) means it is probably best to just use
	// that one.
	int data_window_height = struct_size;
	int data_window_width = blockDim.x + struct_size - 1; // Width of data_window;
	int struct_extent = (struct_size - 1) / 2; // How far the struct goes in each direction beyond the base location.
	extern __shared__ unsigned char data_window[]; // The rolling window of shared data for calculations.
	int x_offset = blockIdx.x * blockDim.x;
	int img_x = x_offset + threadIdx.x; // X location within image.
	int wx = threadIdx.x + struct_extent; // X location within data window.
	int histogram[256]; // A per-thread histogram of values in the current window.
	int swap_y = 0;

	// First, read in data.
	// Start with main set of data.
	for (int wy = struct_extent; wy < data_window_height; ++wy) // wy is the y location within the data window.
	{
		int img_y = wy - struct_extent;
		if ((img_y >= 0) && (img_y < y))
		{
			// Stride across data_window, with w_x being the x location within the window.
			for (int w_x = threadIdx.x; w_x < data_window_width; w_x += blockDim.x)
			{
				int i_x = w_x - struct_extent + x_offset; // The corresponding x location in the image.
				if ((i_x >= 0) && (i_x < x))
				{
					data_window[wy * data_window_width + w_x] = input[i_x + img_y * x];
				}
			}
		}
	}

	// Now need to sync up, so all threads wait until the above is finished.
	__syncthreads();
	if ((img_x >= 0) && (img_x < x)) // Ensure this column is within bounds.
	{
		// Data window is populated for structuring elements at y=0.
		// Need to fill in histogram with these initial values.
		for (int i = 0; i < 256; ++i)
		{
			histogram[i] = 0;
		}
		for (int img_y = 0; img_y <= struct_extent; ++img_y)  // wy is the y value in data_window.
		{
			int wy = img_y + struct_extent;
			for (int i = wx - struct_extent; i <= wx + struct_extent; ++i) // i is the x value in data_window.
			{
				int ix = i - struct_extent + x_offset; // ix is in image-space.
				if ((ix >= 0) && (ix < x))
				{
					histogram[data_window[i + wy * data_window_width]] += 1;
				}
			}
		}
		// Write value for img_y=0.
		if (isdilate)
		{
			// Look for maximum value.
			for (int i = 255; i >= 0; --i)
			{
				if (histogram[i] > 0)
				{
					output[img_x] = i; // Simplified addressing, because img_y = 0.
					break;
				}
			}
		}
		else {
			// Look for minimum value.
			for (int i = 0; i <= 255; ++i)
			{
				if (histogram[i] > 0)
				{
					output[img_x] = i;
					break;
				}
			}
		}
	}

	// Begin marching forward with larger y values.
	for (int img_y = 1; img_y < y; ++img_y)
	{
		// swap_y will be the row that has the data to be removed, and then
		// will also be used to add in the new data.
		swap_y = (img_y - 1) % data_window_height;
		if ((img_x >= 0) && (img_x < x)) // Ensure this column is within bounds.
		{
			// First, adjust histogram based on new window.
			// Start by subtracting out the data that is no longer in scope.

			if (img_y > struct_extent)
			{
				for (int i = wx - struct_extent; i <= wx + struct_extent; ++i) // i is the x value in data_window.
				{
					int ix = i - struct_extent + x_offset; // ix is in image-space.
					if ((ix >= 0) && (ix < x))
					{
						histogram[data_window[i + swap_y * data_window_width]] -= 1;
					}
				}
			}
		}
		__syncthreads();
		int top_y = img_y + struct_extent;
		if (top_y < y)
		{
			// Now, replace window data at swap_y with new data from image.
			// Stride across data_window, with w_x being the x location within the window.
			for (int w_x = threadIdx.x; w_x < data_window_width; w_x += blockDim.x)
			{
				int i_x = w_x - struct_extent + x_offset; // The corresponding x location in the image.
				if ((i_x >= 0) && (i_x < x))
				{
					data_window[swap_y * data_window_width + w_x] = input[i_x + top_y * x];
				}
			}
		}
		// Sync threads, to make sure that histogram is updated with new data.
		__syncthreads();

		if ((img_x >= 0) && (img_x < x)) // Ensure this column is within bounds.
		{
			if (top_y < y)
			{
				// Now, update histogram by adding in the new data.
				for (int i = wx - struct_extent; i <= wx + struct_extent; ++i) // i is the x value in data_window.
				{
					int ix = i - struct_extent + x_offset; // ix is in image-space.
					if ((ix >= 0) && (ix < x))
					{
						histogram[data_window[i + swap_y * data_window_width]] += 1;
					}
				}
			}
			// Finally, write the output value for this img_y.
			if (isdilate)
			{
				// Look for maximum value.
				for (int i = 255; i >= 0; --i)
				{
					if (histogram[i] > 0)
					{
						output[img_x + img_y * x] = i;
						break;
					}
				}
			}
			else {
				// Look for minimum value.
				for (int i = 0; i <= 255; ++i)
				{
					if (histogram[i] > 0)
					{
						output[img_x + img_y * x] = i;
						break;
					}
				}
			}
		}
	}
}
__global__ void process_dilate_erode_rect_2(int x, int y, unsigned char* input, unsigned char* output, bool isdilate, int struct_size)
{
	// This function is only for a rectangular structuring element.
	// This kernel takes in:
	// x, y - image width and height
	// input - image data in row-major format, one unsigned char per pixel.
	// output - image data for the transformed image, after undergoing dilation or erosion.
	// isdilate - true for dilation, false for erosion.
	// struct_size - horizontal and vertical size of the structuring element.
	// 
	// When calling this kernel, a third calling parameter is needed, which is the amount of shared
	// memory needed.  It requires 2*(blockDim.x + struct_size - 1) unsigned chars.
	// 
	// The lower shared memory requirement of this version versus _1 makes it simpler to use.  It
	// performs almost as fast.

	int data_window_width = blockDim.x + struct_size - 1; // Width of data_window;
	int struct_extent = (struct_size - 1) / 2; // How far the struct goes in each direction beyond the base location.
	extern __shared__ unsigned char data_window[]; // Two rows of the rolling window of shared data for calculations. 0 for removal, 1 for adding.
	int x_offset = blockIdx.x * blockDim.x;
	int img_x = x_offset + threadIdx.x; // X location within image.
	int wx = threadIdx.x + struct_extent; // X location within data window.
	int histogram[256]; // A per-thread histogram of values in the current window.

	// Need to fill in histogram with initial values.
	for (int i = 0; i < 256; ++i)
	{
		histogram[i] = 0;
	}

	for (int img_y = 0; img_y <= struct_extent; ++img_y)  // wy is the y value in data_window.
	{
		// Stride across data_window, with w_x being the x location within the window.
		for (int w_x = threadIdx.x; w_x < data_window_width; w_x += blockDim.x)
		{
			int i_x = w_x - struct_extent + x_offset; // The corresponding x location in the image.
			if ((i_x >= 0) && (i_x < x))
			{
				data_window[data_window_width + w_x] = input[i_x + img_y * x];
			}
		}
		__syncthreads();
		if ((img_x >= 0) && (img_x < x))
		{
			for (int i = wx - struct_extent; i <= wx + struct_extent; ++i)
			{
				int ix = i - struct_extent + x_offset;
				if ((ix >= 0) && (ix < x)) {
					histogram[data_window[data_window_width + i]] += 1;
				}
			}
		}
		__syncthreads();
	}
	// Write value for img_y=0.
	if ((img_x >= 0) && (img_x < x))
	{
		if (isdilate)
		{
			// Look for maximum value.
			for (int i = 255; i >= 0; --i)
			{
				if (histogram[i] > 0)
				{
					output[img_x] = i; // Simplified addressing, because img_y = 0.
					break;
				}
			}
		}
		else {
			// Look for minimum value.
			for (int i = 0; i <= 255; ++i)
			{
				if (histogram[i] > 0)
				{
					output[img_x] = i;
					break;
				}
			}
		}
	}

	// Begin marching forward with larger y values.
	for (int img_y = 1; img_y < y; ++img_y)
	{
		bool out_flag = false; // Flag indicates whether to remove data from histogram.
		bool in_flag = false; // Flag indicates whether to add data to histogram.

		// First, fill out data_window.
		int out_y = img_y - struct_extent - 1;
		int in_y = img_y + struct_extent;
		if (out_y >= 0)
		{
			out_flag = true;
			// Stride across data_window, with w_x being the x location within the window.
			for (int w_x = threadIdx.x; w_x < data_window_width; w_x += blockDim.x)
			{
				int i_x = w_x - struct_extent + x_offset; // The corresponding x location in the image.
				if ((i_x >= 0) && (i_x < x))
				{
					data_window[w_x] = input[i_x + out_y * x];
				}
			}
		}
		if (in_y < y)
		{
			in_flag = true;
			// Stride across data_window, with w_x being the x location within the window.
			for (int w_x = threadIdx.x; w_x < data_window_width; w_x += blockDim.x)
			{
				int i_x = w_x - struct_extent + x_offset; // The corresponding x location in the image.
				if ((i_x >= 0) && (i_x < x))
				{
					data_window[data_window_width + w_x] = input[i_x + in_y * x];
				}
			}
		}
		__syncthreads();
		if ((img_x >= 0) && (img_x < x)) // Ensure this column is within bounds.
		{
			// Now, adjust histogram.
			if (out_flag)
			{
				for (int i = wx - struct_extent; i <= wx + struct_extent; ++i)
				{
					int ix = i - struct_extent + x_offset;
					if ((ix >= 0) && (ix < x)) {
						histogram[data_window[i]] -= 1;
					}
				}
			}
			if (in_flag)
			{
				for (int i = wx - struct_extent; i <= wx + struct_extent; ++i)
				{
					int ix = i - struct_extent + x_offset;
					if ((ix >= 0) && (ix < x)) {
						histogram[data_window[data_window_width + i]] += 1;
					}
				}
			}
			// Finally, write the output value for this img_y.
			if (isdilate)
			{
				// Look for maximum value.
				for (int i = 255; i >= 0; --i)
				{
					if (histogram[i] > 0)
					{
						output[img_x + img_y * x] = i;
						break;
					}
				}
			}
			else {
				// Look for minimum value.
				for (int i = 0; i <= 255; ++i)
				{
					if (histogram[i] > 0)
					{
						output[img_x + img_y * x] = i;
						break;
					}
				}
			}
		}
		__syncthreads();
	}
}

__global__ void process_dilate_erode_disc(int x, int y, unsigned char* input, unsigned char* output, bool isdilate, int struct_size)
{
	// This function is only for a disc structuring element.
	// This kernel takes in:
	// x, y - image width and height
	// input - image data in row-major format, one unsigned char per pixel.
	// output - image data for the transfored image, after undering dilation or erosion.
	// isdilate - true for dilation, false for erosion.
	// struct_size - horizontal and vertical size of the structing element.
	// 
	// When calling this kernel, a third calling parameter is needed, which is the amount of shared
	// memory needed.  It requires struct_extent + 1 ints.
	// 

	int struct_extent = (struct_size - 1) / 2; // How far the struct goes in each direction beyond the base location.
	extern __shared__ int disc_chord[]; // Defines shape of the disc.
	int x_offset = blockIdx.x * blockDim.x;
	int img_x = x_offset + threadIdx.x; // X location within image.
	int histogram[256]; // A per-thread histogram of values in the current window.
	int disc_index = 0;

	// Need to fill in histogram with initial values.
	for (int i = 0; i < 256; ++i)
	{
		histogram[i] = 0;
	}

	// Calculate disc chord.
	// Stride across half the extent (plus origin), with w_x being the x location within the window.
	for (int i = threadIdx.x; i <= struct_extent; i += blockDim.x)
	{
		//disc_chord[i] = struct_extent;
		disc_chord[i] = sqrt((float)(struct_extent * struct_extent - i * i));
	}
	__syncthreads();

	if ((img_x >= 0) && (img_x < x))
	{
		// Begin marching through y values.
		for (int img_y = -struct_extent; img_y < y; ++img_y)
		{
			for (int i = img_x - struct_extent; i <= img_x + struct_extent; ++i)
			{
				if ((i >= 0) && (i < x))
				{
					disc_index = abs(img_x - i);
					int j = img_y - disc_chord[disc_index] - 1;
					if ((j >= 0) && (j < y)) {
						histogram[input[i + j * x]] -= 1;
					}
					j = img_y + disc_chord[disc_index];
					if ((j >= 0) && (j < y)) {
						histogram[input[i + j * x]] += 1;
					}
				}
			}
			// Finally, write the output value for this img_y.
			if ((img_y >= 0) && (img_y < y))
			{
				if (isdilate)
				{
					// Look for maximum value.
					for (int i = 255; i >= 0; --i)
					{
						if (histogram[i] > 0)
						{
							output[img_x + img_y * x] = i;
							break;
						}
					}
				}
				else {
					// Look for minimum value.
					for (int i = 0; i <= 255; ++i)
					{
						if (histogram[i] > 0)
						{
							output[img_x + img_y * x] = i;
							break;
						}
					}
				}
			}
		}
	}
}
__global__ void process_gen_grayscale(unsigned char* c_device_data_input, int N, unsigned char* c_device_data_output)
{
	// Generic grayscale calculation for 3 channel data.
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < N)
	{
		float r = (float)c_device_data_input[3 * idx];
		float g = (float)c_device_data_input[3 * idx + 1];
		float b = (float)c_device_data_input[3 * idx + 2];
		c_device_data_output[idx] = (unsigned char)(sqrtf(r * r + g * g + b * b) * 0.57735);
	}
}

__global__ void process_gen_grayscale_c(unsigned char* c_device_data_input, int N, int c, unsigned char* c_device_data_output)
{
	// Grayscale calcuation using only one channel, indicated by c, where c is 1, 2, or 3.
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if ((c > 0) && (c < 4) && (idx < N))
	{
		c_device_data_output[idx] = c_device_data_input[3 * idx + c - 1];
	}
}

__global__ void process_gen_grayscale_c_nc(unsigned char* c_device_data_input, int N, int c, int nc, unsigned char* c_device_data_output)
{
	// Grayscale calculation with a positive channel (c) and a negative channel (nc).
	// Each of c and nc is 1, 2, or 3.
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if ((c > 0) && (c < 4) && (nc > 0) && (nc < 4) && (idx < N))
	{
		int combined = (int)c_device_data_input[3 * idx + c - 1];
		combined = combined - (int)c_device_data_input[3 * idx + nc - 1];
		c_device_data_output[idx] = (unsigned char)((combined + 255) / 2);
	}
}


unsigned char* c_gen_dilate_erode(int x, int y, unsigned char* h_in, bool isdilate, int mode, int struct_size) {
	// CUDA function to dilate or erode an image with dimensions x and y.
	// isdilate is true if this is dilation, false if this is erosion.
	// mode is 0 if the structuring element is rectangular, and 1 if it is a disc.
	// struct_size is the vertical and horizontal size of the structuring element (non-square sizes are not supported).
	// The host version of the input data is h_in, the device version is d_in.
	// The host version of the output data is h_out, the device version is d_out.

	unsigned char* d_in = NULL;
	unsigned char* d_out = NULL;
	unsigned char* h_out = NULL; // The return value from this function.
	int N = x * y;
	cudaError_t cudaStatus;

	d_in = UCharArray(x, y, false);
	d_out = UCharArray(x, y, false);
	h_out = (unsigned char*)malloc(sizeof(unsigned char) * N);
	if (NULL == h_out)
	{
		throw (std::runtime_error("Unable to allocate memory for image in c_gen_dilate_erode.\n"));
		return NULL;
	}
	if (false == CopyFromHost(h_in, N, d_in))
	{
		throw (std::runtime_error("Failed to transfer memory to device in c_gen_dilate_erode.\n"));
		FreeUCharArray(d_in);
		FreeUCharArray(d_out);
		free(h_out);
		return NULL;
	}

	int struct_extents = (struct_size - 1) / 2; // This if the extent to which the structuring element extends beyond the origin in each direction.
	struct_size = 1 + 2 * struct_extents; // Adjusting as needed if struct size is even.
	int blocks_per_grid = (x + threads_per_block - 1) / threads_per_block;
	if (0 == mode)
	{
		int shared_size = 2 * (struct_size + threads_per_block - 1);
		process_dilate_erode_rect_2 << <blocks_per_grid, threads_per_block, shared_size * sizeof(unsigned char) >> > (x, y, d_in, d_out, isdilate, struct_size);
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess)
		{
			throw std::runtime_error("Failed to process in process_dilate_erode_rect_2.\n");
			return NULL;
		}
	}
	else if (1 == mode)
	{
		int shared_size = struct_extents + 1;
		process_dilate_erode_disc << <blocks_per_grid, threads_per_block, shared_size * sizeof(int) >> > (x, y, d_in, d_out, isdilate, struct_size);
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess)
		{
			throw std::runtime_error("Failed to process in process_dilate_erode_disc.\n");
			return NULL;
		}
	}
	if (false == CopyToHost(d_out, N, h_out))
	{
		throw (std::runtime_error("Failed to transfer memory from device in c_gen_dilate_erode.\n"));
		FreeUCharArray(d_in);
		FreeUCharArray(d_out);
		free(h_out);
		return NULL;
	}
	FreeUCharArray(d_in);
	FreeUCharArray(d_out);
	return h_out;
}

bool c_gen_gray(unsigned char* c_device_data_input, int x, int y, int channels, int c, int nc, unsigned char* c_device_data_output)
{
	// Calculate a grayscale image from another image.
	// x,y are the width and height of the images.
	// channels is the number of color channels (either 1 or 3).
	// c is the positive color channel to use (if non-zero).
	// nc is the negative color channel to use (if non-zero).
	// If c and nc are given, then the grayscale image is the value of c-nc everywhere.
	// c_device_data_[input|output] arguments are the input and output data (as specified in the name).
	// The input data is x*y*channels bytes in length.  The output data is x*y bytes in length.
	int N = x * y;
	cudaError_t cudaStatus;
	int blocks_per_grid = (N + threads_per_block - 1) / threads_per_block;

	if (channels > 2) // Multi-channel image on input.
	{
		if (0 == c) // No specific color channels, so generic grayscale.
		{
			process_gen_grayscale << < blocks_per_grid, threads_per_block >> > (c_device_data_input, N, c_device_data_output);
			cudaDeviceSynchronize();
			cudaStatus = cudaGetLastError();
			if (cudaStatus != cudaSuccess)
			{
				throw std::runtime_error("Failed in process_gen_grayscale.\n");
				return false;
			}
		}
		else if ((0 < c) && (c <= 3) && (0 == nc)) {
			process_gen_grayscale_c << < blocks_per_grid, threads_per_block >> > (c_device_data_input, N, c, c_device_data_output);
			cudaDeviceSynchronize();
			cudaStatus = cudaGetLastError();
			if (cudaStatus != cudaSuccess)
			{
				throw std::runtime_error("Failed in process_gen_grayscale_c.\n");
				return false;
			}
		}
		else if ((0 < c) && (c <= 3) && (0 < nc) && (nc <= 3) && (c != nc)) // channel and nchannel both nonzero and in range.
		{
			process_gen_grayscale_c_nc << < blocks_per_grid, threads_per_block >> > (c_device_data_input, N, c, nc, c_device_data_output);
			cudaDeviceSynchronize();
			cudaStatus = cudaGetLastError();
			if (cudaStatus != cudaSuccess)
			{
				throw std::runtime_error("Failed in process_gen_grayscale_c_nc.\n");
				return false;
			}
		}
		else {
			throw std::runtime_error("Values for channel and nchannel out of range or conflict.\n");
			return false;
		}
	}
	else {
		// Calculate grayscale from a monochrome image.
		return OnDeviceCopy(c_device_data_input, N, c_device_data_output); // Just a copy of the data.
	}
	return true;
}
