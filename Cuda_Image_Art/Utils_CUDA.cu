// Copyright (c) 2023-2025 Steve Young
// Licensed under the MIT License

#include "Utils_CUDA.cuh"

__global__ void initialize_Float_Array(float* array, float value, int N)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < N)
	{
		array[idx] = value;
	}
}

__global__ void initialize_Int_Array(int* array, int value, int N)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < N)
	{
		array[idx] = value;
	}
}

__global__ void initialize_UChar_Array(unsigned char* array, unsigned char value, int N)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < N)
	{
		array[idx] = value;
	}
}

__global__ void renormalize_Float_Array(float* array, int N, float min_value, float max_value)
{
	float low_val = array[0];
	float high_val = array[0];
	for (int i = 1; i < N; ++i)
	{
		float val = array[i];
		if (val > high_val)
		{
			high_val = val;
		}
		if (val < low_val)
		{
			low_val = val;
		}
	}
	float old_range = high_val - low_val;
	float new_range = max_value - min_value;
	float scale = 1.0f;
	if (0.0f != old_range)
	{
		scale = new_range / old_range;
	}
	for (int i = 0; i < N; ++i)
	{
		array[i] = min_value + (array[i] - low_val) * scale;
	}
}

__global__ void initialize_Bool_Array(bool* array, bool value, int N)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < N)
	{
		array[idx] = value;
	}
}

__global__ void test_Bool_Array(bool* array, int w, int h, int wx, int wy, bool* result)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int N = chunk_size * chunk_size;
	if (idx < N)
	{
		int y = wy + idx % chunk_size;
		int x = wx + idx / chunk_size;
		if ((x < w) && (y < h))
		{
			int pos = x + w * y;
			if (array[pos])
			{
				*result = true;
			}
		}
	}
}

__global__ void process_Copy_Unsigned_Char_Array(unsigned char* source, int N, unsigned char* target)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < N)
	{
		target[idx] = source[idx];
	}
}

__global__ void process_Copy_Int_Array(int* source, int N, int* target)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < N)
	{
		target[idx] = source[idx];
	}
}

__global__ void process_add_partial_Float_Array(float* source, int source_width, int source_height, float* target, int target_width, int target_height, int target_x_offset, int target_y_offset)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int src_N = source_width * source_height;
	if (idx < src_N)
	{
		int dst_N = target_width * target_height;
		int src_y = idx / source_width;
		int src_x = idx % source_width;
		int dst_x = target_x_offset + src_x;
		int dst_y = target_y_offset + src_y;
		int dst_idx = dst_y * target_width + dst_x;
		if ((dst_x < target_width) && (dst_y < target_height))
		{
			target[dst_idx] += source[idx];
		}
	}
}

__global__ void process_copy_partial_Float_Array(float* source, int source_width, int source_height, float* target, int target_width, int target_height, int target_x_offset, int target_y_offset)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int src_N = source_width * source_height;
	if (idx < src_N)
	{
		int dst_N = target_width * target_height;
		int src_y = idx / source_width;
		int src_x = idx % source_width;
		int dst_x = target_x_offset + src_x;
		int dst_y = target_y_offset + src_y;
		int dst_idx = dst_y * target_width + dst_x;
		if (dst_idx < dst_N)
		{
			target[dst_idx] = source[idx];
		}
	}
}

__global__ void process_copy_Float_Array_portion(float* source, int source_image_width, int source_image_height, int source_x_offset, int source_y_offset, int portion_width, int portion_height, float* target, int target_width, int target_height, int target_x_offset, int target_y_offset)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int src_N = portion_width * portion_height;
	if (idx < src_N)
	{
		int dst_N = target_width * target_height;
		int src_y = idx / portion_width + source_y_offset;
		int src_x = idx % portion_width + source_x_offset;
		int src_idx = src_x + src_y * source_image_width;
		int dst_x = target_x_offset + src_x - source_x_offset;
		int dst_y = target_y_offset + src_y - source_y_offset;
		int dst_idx = dst_y * target_width + dst_x;
		if ((dst_idx < dst_N) && (src_idx < source_image_width * source_image_height))
		{
			target[dst_idx] = source[src_idx];
		}
	}
}

__global__ void process_add_Float_Array_portion(float* source, int source_image_width, int source_image_height, int source_x_offset, int source_y_offset, int portion_width, int portion_height, float* target, int target_width, int target_height, int target_x_offset, int target_y_offset)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int src_N = portion_width * portion_height;
	if (idx < src_N)
	{
		int dst_N = target_width * target_height;
		int src_y = idx / portion_width + source_y_offset;
		int src_x = idx % portion_width + source_x_offset;
		int src_idx = src_x + src_y * source_image_width;
		int dst_x = target_x_offset + src_x - source_x_offset;
		int dst_y = target_y_offset + src_y - source_y_offset;
		int dst_idx = dst_y * target_width + dst_x;
		if ((dst_idx < dst_N) && (src_idx < source_image_width * source_image_height))
		{
			target[dst_idx] += source[src_idx];
		}
	}
}

__global__ void initialize_gauss_kernel(float* matrix, int g_radius, float thin_factor)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int side_size = g_radius + 1;
	int N = side_size * side_size;
	if (idx < N)
	{
		int x = idx % side_size;
		int y = idx / side_size;
		matrix[idx] = expf(-sqrtf(x * x + y * y) * thin_factor);
	}
}

__global__ void set_single_value(float* matrix, int w, int h, int x, int y, float value, bool sum)
{
	if ((x >= 0) && (y >= 0) && (x < w) && (y < h))
	{
		int idx = x + y * w;
		if (sum)
		{
			matrix[idx] += value;
		}
		else {
			matrix[idx] = value;
		}
	}
}

float* GaussKernel(int gauss_radius, float gauss_thin_factor)
{
	// Gaussian smoothing.  
	// gauss_kernel[index] = exp(-sqrt(i * i + j * j) * gauss_thin_factor) where index = i + (gauss_radius+1)*j;
	// One quadrant represented of origin + gauss_radius in each dimension.

	float* ret = NULL;
	int N = (gauss_radius + 1) * (gauss_radius + 1);
	int blocksPerGrid = (N + threads_per_block - 1) / threads_per_block;
	cudaError_t cudaStatus;
	cudaStatus = cudaMalloc(&ret, sizeof(float) * N);
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed to allocate memory in GaussKernel.\n");
		return NULL;
	}
	initialize_gauss_kernel << <blocksPerGrid, threads_per_block >> > (ret, gauss_radius, gauss_thin_factor);
	cudaDeviceSynchronize();
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed to initialize matrix in GaussKernel.\n");
		return NULL;
	}
	return ret;
}

float GaussTotal(float* gauss_kernel, int gauss_radius)
{
	int N = (gauss_radius + 1) * (gauss_radius + 1);
	float* host_matrix = NULL;
	float gauss_total = 0.0;
	cudaError_t cudaStatus;

	host_matrix = (float*)malloc(N * sizeof(float));
	if (NULL == host_matrix)
	{
		throw std::runtime_error("Failed to allocate memory for host_matrix in GaussTotal.\n");
		return -1.0f;
	}
	cudaStatus = cudaMemcpy(host_matrix, gauss_kernel, N * sizeof(float), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess)
	{
		std::cout << "Error copying memory from device to host in GaussTotal.\n";
		return -1.0f;
	}
	int pos = 0;
	for (int j = 0; j <= gauss_radius; ++j)
	{
		for (int i = 0; i <= gauss_radius; ++i)
		{
			if ((i == 0) || (j == 0))
			{
				if ((i == 0) && (j == 0))
				{
					gauss_total += host_matrix[0];  // Special case at the origin.  There is only one instance of this.
				}
				else {
					gauss_total += 2.0 * host_matrix[pos]; // Case where either i or j (but not both) are on the axis.  Account for the instance directly across from the origin.
				}
			}
			else {
				gauss_total += 4.0 * host_matrix[pos]; // General case.  Account for all four quadrants.
			}
			pos++;
		}
	}
	if (gauss_total <= 0.0f)
	{
		gauss_total = 1.0f;
	}
	free(host_matrix);
	return gauss_total;
}

bool FreeGaussKernel(float* kernel)
{
	bool ret = true;
	if (NULL != kernel)
	{
		cudaFree(kernel);
	}
	return ret;
}

float* FloatArray(int x, int y, bool initialize_value, float value)
{
	// FloatArray is a 2D array stored in a 1D memory space.  It is laid out by rows (first row, then second, etc.).
	float* ret = NULL;
	int N = x * y;
	cudaError_t cudaStatus;
	cudaStatus = cudaMalloc(&ret, sizeof(float) * N);
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed to allocate memory in FloatArray.\n");
		return NULL;
	}
	if (initialize_value)
	{
		int blocksPerGrid = (N + threads_per_block - 1) / threads_per_block;
		initialize_Float_Array << <blocksPerGrid, threads_per_block >> > (ret, value, N);
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess)
		{
			throw std::runtime_error("Failed to initialize memory in FloatArray.\n");
			return NULL;
		}
	}
	return ret;
}

bool RenormalizeFloatArray(float* matrix, int x, int y, float min_value, float max_value)
{
	bool ret = true;
	int N = x * y;
	cudaError_t cudaStatus;
	renormalize_Float_Array << <1, 1 >> > (matrix, N, min_value, max_value);
	cudaDeviceSynchronize();
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed in RenormalizeFloatArray.\n");
		ret = false;
	}
	return ret;
}

bool* BoolArray(int x, int y, bool value)
{
	// BoolArray is a 2D array stored in a 1D memory space.  It is laid out by rows (first row, then second, etc.).
	bool* ret = NULL;
	int N = x * y;
	cudaError_t cudaStatus;
	cudaStatus = cudaMalloc(&ret, sizeof(bool) * N);
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed to allocate memory in BoolArray.\n");
		return NULL;
	}
	int blocksPerGrid = (N + threads_per_block - 1) / threads_per_block;
	initialize_Bool_Array << <blocksPerGrid, threads_per_block >> > (ret, value, N);
	cudaDeviceSynchronize();
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed to initialize memory in BoolArray.\n");
		return NULL;
	}
	return ret;
}

bool FreeFloatArray(float* a)
{
	bool ret = true;
	if (NULL == a)
	{
		ret = false;
	}
	else {
		cudaFree(a);
	}
	return ret;
}

bool FreeBoolArray(bool* a)
{
	bool ret = true;
	if (NULL == a)
	{
		ret = false;
	}
	else {
		cudaFree(a);
	}
	return ret;
}

bool ResetBoolArray(bool* matrix, int x, int y, bool value)
{
	bool ret = true;
	cudaError_t cudaStatus;
	int blocksPerGrid = (x * y + threads_per_block - 1) / threads_per_block;
	initialize_Bool_Array << <blocksPerGrid, threads_per_block >> > (matrix, value, x * y);
	cudaDeviceSynchronize();
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed to reset memory in ResetBoolArray.\n");
		ret = false;
	}
	return ret;
}

unsigned char* UCharArray(int x, int y, bool initialize_value, unsigned char value)
{
	unsigned char* ret = NULL;
	int N = x * y;
	cudaError_t cudaStatus;
	cudaStatus = cudaMalloc(&ret, sizeof(unsigned char) * N);
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed to allocate memory in UCharArray.\n");
		return NULL;
	}
	if (initialize_value)
	{
		int blocksPerGrid = (N + threads_per_block - 1) / threads_per_block;
		initialize_UChar_Array << <blocksPerGrid, threads_per_block >> > (ret, value, N);
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess)
		{
			throw std::runtime_error("Failed to initialize memory in UCharArray.\n");
			return NULL;
		}
	}
	return ret;
}

bool FreeUCharArray(unsigned char* a)
{
	bool ret = true;
	if (NULL == a)
	{
		ret = false;
	}
	else {
		cudaFree(a);
	}
	return ret;
}

bool ResetUCharArray(unsigned char* a, int w, int h, unsigned char value)
{
	bool ret = true;
	int N = w * h;
	cudaError_t cudaStatus;
	int blocksPerGrid = (N + threads_per_block - 1) / threads_per_block;
	initialize_UChar_Array << <blocksPerGrid, threads_per_block >> > (a, value, N);
	cudaDeviceSynchronize();
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed to initialize memory in ResetUCharArray.\n");
		ret = false;
	}
	return ret;
}

bool CopyFromHost(unsigned char* source, int N, unsigned char* dest)
{
	bool ret = true;
	cudaError_t cudaStatus;
	cudaStatus = cudaMemcpy(dest, source, N * sizeof(unsigned char), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed to copy memory to device in CopyFromHost.\n");
		ret = false;
	}
	return ret;
}

bool CopyToHost(unsigned char* source, int N, unsigned char* dest)
{
	bool ret = true;
	cudaError_t cudaStatus;
	cudaStatus = cudaMemcpy(dest, source, N * sizeof(unsigned char), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed to copy memory from device in CopyToHost.\n");
		ret = false;
	}
	return ret;
}

bool OnDeviceCopy(unsigned char* source, int N, unsigned char* dest)
{
	// Copy a matrix of unsigned char values from one on-device location to another.
	bool ret = true;
	cudaError_t cudaStatus;
	int blocksPerGrid = (N + threads_per_block - 1) / threads_per_block;
	process_Copy_Unsigned_Char_Array << <blocksPerGrid, threads_per_block >> > (source, N, dest);
	cudaDeviceSynchronize();
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed to copy memory in OnDeviceCopy.\n");
		ret = false;
	}
	return ret;
}

int* IntArray(int x, int y, bool initialize_value, int value)
{
	int* ret = NULL;
	int N = x * y;
	cudaError_t cudaStatus;
	cudaStatus = cudaMalloc(&ret, sizeof(int) * N);
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed to allocate memory in IntArray.\n");
		return NULL;
	}
	if (initialize_value)
	{
		int blocksPerGrid = (N + threads_per_block - 1) / threads_per_block;
		initialize_Int_Array << <blocksPerGrid, threads_per_block >> > (ret, value, N);
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess)
		{
			throw std::runtime_error("Failed to initialize memory in IntArray.\n");
			return NULL;
		}
	}
	return ret;
}

bool FreeIntArray(int* a)
{
	bool ret = true;
	if (NULL == a)
	{
		ret = false;
	}
	else {
		cudaFree(a);
	}
	return ret;
}

bool ResetIntArray(int* a, int w, int h, int value)
{
	bool ret = true;
	int N = w * h;
	cudaError_t cudaStatus;
	int blocksPerGrid = (N + threads_per_block - 1) / threads_per_block;
	initialize_Int_Array << <blocksPerGrid, threads_per_block >> > (a, value, N);
	cudaDeviceSynchronize();
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed to initialize memory in ResetIntArray.\n");
		ret = false;
	}
	return ret;
}

bool CopyFromHost(int* source, int N, int* dest)
{
	bool ret = true;
	cudaError_t cudaStatus;
	cudaStatus = cudaMemcpy(dest, source, N * sizeof(int), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed to copy memory to device in CopyFromHost.\n");
		ret = false;
	}
	return ret;
}

bool CopyToHost(int* source, int N, int* dest)
{
	bool ret = true;
	cudaError_t cudaStatus;
	cudaStatus = cudaMemcpy(dest, source, N * sizeof(int), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed to copy memory from device in CopyToHost.\n");
		ret = false;
	}
	return ret;
}

bool OnDeviceCopy(int* source, int N, int* dest)
{
	// Copy a matrix of unsigned char values from one on-device location to another.
	bool ret = true;
	cudaError_t cudaStatus;
	int blocksPerGrid = (N + threads_per_block - 1) / threads_per_block;
	process_Copy_Int_Array << <blocksPerGrid, threads_per_block >> > (source, N, dest);
	cudaDeviceSynchronize();
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed to copy memory in OnDeviceCopy.\n");
		ret = false;
	}
	return ret;
}

bool CopyFloatArray(float* source, float* target, int x, int y)
{
	bool ret = true;
	cudaError_t cudaStatus;
	if ((NULL == source) || (NULL == target))
	{
		throw std::runtime_error("NULL source or target passed to CopyFloatArray.\n");
		ret = false;
	}
	else {
		cudaStatus = cudaMemcpy(target, source, x * y * sizeof(float), cudaMemcpyDeviceToDevice);
		if (cudaStatus != cudaSuccess)
		{
			throw std::runtime_error("Failed to copy memory in CopyFloatArray.\n");
			ret = false;
		}
	}
	return ret;
}

bool AddPartialFloatArray(float* source, int source_width, int source_height, float* target, int target_width, int target_height, int target_x_offset, int target_y_offset)
{
	// Add full source array to sub-section of target array, at offset.
	bool ret = true;
	cudaError_t cudaStatus;
	int blocksPerGrid = (source_width * source_height + threads_per_block - 1) / threads_per_block;
	process_add_partial_Float_Array << <blocksPerGrid, threads_per_block >> > (source, source_width, source_height, target, target_width, target_height, target_x_offset, target_y_offset);
	cudaDeviceSynchronize();
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed to add memory in AddPartialFloatArray.\n");
		ret = false;
	}
	return ret;
}

bool CopyPartialFloatArray(float* source, int source_width, int source_height, float* target, int target_width, int target_height, int target_x_offset, int target_y_offset)
{
	// Copy full source array to sub-section of target array, at offset.
	bool ret = true;
	cudaError_t cudaStatus;
	int blocksPerGrid = (source_width * source_height + threads_per_block - 1) / threads_per_block;
	process_copy_partial_Float_Array << <blocksPerGrid, threads_per_block >> > (source, source_width, source_height, target, target_width, target_height, target_x_offset, target_y_offset);
	cudaDeviceSynchronize();
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed to copy memory in CopyPartialFloatArray.\n");
		ret = false;
	}
	return ret;
}

bool CopyFloatArrayPortion(float* source, int source_image_width, int source_image_height, int source_x_offset, int source_y_offset, int portion_width, int portion_height, float* target, int target_width, int target_height, int target_x_offset, int target_y_offset)
{
	// Copy a portion of the source array to a sub-section of target array, at offset.
	bool ret = true;
	cudaError_t cudaStatus;
	int blocksPerGrid = (portion_width * portion_height + threads_per_block - 1) / threads_per_block;
	process_copy_Float_Array_portion << <blocksPerGrid, threads_per_block >> > (source, source_image_width, source_image_height, source_x_offset, source_y_offset, portion_width, portion_height, target, target_width, target_height, target_x_offset, target_y_offset);
	cudaDeviceSynchronize();
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed to copy memory in CopyFloatArrayPortion.\n");
		ret = false;
	}
	return ret;
}

bool AddFloatArrayPortion(float* source, int source_image_width, int source_image_height, int source_x_offset, int source_y_offset, int portion_width, int portion_height, float* target, int target_width, int target_height, int target_x_offset, int target_y_offset)
{
	// Add a portion of the source array to a sub-section of target array, at offset.
	bool ret = true;
	cudaError_t cudaStatus;
	int blocksPerGrid = (portion_width * portion_height + threads_per_block - 1) / threads_per_block;
	process_add_Float_Array_portion << <blocksPerGrid, threads_per_block >> > (source, source_image_width, source_image_height, source_x_offset, source_y_offset, portion_width, portion_height, target, target_width, target_height, target_x_offset, target_y_offset);
	cudaDeviceSynchronize();
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed to copy memory in AddFloatArrayPortion.\n");
		ret = false;
	}
	return ret;
}

bool ResetFloatArray(float* matrix, int x, int y, float value)
{
	bool ret = true;
	cudaError_t cudaStatus;
	int blocksPerGrid = (x * y + threads_per_block - 1) / threads_per_block;
	initialize_Float_Array << <blocksPerGrid, threads_per_block >> > (matrix, value, x * y);
	cudaDeviceSynchronize();
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed to reset memory in ResetFloatArray.\n");
		ret = false;
	}
	return ret;
}

bool WriteOutFloatArray(float* source, int x, int y, std::string name, float min, float max)
{
	bool ret = true;
	cudaError_t cudaStatus;
	unsigned char* data;
	float* h_source;
	float delta = std::abs(max - min);
	data = (unsigned char*)malloc(sizeof(unsigned char) * x * y * 3);
	if (NULL == data)
	{
		throw std::runtime_error("Failed to allocate data for image.\n");
		ret = false;
	}
	if (ret)
	{
		h_source = (float*)malloc(sizeof(float) * x * y);
		if (NULL == h_source)
		{
			throw std::runtime_error("Failed to allocate host data for WriteOutFloatArray.\n");
			ret = false;
		}
	}
	if (ret)
	{
		cudaStatus = cudaMemcpy(h_source, source, x * y * sizeof(float), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess)
		{
			ret = false;
			std::cout << "Error copying memory from device to host in WriteOutFloatArray.\n";
		}
		long pos = 0;
		for (int j = 0; j < y; ++j)
		{
			for (int i = 0; i < x; ++i)
			{
				long image_pos = 3 * (j * x + i);
				unsigned char value;
				if (h_source[pos] > max)
				{
					value = 255;
				}
				else if (h_source[pos] < min)
				{
					value = 0;
				}
				else {
					value = 255 * ((h_source[pos] - min) / delta);
				}
				data[image_pos] = value;
				data[image_pos + 1] = value;
				data[image_pos + 2] = value;
				pos++;
			}
		}
		if (0 == stbi_write_png(name.c_str(), x, y, 3, data, x * 3))
		{
			throw std::runtime_error("Unable to write out matrix image in WriteOutFloatArray.\n");
			ret = false;
		}
		free(data);
		free(h_source);
	}
	return ret;
}

bool CopyFromHost(float* source, int N, float* dest)
{
	bool ret = true;
	cudaError_t cudaStatus;
	cudaStatus = cudaMemcpy(dest, source, N * sizeof(float), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed to copy memory to device in CopyFromHost.\n");
		ret = false;
	}
	return ret;
}

bool CopyToHost(float* source, int N, float* dest)
{
	bool ret = true;
	cudaError_t cudaStatus;
	cudaStatus = cudaMemcpy(dest, source, N * sizeof(float), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed to copy memory from device in CopyToHost.\n");
		ret = false;
	}
	return ret;
}

bool CopyFromHost(bool* source, int N, bool* dest)
{
	bool ret = true;
	cudaError_t cudaStatus;
	cudaStatus = cudaMemcpy(dest, source, N * sizeof(bool), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed to copy memory to device in CopyFromHost.\n");
		ret = false;
	}
	return ret;
}

bool CopyToHost(bool* source, int N, bool* dest)
{
	bool ret = true;
	cudaError_t cudaStatus;
	cudaStatus = cudaMemcpy((void*)dest, (const void*)source, N * sizeof(bool), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed to copy memory from device in CopyToHost.\n");
		ret = false;
	}
	return ret;
}

bool PartialCopyToHost(float* source, int source_width, int source_height, int source_x_offset, int source_y_offset, float* dest, int dest_width, int dest_height, int dest_x_offset, int dest_y_offset, int copy_width, int copy_height)
{
	// Copy all or part of device source array to a host target array.
	bool ret = true;
	cudaError_t cudaStatus;
	if ((NULL == source) || (NULL == dest))
	{
		throw std::runtime_error("NULL pointer passed to PartialCopyToHost.\n");
		ret = false;
	}

	copy_width = std::min(std::min(copy_width, dest_width - dest_x_offset), source_width - source_x_offset);
	copy_height = std::min(std::min(copy_height, dest_height - dest_y_offset), source_height - source_y_offset);

	for (int j = 0; ret && (j < copy_height); ++j)
	{
		int src_pos = (j + source_y_offset) * source_width + source_x_offset;
		int dest_pos = (j + dest_y_offset) * dest_width + dest_x_offset;

		cudaStatus = cudaMemcpy(&dest[dest_pos], &source[src_pos], copy_width * sizeof(float), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess)
		{
			throw std::runtime_error("Failed to copy memory from device in PartialCopyToHost.\n");
			ret = false;
		}
	}
	return ret;
}

bool WriteOutSparseFloatMatrix(SparseFloatMatrix* source, int x, int y, std::string name, float min, float max)
{
	bool ret = true;
	int x_chunks = source->GetXChunks();
	int y_chunks = source->GetYChunks();
	float background = source->GetBackround();
	int num_chunks = x_chunks * y_chunks;
	float* matrix = FloatArray(x, y, false);
	if (NULL == matrix)
	{
		throw std::runtime_error("Failed to allocate matrix for WriteOutSparseFloatMatrix.\n");
		ret = false;
	}
	float* background_array = FloatArray(chunk_size, chunk_size, true, background);
	if (NULL == background_array)
	{
		throw std::runtime_error("Failed to allocate background_array for WriteOutSparseFloatMatrix.\n");
		ret = false;
	}
	if (ret)
	{
		for (int chunk = 0; ret && (chunk < num_chunks); ++chunk)
		{
			int wx = chunk % x_chunks;
			int wy = chunk / x_chunks;
			wx = wx * chunk_size;
			wy = wy * chunk_size;
			if (source->CheckChunkNumber(chunk))
			{
				ret = CopyPartialFloatArray(source->GetChunk(chunk), chunk_size, chunk_size, matrix, x, y, wx, wy);
			}
			else {
				ret = CopyPartialFloatArray(background_array, chunk_size, chunk_size, matrix, x, y, wx, wy);
			}
		}
		ret = WriteOutFloatArray(matrix, x, y, name, min, max);
	}
	FreeFloatArray(background_array);
	FreeFloatArray(matrix);
	return ret;
}

SparseFloatMatrix::SparseFloatMatrix(int width, int height, float value)
{
	w = width;
	h = height;
	background_value = value;

	x_chunks = (w + chunk_size - 1) / chunk_size;
	y_chunks = (h + chunk_size - 1) / chunk_size;

	chunk_set.clear(); // Start with no allocated chunks.

	chunk = (float**)malloc(sizeof(float*) * x_chunks * y_chunks);
	if (NULL == chunk)
	{
		throw std::runtime_error("Failed to allocate chunk in SparseFloatMatrix.\n");
		return;
	}
	for (int chunk_count = 0; chunk_count < (x_chunks * y_chunks); ++chunk_count)
	{
		chunk[chunk_count] = NULL;
	}
	b_result = BoolArray(1, 1, false);
	d_chunk_updates = BoolArray(x_chunks, y_chunks, false);
	h_chunk_updates = (bool*)malloc(sizeof(bool) * x_chunks * y_chunks);
	if (NULL == h_chunk_updates)
	{
		throw std::runtime_error("Failed to allocate h_chunk_updates in SparseFloatMatrix.\n");
		return;
	}
}

SparseFloatMatrix::~SparseFloatMatrix()
{
	std::set<int>::iterator chunk_iterator;
	for (chunk_iterator = chunk_set.begin(); chunk_iterator != chunk_set.end(); ++chunk_iterator)
	{
		int chunk_index = *chunk_iterator;
		FreeFloatArray(chunk[chunk_index]);
	}
	free(chunk);
	FreeBoolArray(b_result);
	FreeBoolArray(d_chunk_updates);
	free(h_chunk_updates);
}

bool SparseFloatMatrix::Set_Value(int x, int y, float value, bool add)
{
	// Sets value in matrix, unless add is true in which case the value is added to existing value in the matrix.
	bool ret = true;
	cudaError_t cudaStatus;
	int chunk_index = GetChunkNumber(x, y);
	if (chunk_index >= 0)
	{
		int x_offset, y_offset, offset;
		if (false == CheckChunkNumber(chunk_index))
		{
			chunk[chunk_index] = FloatArray(chunk_size, chunk_size, true, background_value);
			chunk_set.insert(chunk_index);
		}
		x_offset = x % chunk_size;
		y_offset = y % chunk_size;
		offset = x_offset + y_offset * chunk_size;
		float device_value;
		cudaStatus = cudaMemcpy(&device_value, &chunk[chunk_index][offset], sizeof(float), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess)
		{
			ret = false;
			std::cout << "Error copying memory from device to host in SparseFloatMatrix::Set_Value.\n";
		}
		if (add)
		{
			device_value += value;
		}
		else {
			device_value = value;
		}
		cudaStatus = cudaMemcpy(&chunk[chunk_index][offset], &device_value, sizeof(float), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess)
		{
			ret = false;
			std::cout << "Error copying memory from host to device in SparseFloatMatrix::Set_Value.\n";
		}
	}
	else {
		throw std::runtime_error("Set_Value called with thick_width or y out of bounds.\n");
		ret = false;
	}
	return ret;
}

bool SparseFloatMatrix::Reset_Values(bool free_memory)
{
	// Reset to background, by removing allocated chunks.
	// free_memory indicates whether to fully free the memory (because it won't likely be needed again) or whether the values should just be reset to background_value.
	bool ret = true;
	std::set<int>::iterator chunk_it;
	if (free_memory)
	{
		for (chunk_it = chunk_set.begin(); chunk_it != chunk_set.end(); ++chunk_it)
		{
			int chunk_index = *chunk_it;
			ret = ret && FreeFloatArray(chunk[chunk_index]);
			chunk[chunk_index] = NULL;
		}
		chunk_set.clear();
	}
	else {
		for (chunk_it = chunk_set.begin(); chunk_it != chunk_set.end(); ++chunk_it)
		{
			int chunk_index = *chunk_it;
			ret = ret && ResetFloatArray(chunk[chunk_index], chunk_size, chunk_size, background_value);
		}
	}
	return ret;
}

float SparseFloatMatrix::Get_Value(int x, int y)
{
	// Return the value at x, y.
	float ret = 0.0f;
	cudaError_t cudaStatus;
	int chunk_index = GetChunkNumber(x, y);
	if ((false == CheckChunkNumber(chunk_index)) || (x < 0) || (x >= w) || (y < 0) || (y >= h))
	{
		ret = background_value;
	}
	else {
		int x_offset = x % chunk_size;
		int y_offset = y % chunk_size;
		int offset = x_offset + y_offset * chunk_size;

		cudaStatus = cudaMemcpy(&ret, &chunk[chunk_index][offset], sizeof(float), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess)
		{
			ret = false;
			std::cout << "Error copying memory from device to host in SparseFloatMatrix::Get_Value.\n";
		}
	}
	return ret;
}

bool SparseFloatMatrix::Copy(SparseFloatMatrix* src)
{
	// Copy a source SparseFloatMatrix to this one.
	bool ret = true;
	if ((src->w != w) || (src->h != h))
	{
		throw std::runtime_error("Copy called on SparseFloatMatrix with different size.\n");
		ret = false;
	}
	else {
		//  Go through the src->chunk_set, create chunks needed, and then copy from src chunks to local chunks.  Remove extra chunks.
		for (int chunk_index = 0; chunk_index < (x_chunks * y_chunks); ++chunk_index)
		{
			bool src_chunk, dest_chunk;
			src_chunk = src->CheckChunkNumber(chunk_index);
			dest_chunk = CheckChunkNumber(chunk_index);
			if (src_chunk) // The chunk is present in the source matrix.
			{
				if (!dest_chunk)
				{
					chunk[chunk_index] = FloatArray(chunk_size, chunk_size, true, background_value);
					chunk_set.insert(chunk_index);
				}
				ret = ret && CopyFloatArray(src->chunk[chunk_index], chunk[chunk_index], chunk_size, chunk_size);
			}
			else { // The chunk is not present in the source matrix.
				if (dest_chunk)
				{
					FreeFloatArray(chunk[chunk_index]);
					chunk_set.erase(chunk_index);
				}
			}
		}
	}
	background_value = src->background_value;
	return ret;
}

int SparseFloatMatrix::GetChunkNumber(int x, int y)
{
	int ret = 0;
	if ((x < 0) || (x >= w) || (y < 0) || (y >= h))
	{
		ret = -1;
	}
	else {
		ret = (x / chunk_size) + (y / chunk_size) * x_chunks;
	}
	return ret;
}

bool SparseFloatMatrix::CheckChunkNumber(int chunk_index)
{
	bool ret = (chunk_set.find(chunk_index) != chunk_set.end());
	return ret;
}

bool SparseFloatMatrix::Check(int x, int y)
{
	bool ret = false;
	if ((x >= 0) && (x < w) && (y >= 0) && (y < h))
	{
		int chunk_index = (x / chunk_size) + (y / chunk_size) * x_chunks;
		ret = (chunk_set.find(chunk_index) != chunk_set.end());
	}
	return ret;
}

std::set<int> SparseFloatMatrix::GetChunkSet()
{
	return chunk_set;
}

int SparseFloatMatrix::GetXChunks()
{
	return x_chunks;
}

int SparseFloatMatrix::GetYChunks()
{
	return y_chunks;
}

int SparseFloatMatrix::GetWidth()
{
	return w;
}

int SparseFloatMatrix::GetHeight()
{
	return h;
}

float* SparseFloatMatrix::GetChunk(int chunk_index)
{
	return chunk[chunk_index];
}

bool SparseFloatMatrix::TestMask(bool* M, int chunk_index)
{
	// Return true if any element of M is true within the identified chunk.
	cudaError_t cudaStatus;
	int blocksPerGrid = (chunk_size * chunk_size + threads_per_block - 1) / threads_per_block;
	bool ret = false;
	cudaStatus = cudaMemcpy(b_result, &ret, sizeof(bool), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess)
	{
		ret = false;
		std::cout << "Error copying memory from host to device in SparseFloatMatrix::TestMask.\n";
	}
	int wy = chunk_index / x_chunks;  // Y component of beginning of window represented by chunk.
	int wx = chunk_index % x_chunks; // X component of the beginning of windows represented by chunk. 
	wx = wx * chunk_size;
	wy = wy * chunk_size;

	test_Bool_Array << <blocksPerGrid, threads_per_block >> > (M, w, h, wx, wy, &ret);  // *** Shouldn't the passed in variable be b_result? ***

	cudaStatus = cudaMemcpy(&ret, b_result, sizeof(bool), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess)
	{
		ret = false;
		std::cout << "Error copying memory from host to device in SparseFloatMatrix::TestMask.\n";
	}
	return ret;
}

float SparseFloatMatrix::GetBackround()
{
	return background_value;
}

bool* SparseFloatMatrix::GetDChunkUpdates()
{
	return d_chunk_updates;
}

bool* SparseFloatMatrix::GetHChunkUpdates()
{
	return h_chunk_updates;
}

bool SparseFloatMatrix::UpdateChunks()
{
	bool ret = true;
	cudaError_t cudaStatus;
	// First, need to bring the vector of potential updates to the host.
	cudaStatus = cudaMemcpy(h_chunk_updates, d_chunk_updates, sizeof(bool) * x_chunks * y_chunks, cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess)
	{
		ret = false;
		std::cout << "Error copying memory from device to host in SparseFloatMatrix::UpdateChunks.\n";
	}
	else {
		for (int chunk_num = 0; chunk_num < (x_chunks * y_chunks); ++chunk_num)
		{
			if (h_chunk_updates[chunk_num]) // This chunk was identified as one that should be allocated.
			{
				if (!CheckChunkNumber(chunk_num))  // This chunk was identified as not yet allocated.
				{
					chunk[chunk_num] = FloatArray(chunk_size, chunk_size, true, background_value);
					chunk_set.insert(chunk_num);
				}
			}
		}
	}
	return ret;
}

bool SparseFloatMatrix::ExpandToFloatArray(float* d_dst)
{
	bool ret = true;
	int N = w * h;
	// First, set array to background values.
	ret = ResetFloatArray(d_dst, w, h, background_value);
	if (!ret)
	{
		throw std::runtime_error("Failed to reset background values in SparseFloatMatrix::ExpandToFloatArray.\n");
	}
	std::set<int>::iterator k;
	for (k = chunk_set.begin(); (ret && (k != chunk_set.end())); ++k)  // Go through allocated chunks and copy data to d_dst.
	{
		int chunk_num = *k;
		int x_offset = chunk_size * (chunk_num % x_chunks);
		int y_offset = chunk_size * (chunk_num / x_chunks);
		ret = ret && CopyPartialFloatArray(chunk[chunk_num], chunk_size, chunk_size, d_dst, w, h, x_offset, y_offset);
	}
	return ret;
}

bool SparseFloatMatrix::CompressFromFloatArray(float* d_src)
{
	bool ret = true;
	cudaError_t cudaStatus;
	// First, need to bring the vector of potential updates to the host.
	cudaStatus = cudaMemcpy(h_chunk_updates, d_chunk_updates, sizeof(bool) * x_chunks * y_chunks, cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess)
	{
		ret = false;
		std::cout << "Error copying memory from device to host in SparseFloatMatrix::CompressFromFloatArray.\n";
	}
	else {
		for (int chunk_num = 0; chunk_num < (x_chunks * y_chunks); ++chunk_num)
		{
			if (h_chunk_updates[chunk_num]) // This chunk was identified as one that was updated.
			{
				if (!CheckChunkNumber(chunk_num))  // This chunk was identified as not yet allocated.
				{
					chunk[chunk_num] = FloatArray(chunk_size, chunk_size, true, background_value);
					chunk_set.insert(chunk_num);
				}
				int source_x_offset = chunk_size * (chunk_num % x_chunks);
				int source_y_offset = chunk_size * (chunk_num / x_chunks);
				ret = ret && CopyFloatArrayPortion(d_src, w, h, source_x_offset, source_y_offset, chunk_size, chunk_size, chunk[chunk_num], chunk_size, chunk_size, 0, 0);
			}
		}
	}
	return ret;
}

bool SparseFloatMatrix::SyncChunks(SparseFloatMatrix* src)
{
	bool ret = true;
	std::set<int> src_set = src->GetChunkSet();
	for (std::set<int>::iterator chunk_it = src_set.begin(); chunk_it != src_set.end(); ++chunk_it)
	{
		int chunk_num = *chunk_it;
		if (!CheckChunkNumber(chunk_num))  // This chunk was identified as not yet allocated.
		{
			chunk[chunk_num] = FloatArray(chunk_size, chunk_size, true, background_value);
			chunk_set.insert(chunk_num);
		}
	}
	return ret;
}
