// Copyright (c) 2023-2025 Steve Young
// Licensed under the MIT License

#include "..\General.h"
#include "..\Paper.h"
#include "Brush_CUDA.cuh"
//#include "Brush.h"
#include "..\SPixelData.h"
#include <iostream>

__global__ void process_relax_divergence_step1(bool* M, float* u, float* v, float* delta_matrix, int w, int h)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int N = w * h;
	if (idx < N)
	{
		//if (M[idx])
		//{
			int img_x = idx % w;
			int img_y = idx / w;
			int u_pos = img_y * (w + 1) + img_x;  // Position for u is affected by the fact that its width is w + 1.
			int v_pos = idx; // Position for v is not affected by height being h + 1;
			delta_matrix[idx] = -xi * (u[u_pos + 1] - u[u_pos] + v[v_pos + w] - v[v_pos]);
		//}
		//else {
		//	delta_matrix[idx] = 0.0f;
		//}
	}
}

__global__ void process_relax_divergence_step2(bool* M, float* u, float* v, float* p, float* delta_matrix, int w, int h)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int N = w * h;
	if (idx < N)
	{
		//if (M[idx])
		//{
			int img_x = idx % w;
			int img_y = idx / w;
			int u_pos = img_y * (w + 1) + img_x;  // Position for u is affected by the fact that its width is w + 1.
			int v_pos = idx; // Position for v is not affected by height being h + 1;
			float delta = delta_matrix[idx];
			float delta_left;
			float delta_top;
			if (img_x > 0)
			{
				delta_left = delta_matrix[idx - 1];
			}
			else {
				delta_left = 0.0f;
			}
			if (img_y > 0)
			{
				delta_top = delta_matrix[idx - w];
			}
			else {
				delta_top = 0.0f;
			}
			p[idx] += delta;

			u[u_pos] += (delta_left - delta);
			v[v_pos] += (delta_top - delta);

			if (img_x == w - 1)
			{
				u[u_pos + 1] += delta;
			}
			if (img_y == h - 1)
			{
				v[v_pos + w] += delta;
			}
		//}
	}
}

__global__ void process_relax_divergence_step2_b(bool* M, float* u, float* v, float* p, float* delta_matrix, int w, int h)
{
	// Block is block_dimension x block_dimension in size (currently 16 x 16).
	// threadIdx.x is interpreted to be the rows of the block.  0-15 is the first row, 16-31 is the second, etc.
	// Therefore, the correspondence to idx and the img_x and img_y will be more complicated.

	int img_x = blockIdx.x * blockDim.x + threadIdx.x;
	int img_y = blockIdx.y * blockDim.y + threadIdx.y;
	int img_idx = img_x + img_y * w;
	bool in_bounds = false;
	__shared__ float delta_shared[block_dimension][block_dimension];
	if ((img_x < w) && (img_y < h))
	{
		delta_shared[threadIdx.x][threadIdx.y] = delta_matrix[img_idx];
		in_bounds = true;
	}
	else {
		delta_shared[threadIdx.x][threadIdx.y] = 0.0;

	}
	__syncthreads();

	if (in_bounds)
	{
		if (M[img_idx])
		{
			int u_pos = img_y * (w + 1) + img_x;  // Position for u is affected by the fact that its width is w + 1.
			int v_pos = img_idx; // Position for v is not affected by height being h + 1;
			float delta = delta_shared[threadIdx.x][threadIdx.y];
			float delta_left;
			float delta_top;
			// Handle horizontal values first.
			if (threadIdx.x > 0)
			{
				delta_left = delta_shared[threadIdx.x - 1][threadIdx.y];
			}
			else {
				if (img_x > 0)
				{
					delta_left = delta_matrix[img_idx - 1];
				}
				else {
					delta_left = 0.0f;
				}
			}

			// Handle vertical values.
			if (threadIdx.y > 0)
			{
				delta_top = delta_shared[threadIdx.x][threadIdx.y - 1];
			}
			else {
				if (img_y > 0)
				{
					delta_top = delta_matrix[img_idx - w];
				}
				else {
					delta_top = 0.0f;
				}
			}

			p[img_idx] += delta;

			u[u_pos] += (delta_left - delta);
			v[v_pos] += (delta_top - delta);

			if (img_x == w - 1)
			{
				u[u_pos + 1] += delta;
			}
			if (img_y == h - 1)
			{
				v[v_pos + w] += delta;
			}
		}
	}
}

__global__ void process_convolution(float* C_input, float* C_output, float* C_kernel, float scale, int width, int height, int padding, int kernel_radius)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int N = width * height;
	if (idx < N)
	{
		int x = width + 2 * padding;
		int kernel_size = kernel_radius + 1;
		int K_N = kernel_size * kernel_size;  // Only one quadrant of the kernel is included, since all are reflected images of it.
		int top_pad = padding * x;  // The offset into C_output due to the padding of columns on the top side.
		int image_y = idx / width;
		int image_x = idx % width;
		int adusted_i = top_pad + x * image_y + padding + image_x;
		C_output[adusted_i] = 0.0f;
		// Move sequentially through the kernel.
		int K_x = 0;
		int K_y = 0;
		int I1 = adusted_i;
		int I2 = adusted_i;
		int I3 = adusted_i;
		int I4 = adusted_i;
		for (int K_i = 0; K_i < K_N; ++K_i)
		{
			// First, address (positive, positive) quadrant, using I1 for index into input image.
			C_output[adusted_i] += C_input[I1] * C_kernel[K_i];

			// Next, address (positive, negative) quadrant, using I2.  Note that K_x = 0 and K_y > 0 is included here.
			if (K_y > 0)
			{
				C_output[adusted_i] += C_input[I2] * C_kernel[K_i];
			}

			// Next address (negative, positive) quadrant, using I3.  Note that K_x > 0 and K_y = 0 is included here.
			if (K_x > 0)
			{
				C_output[adusted_i] += C_input[I3] * C_kernel[K_i];
			}

			// Lastly, address (negative, negative) quadrant, using I4.
			if ((K_x > 0) && (K_y > 0))
			{
				C_output[adusted_i] += C_input[I4] * C_kernel[K_i];
			}

			if (K_x < kernel_radius)
			{
				K_x++;
				I1++;
				I2++;
				I3--;
				I4--;
			}
			else {
				K_x = 0;
				K_y++;
				I1 = I1 - kernel_radius + x;
				I2 = I2 - kernel_radius - x;
				I3 = I3 + kernel_radius + x;
				I4 = I4 + kernel_radius - x;
			}
		}
		C_output[adusted_i] = C_output[adusted_i] / scale;
	}
}

__global__ void process_slope_velocities(float* thickness, float* u, float* v, int w, int h, int padding)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int N = w * h;
	if (idx < N)
	{
		int thick_width = w + 2 * padding;
		int top_pad = padding * thick_width;  // The offset into thickness due to the padding of columns on the top side.

		int img_x = idx % w;
		int img_y = idx / w;
		int u_pos = img_y * (w + 1) + img_x;  // Position for u is affected by the fact that its width is w + 1.
		int thick_pos = top_pad + thick_width * img_y + padding + img_x;  // Position for thickness is calculated based on padding.
		if (img_x > 0)
		{
			u[u_pos] -= slope_factor * (thickness[thick_pos] - thickness[thick_pos - 1]);
			if (w - 1 == img_x) // Pick up right side of last column.
			{
				u[u_pos + 1] -= slope_factor * (thickness[thick_pos + 1] - thickness[thick_pos]);
			}
		}
		if (img_y > 0)
		{
			v[idx] -= slope_factor * (thickness[thick_pos] - thickness[thick_pos - thick_width]);
			if (h - 1 == img_y) // Pick up bottom of last row.
			{
				v[idx + w] -= slope_factor * (thickness[thick_pos + thick_width] - thickness[thick_pos]);
			}
		}
	}
}

__global__ void process_max_velocity(float* u, float* v, float* results, int w, int h)
{
	int index = threadIdx.x;
	int stride = threads_per_block;
	int N = w * h;
	float value = -FLT_MAX;
	results[index] = -FLT_MAX;
	for (int i = index; i < N; i += stride)
	{
		int x = i % w;
		int y = i / w;

		int u_i = x + y * (w + 1);
		value = u[u_i];
		if (value > results[index])
		{
			results[index] = value;
		}
		if (w - 1 == x)
		{
			value = u[u_i + 1];
			if (value > results[index])
			{
				results[index] = value;
			}
		}
		value = v[i];
		if (value > results[index])
		{
			results[index] = value;
		}
		if (h - 1 == y)
		{
			value = v[i + w];
			if (value > results[index])
			{
				results[index] = value;
			}
		}
	}
}

__global__ void process_max_velocity_2(int w, int h, float* u, float* v, float* result)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int N = w * h;
	__shared__ float val[threads_per_block];
	if (idx < N)
	{
		float u_val = u[idx];
		float v_val = v[idx];
		if (u_val > v_val) {
			val[threadIdx.x] = u_val;
		}
		else {
			val[threadIdx.x] = v_val;
		}
	}
	else {
		val[threadIdx.x] = -FLT_MAX;
	}
	__syncthreads();

	for (unsigned int stride = blockDim.x / 2; stride > 0; stride >>= 1) {
		if ((threadIdx.x < stride) && (idx + stride < N))
		{
			if (val[threadIdx.x + stride] > val[threadIdx.x])
			{
				val[threadIdx.x] = val[threadIdx.x + stride];
			}
		}
		__syncthreads();
	}
	if (0 == threadIdx.x)
	{
		result[blockIdx.x] = val[0];
	}
}

__global__ void get_max_results(float* results)
{
	// Put max into first position in results.
	float max = results[0];
	for (int i = 1; i < threads_per_block; ++i)
	{
		if (results[i] > max)
		{
			max = results[i];
		}
	}
	results[0] = max;
}

__global__ void get_max_results_2(float* input, int input_length, float* output)
{
	// Put max of input into first position in output.
	// Position of the max of each input block is into the output position of that block number.

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	__shared__ float val[threads_per_block];
	if (idx < input_length)
	{
		val[threadIdx.x] = input[idx];
	}
	else {
		val[threadIdx.x] = -FLT_MAX;
	}
	__syncthreads();

	for (unsigned int stride = blockDim.x / 2; stride > 0; stride >>= 1) {
		if ((threadIdx.x < stride) && (idx + stride < input_length))
		{
			if (val[threadIdx.x + stride] > val[threadIdx.x])
			{
				val[threadIdx.x] = val[threadIdx.x + stride];
			}
		}
		__syncthreads();
	}
	if (0 == threadIdx.x)
	{
		output[blockIdx.x] = val[0];
	}
}

__global__ void process_calc_velocities(float dt, int w, int h, float* u, float* v, float* u_prime, float* v_prime, float* p)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int N = w * h;
	if (idx < N)
	{
		float A, B;
		float u_ij;  // Value of u at center of i, j.
		float u_i1j;  // Value of u at center of i+1, j.
		float v_ij; // Value of v at center of i, j.
		float v_ij1; // Value of v at center of i, j+1.

		int img_x = idx % w;
		int img_y = idx / w;
		int u_pos = img_y * (w + 1) + img_x;  // Position for u is affected by the fact that its width is w + 1.
		int v_pos = idx; // Position for v is not affected by height being h + 1;

		// First, do calculations for u.

		u_ij = (u[u_pos] + u[u_pos + 1]) / 2.0f;
		if (img_x < w - 1)
		{
			u_i1j = (u[u_pos + 1] + u[u_pos + 2]) / 2.0f;
		}
		else {
			u_i1j = 0.0f;  // *** Should this be u[u_pos+1]/2.0f ? ***
		}

		float u_val, v_val; // Temporary variables.
		float uv_corner_10; // Value of u*v at the i+.5, j-.5 corner.  First digit is 1 if x is advanced, second is 1 if y is advanced.
		if (img_y == 0) // First row
		{
			u_val = u[u_pos + 1];
		}
		else {
			u_val = (u[u_pos - w] + u[u_pos + 1]); // (u_pos + 1) - (w + 1) , (u_pos + 1)
		}
		if (img_x == w - 1) // Last column
		{
			v_val = v[v_pos];
		}
		else {
			v_val = v[v_pos] + v[v_pos + 1];
		}
		uv_corner_10 = u_val * v_val / 4.0f;


		float uv_corner_11; // Value of u*v at the i+.5, j+.5 corner.
		if (img_y == h - 1) // Last row
		{
			u_val = u[u_pos + 1];
		}
		else {
			u_val = u[u_pos + 1] + u[u_pos + w + 2]; // (u_pos + 1), (u_pos + 1) + (w + 1)
		}
		if (img_x == w - 1) // Last column
		{
			v_val = v[v_pos + w];
		}
		else {
			v_val = v[v_pos + w] + v[v_pos + 1 + w]; // (v_pos) + (w), (v_pos + 1) + (w)
		}
		uv_corner_11 = u_val * v_val / 4.0f;

		A = u_ij * u_ij - u_i1j * u_i1j + uv_corner_10 - uv_corner_11;
		B = u[u_pos] - 4.0 * u[u_pos + 1];
		if (img_x < (w - 1))
		{
			B += u[u_pos + 2];
		}
		if (img_y < (h - 1))
		{
			B += u[u_pos + w + 2]; // (u_pos + 1) + (w + 1)
		}
		if (img_y > 0)
		{
			B += u[u_pos - w]; // (u_pos + 1) - (w + 1)
		}

		if (img_x < (w - 1))
		{
			u_prime[u_pos + 1] = u[u_pos + 1] + dt * (A + mu * B + p[idx] - p[idx + 1] - kappa * u[u_pos + 1]);
		}
		else {
			u_prime[u_pos + 1] = u[u_pos + 1] + dt * (A + mu * B + p[idx] - kappa * u[u_pos + 1]);  // Drop out the p terms that rely on information outside the bounds.
		}
		if (u_prime[u_pos + 1] > max_velocity)
		{
			u_prime[u_pos + 1] = max_velocity;
		}
		else if (u_prime[u_pos + 1] < -max_velocity)
		{
			u_prime[u_pos + 1] = -max_velocity;
		}

		// Next, do calculations for v.
		v_ij = (v[idx] + v[idx + w]) / 2.0f;
		if (img_y < h - 1)
		{
			v_ij1 = (v[idx + w] + v[idx + 2 * w]) / 2.0f;
		}
		else {
			v_ij1 = 0.0f;  // *** Consider whether this should be v[i + w]/2.0f.
		}
		float uv_corner_01; // Value of u*v at the i-.5, j+.5 corner.
		if (img_y == h - 1) // Last row.
		{
			u_val = u[u_pos];
		}
		else {
			u_val = (u[u_pos] + u[u_pos + w + 1]);
		}
		if (img_x == 0)
		{
			v_val = v[idx + w];
		}
		else {
			v_val = (v[idx - 1 + w] + v[idx + w]); // (i - 1) + (w), (i) + (w)
		}
		uv_corner_01 = u_val * v_val / 4.0f;

		A = v_ij * v_ij - v_ij1 * v_ij1 + uv_corner_01 - uv_corner_11;
		B = v[idx] - 4.0 * v[idx + w];
		if (img_x < (w - 1))
		{
			B += v[idx + 1 + w];
		}
		if (img_x > 0)
		{
			B += v[idx - 1 + w];
		}
		if (img_y < (h - 1))
		{
			B += v[idx + 2 * w];
		}

		if (img_y < (h - 1))
		{
			v_prime[idx + w] = v[idx + w] + dt * (A + mu * B + p[idx] - p[idx + w] - kappa * v[idx + w]);
		}
		else {
			v_prime[idx + w] = v[idx + w] + dt * (A + mu * B + p[idx] - kappa * v[idx + w]);  // Drop out the p terms that rely on information outside the bounds.
		}
		if (v_prime[idx + w] > max_velocity)
		{
			v_prime[idx + w] = max_velocity;
		}
		else if (v_prime[idx + w] < -max_velocity)
		{
			v_prime[idx + w] = -max_velocity;
		}
	}
}

__global__ void process_enforce_boundaries(int w, int h, float* u, float* v, bool* M)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int N = w * h;
	if (idx < N)
	{
		int img_x = idx % w;
		int img_y = idx / w;
		int u_pos = img_y * (w + 1) + img_x;  // Position for u is affected by the fact that its width is w + 1.
		int v_pos = idx; // Position for v is not affected by height being h + 1;

		// First, update u.
		if (img_x == 0)
		{
			u[u_pos] = 0.0f;
		}
		else {
			if (img_x == w - 1)
			{
				u[u_pos + 1] = 0.0f;
			}
			if (!(M[idx] && M[idx - 1]))
			{
				u[u_pos] = 0.0f;
			}
		}

		// Next, update v.
		if (img_y == 0)
		{
			v[v_pos] = 0.0f;
		}
		else {
			if (img_y == h - 1)
			{
				v[v_pos + w] = 0.0f;
			}
			if (!(M[idx] && M[idx - w]))
			{
				v[v_pos] = 0.0f;
			}
		}
	}
}

__global__ void process_calc_Mprime(bool* M, float* Mprime, float* Mprime_Kernel, int k_radius, float MprimeKernelTotal, int w, int h, int x0, int y0, int x1, int y1)
{
	// Calculates Mprime values using the gaussian Mprime_Kernel, in the window defined by x0, y0 to x1, y1 (inclusive). Assumes x1 > x0, and y1 > y0.
	// Mprime_Kernel is one quadrant of the 2*k_radius+1 by 2*k_radius+1 matrix.

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	int N = (x1 - x0 + 1) * (y1 - y0 + 1);
	if (idx < N)
	{
		int window_width = (x1 - x0 + 1);
		int img_x = idx % window_width + x0;
		int img_y = idx / window_width + y0;
		int img_index = img_x + img_y * w;
		float accumulation = 0.0f;
		for (int wx = img_x - k_radius; wx <= img_x + k_radius; ++wx) // Window for kernel centered on img_x, img_y.
		{
			for (int wy = img_y - k_radius; wy <= img_y + k_radius; ++wy)
			{
				if ((wx >= 0) && (wx < w) && (wy >= 0) && (wy < h))
				{
					if (M[wx + wy * w])
					{
						accumulation += Mprime_Kernel[abs(wx - img_x) + abs(wy - img_y) * (k_radius + 1)];
					}
				}
			}
		}
		Mprime[img_index] = 1.0f - (accumulation / MprimeKernelTotal);
	}
}

__global__ void process_flow_outward(bool* M, float* Mprime, float* p, int w, int h)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int N = w * h;
	if (idx < N)
	{
		if (M[idx])
		{
			//p[i][j] = std::max((float)(p[i][j] - eta * M_prime[i][j]), -5.0f); // Old version.
			p[idx] = (float)(p[idx] - eta * Mprime[idx]);
		}
	}
}

__global__ void process_active_chunks(bool* M, bool* M_chunks, int x_chunks, int blocks_per_chunk, int w, int h)
{
	// This kernel uses x and y for blocks.  Each block is block_dimension by block_dimension in size, and each
	// chunk is chunk_size by chunk_size in size, where chunk_size is a multiple of block_dimension.
	// M_chunks needs to be initialized to all false before this kernel is called.

	int chunk_x = blockIdx.x / blocks_per_chunk;
	int chunk_y = blockIdx.y / blocks_per_chunk;
	int chunk_num = chunk_x + chunk_y * x_chunks;

	if (!M_chunks[chunk_num]) // If the chunk has already been identified as true, no further processing needed.
	{
		int index_x = blockIdx.x * blockDim.x + threadIdx.x;
		int index_y = blockIdx.y * blockDim.y + threadIdx.y;
		if ((index_x < w) && (index_y < h))
		{
			if (M[index_x + index_y * w])
			{
				M_chunks[chunk_num] = true;
			}
		}
	}
}

__global__ void process_dab_step1(bool* M, bool* chunk_updates, float* s, float* p, int w, int h, int x, int y, int wx, int wy, int w_width, int w_height, int radius_squared, float saturation, float pressure)
{
	// Index is based on the window that starts at wx, wy and is w_width wide and w_height high.
	// x,y is the centerpoint of the circle.
	// w,h is the full image width and height (relevant for calculating index into M, s, and p).
	// chunk_updates is a matrix of which chunks in the g SparseMatrix are touched and may need to be allocated if not already allocated.

	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int N = w_width * w_height;
	if (idx < N)
	{
		int img_x = wx + idx % w_width;
		int img_y = wy + idx / w_width;
		int img_idx = img_x + img_y * w;
		int dx = (img_x - x);
		int dy = (img_y - y);
		if ((img_x >= 0) && (img_y >= 0) && (img_x < w) && (img_y < h) && (dx * dx + dy * dy <= radius_squared))
		{
			s[img_idx] += saturation;
			M[img_idx] = true;
			p[img_idx] += pressure;
		}
		int x_chunks = (w + chunk_size - 1) / chunk_size;
		int chunk_num = x_chunks * (img_y / chunk_size) + img_x / chunk_size;
		if (!chunk_updates[chunk_num])
		{
			chunk_updates[chunk_num] = true;
		}
	}
}

__global__ void process_dab_step2(float* chunk, float concentration, int w, int h, int x, int y, int wx, int wy, int w_width, int w_height, int chunk_x0, int chunk_y0, int radius_squared)
{
	// Index, x,y and w,h are based on the same window as process_dab_step1.
	// chunk_x0 and chunk_y0 are the image-based location of the first pixel of the chunk to be processed.
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int N = w_width * w_height;
	if (idx < N)
	{
		int img_x = wx + idx % w_width;
		int img_y = wy + idx / w_width;
		int chunk_x = img_x - chunk_x0;
		int chunk_y = img_y - chunk_y0;
		if ((chunk_x >= 0) && (chunk_y >= 0) && (chunk_x < chunk_size) && (chunk_y < chunk_size)) // Only proceed if point is in this chunk.
		{
			int dx = (img_x - x);
			int dy = (img_y - y);
			if ((img_x >= 0) && (img_y >= 0) && (img_x < w) && (img_y < h) && (dx * dx + dy * dy <= radius_squared))  // Only proceed if point is in the image and in the circle.
			{
				int chunk_idx = chunk_x + chunk_y * chunk_size;
				chunk[chunk_idx] += concentration;
			}
		}
	}
}

__global__ void process_paint_step1(bool* M, bool* chunk_updates, float* s, float* p, int w, int h, int* SPdata, int identifier, int wx, int wy, int w_width, int w_height, float saturation, float pressure)
{
	// Index is based on the window that starts at wx, wy and is w_width wide and w_height high.
	// SPdata is the SuperPixel data with one identifier per pixel.
	// identifier is the number that indicates which pixels are to be painted.
	// w,h is the full image width and height (relevant for calculating index into M, s, and p).
	// chunk_updates is a matrix of which chunks in the g SparseMatrix are touched and may need to be allocated if not already allocated.

	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int N = w_width * w_height;
	if (idx < N)
	{
		int img_x = wx + idx % w_width;
		int img_y = wy + idx / w_width;
		int img_idx = img_x + img_y * w;

		if ((img_x >= 0) && (img_y >= 0) && (img_x < w) && (img_y < h) && (identifier == SPdata[img_idx]))
		{
			s[img_idx] += saturation;
			M[img_idx] = true;
			p[img_idx] += pressure;
			int x_chunks = (w + chunk_size - 1) / chunk_size;
			int chunk_num = x_chunks * (img_y / chunk_size) + img_x / chunk_size;
			if (!chunk_updates[chunk_num])
			{
				chunk_updates[chunk_num] = true;
			}
		}
	}
}

__global__ void process_paint_step2(float* chunk, float concentration, int w, int h, int* SPdata, int identifier, int wx, int wy, int w_width, int w_height, int chunk_x0, int chunk_y0)
{
	// Index, x,y and w,h are based on the same window as process_paint_step1.
	// chunk_x0 and chunk_y0 are the image-based location of the first pixel of the chunk to be processed.
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int N = w_width * w_height;
	if (idx < N)
	{
		int img_x = wx + idx % w_width;
		int img_y = wy + idx / w_width;
		int img_idx = img_x + img_y * w;
		int chunk_x = img_x - chunk_x0;
		int chunk_y = img_y - chunk_y0;
		if ((chunk_x >= 0) && (chunk_y >= 0) && (chunk_x < chunk_size) && (chunk_y < chunk_size)) // Only proceed if point is in this chunk.
		{
			if ((img_x >= 0) && (img_y >= 0) && (img_x < w) && (img_y < h) && (identifier == SPdata[img_idx]))  // Only proceed if point is in the image and has the right identifier.
			{
				int chunk_idx = chunk_x + chunk_y * chunk_size;
				chunk[chunk_idx] += concentration;
			}
		}
	}
}

__global__ void process_bristle_dab(cudaBrush* brush, int num_bristles, FloatPointPair direction, int* mask, int mask_width, int mask_height, int mask_value, float scaled_spot_radius_squared, bool begin, float paint_scale, float saturation, float pressure, float concentration, bool* chunk_updates)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x; // Bristle number.
	if (idx < num_bristles)
	{
		bool new_stroke = false; // Default
		int bristle_kernel_side = 1 + BRISTLE_KERNEL / 2;
		cudaBristle* local_bristle = &brush->bristles[idx];
		float flow_diff = local_bristle->flow_difference;
		float flow = brush->paint_prop.flow;
		float c_value = brush->orient_cosine;
		float s_value = brush->orient_sine;
		int width = brush->width;
		int height = brush->height;

		// Calculate the oriented_location of the bristle.
		FloatPointPair oriented_location;  // Offset location of bristle taking into account wander (relative to brush center).
		oriented_location.x = paint_scale * (c_value * (local_bristle->offset.x + local_bristle->wander.x) - s_value * (local_bristle->offset.y + local_bristle->wander.y));
		oriented_location.y = paint_scale * (c_value * (local_bristle->offset.y + local_bristle->wander.y) + s_value * (local_bristle->offset.x + local_bristle->wander.x));

		// The below values do not take into account paint_scale.
		FloatPointPair raw_offset = { local_bristle->offset.x , local_bristle->offset.y };  // Offset location of bristle without taking into account wander (relative to brush center).
		float x = c_value * raw_offset.x - s_value * raw_offset.y + brush->location.x;  // Absolute location of bristle (not relative to brush center), without wander and without paint_scale.
		float y = c_value * raw_offset.y + s_value * raw_offset.x + brush->location.y;

		// Calculate the actual bristle_location.
		FloatPointPair bristle_location = { paint_scale * brush->location.x + oriented_location.x, paint_scale * brush->location.y + oriented_location.y };  // Absolute location of bristle (not relative to brush center), with wander.

		float movement_distance = 0.0f;  // How far has bristle moved since last painted spot?
		float dot_product = 0.0f; // What is the direction of bristle movement when compared to brush movement?  Positive is same general direction.
		if (begin)
		{
			new_stroke = true;
			local_bristle->last_loc = bristle_location;
			local_bristle->down = true;
		}
		else
		{
			FloatPointPair bristle_movement = { bristle_location.x - local_bristle->last_loc.x, bristle_location.y - local_bristle->last_loc.y }; // Where has the bristle moved since the last Dab call?
			movement_distance = bristle_movement.x * bristle_movement.x + bristle_movement.y * bristle_movement.y;
			dot_product = bristle_movement.x * direction.x + bristle_movement.y * direction.y;
		}

		// Apply masking, if relevant.
		if (((NULL == mask) || ((x >= 0) && (x < mask_width) && (y >= 0) && (y < mask_height) && (mask_value == mask[(int)x + (int)y * mask_width]))) &&
			((scaled_spot_radius_squared < 0) || ((oriented_location.x * oriented_location.x + oriented_location.y * oriented_location.y) < scaled_spot_radius_squared)) &&
			(dot_product >= 0.0f))
		{
			if (false == local_bristle->down)  // Newly in writeable area.
			{
				new_stroke = true;
				local_bristle->last_loc = bristle_location;
				local_bristle->down = true;
			}
		}
		else {  // Outside of mask or writeable area.
			local_bristle->down = false;
		}

		// Paint bristle, if appropriate.
		if (local_bristle->down &&  // First, make sure bristle is down.
			(new_stroke || ((movement_distance >= 1.0))))  // But we also need to either have a new stroke or sufficient movement distance.
		{
			local_bristle->last_loc = bristle_location;
			float* local_bristle_kernel = brush->bristle_kernel;
			int x_chunks = (width + chunk_size - 1) / chunk_size;
			for (int i = -(bristle_kernel_side - 1); i < bristle_kernel_side; ++i)
			{
				for (int j = -(bristle_kernel_side - 1); j < bristle_kernel_side; ++j)
				{
					float adjustment;
					int img_x = (int)bristle_location.x + i;
					int img_y = (int)bristle_location.y + j;
					adjustment = local_bristle_kernel[abs(i) + bristle_kernel_side * abs(j)];

					float adjusted_flow = (flow_diff + flow) * adjustment;
					if ((img_x >= 0) && (img_x < width) && (img_y >= 0) && (img_y < height))
					{
						int img_idx = img_x + img_y * width;
						atomicAdd(&brush->s[img_idx], saturation);  // Need atomicAdd because other threads may be writing to the same pixel.
						//brush->s[img_idx] += saturation;
						brush->M[img_idx] = true;
						atomicAdd(&brush->p[img_idx], pressure);
						//brush->p[img_idx] += pressure;
						atomicAdd(&brush->full_g[img_idx], concentration);
						//brush->full_g[img_idx] += concentration;
						int chunk_num = x_chunks * (img_y / chunk_size) + img_x / chunk_size;
						if (!chunk_updates[chunk_num])
						{
							chunk_updates[chunk_num] = true;
						}
					}
				}
			}

		}
		//if (Adjust(gen) < 0.1)
		//{
		//	offset.x = Adjust_x(gen);
		//	offset.y = Adjust_y(gen);
		//	bristle->AdjustWander(offset);
		//}

	}
}

__global__ void calc_render_step1(float* reflected, float* upward_transmission, float* downward_transmission, float* g_chunk, float* d_chunk, float K, float S)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int N = chunk_size * chunk_size;
	if (idx < N)
	{
		if (upward_transmission[idx] > EFFECTIVE_ZERO)  // When no more light is being reflected back, no more calculation is needed.
		{
			float thickness = 0.0f;
			if (NULL != d_chunk)
			{
				thickness = d_chunk[idx];
			}
			if (NULL != g_chunk)
			{
				thickness += g_chunk[idx] * render_g_value;
			}
			if (thickness > EFFECTIVE_ZERO)
			{
				// Calculate R.
				float R;
				float a = 1.0 + K / S;
				float b = sqrt(a * a - 1);
				float d = b * S * thickness;
				float c = a * sinh(d) + b * cosh(d);
				if (c > EFFECTIVE_ZERO)
				{
					R = sinh(d) / c;
				}
				else {
					R = 0.0;
				}
				if (R > 1.0)
				{
					R = 1.0;
				}
				else if (R < 0.0) {
					R = 0.0;
				}
				// Calculate T.
				float T;
				a = 1.0 + K / S;
				b = sqrt(a * a - 1);
				d = b * S * thickness;
				c = a * sinh(d) + b * cosh(d);
				if (c > EFFECTIVE_ZERO)
				{
					T = b / c;
				}
				else {
					T = 0.0;
				}
				if (T > 1.0)
				{
					T = 1.0;
				}
				else if (T < 0.0) {
					T = 0.0;
				}
				reflected[idx] += downward_transmission[idx] * R * upward_transmission[idx];
				downward_transmission[idx] = downward_transmission[idx] * (1.0 - R) * T;
				upward_transmission[idx] = upward_transmission[idx] * T;
			}
		}
	}
}

__global__ void calc_render_step2(float* reflected, float* upward_transmission, float* downward_transmission, float substrate_color, unsigned char* d_image, int w, int h, int x0, int y0, int channel)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int N = chunk_size * chunk_size;
	if (idx < N)
	{
		float value = reflected[idx] + downward_transmission[idx] * upward_transmission[idx] * substrate_color;
		if (value > 1.0)
		{
			value = 1.0;
		}
		else if (value < 0.0)
		{
			value = 0.0;
		}
		int img_x = x0 + idx % chunk_size;
		int img_y = y0 + idx / chunk_size;
		if ((img_x >= 0) && (img_x < w) && (img_y >= 0) && (img_y < h))
		{
			int adjusted_idx = 3 * (img_x + img_y * w) + channel;
			d_image[adjusted_idx] = 255 * value;
		}
	}
}

__global__ void calc_move_pigment(int w, int h, bool* M, int chunk_num, float* u, float* v, float* g, float* g_prime, bool* d_chunk_updates, float step_size)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int N = chunk_size * chunk_size;
	if (idx < N)
	{
		int x_chunks = (w + chunk_size - 1) / chunk_size;
		int wx = (chunk_num % x_chunks) * chunk_size;
		int wy = (chunk_num / x_chunks) * chunk_size;
		int x = wx + idx % chunk_size;
		int y = wy + idx / chunk_size;
		bool left_flag = (0 == idx % chunk_size);
		bool right_flag = ((chunk_size - 1) == idx % chunk_size);
		bool top_flag = (0 == idx / chunk_size);
		bool bottom_flag = ((chunk_size - 1) == idx / chunk_size);

		if ((x >= 0) && (x < w) && (y >= 0) && (y < h))
		{
			int adjusted_idx = x + y * w;
			if (M[adjusted_idx])
			{
				int u_pos = x + y * (w + 1); // Because u is 1 pixel wider than w.
				int v_pos = adjusted_idx;
				if (!d_chunk_updates[chunk_num])  // *** Is this still needed, if only intersection chunks are called with this kernel? ***
				{
					d_chunk_updates[chunk_num] = true;
				}
				float g_move_left = 0.0f;
				float g_move_right = 0.0f;
				float g_move_top = 0.0f;
				float g_move_bottom = 0.0f;

				if (x > 0) // Left pixel available
				{
					g_move_left = -u[u_pos];
					if (g_move_left < 0.0)
					{
						g_move_left = 0.0;
					}
				}
				if (x < (w - 1)) // Right pixel available
				{
					g_move_right = u[u_pos + 1];
					if (g_move_right < 0.0)
					{
						g_move_right = 0.0;
					}
				}
				if (y > 0) // Top pixel available
				{
					g_move_top = -v[v_pos];
					if (g_move_top < 0.0)
					{
						g_move_top = 0.0;
					}
				}
				if (y < (h - 1)) // Bottom pixel available
				{
					g_move_bottom = v[v_pos + w];
					if (g_move_bottom < 0.0)
					{
						g_move_bottom = 0.0;
					}
				}
				float outflow_factor = g_move_left + g_move_right + g_move_top + g_move_bottom;  // Detect if total outflow would be higher than 1.0.
				if (outflow_factor <= 1.0f)
				{
					outflow_factor = pigment_lag * g[adjusted_idx] * step_size;
				}
				else {
					outflow_factor = pigment_lag * g[adjusted_idx] * step_size / outflow_factor;
				}
				if (outflow_factor > EFFECTIVE_ZERO)
				{
					if (x > 0)
					{
						atomicAdd(&g_prime[adjusted_idx - 1], g_move_left * outflow_factor);
						if (left_flag && (!d_chunk_updates[chunk_num - 1])) // If on the left edge of the chunk, we need to flag the previous chunk for updating.
						{
							d_chunk_updates[chunk_num - 1] = true;
						}
					}
					if (x < (w - 1))
					{
						atomicAdd(&g_prime[adjusted_idx + 1], g_move_right * outflow_factor);
						if (right_flag && (!d_chunk_updates[chunk_num + 1])) // If on the right edge of the chunk, we need to flag the next chunk for updating.
						{
							d_chunk_updates[chunk_num + 1] = true;
						}
					}
					if (y > 0)
					{
						atomicAdd(&g_prime[adjusted_idx - w], g_move_top * outflow_factor);
						if (top_flag && (!d_chunk_updates[chunk_num - x_chunks])) // If on the top edge of the chunk, we need to flag the chunk in the previous row for updating.
						{
							d_chunk_updates[chunk_num - x_chunks] = true;
						}
					}
					if (y < (h - 1))
					{
						atomicAdd(&g_prime[adjusted_idx + w], g_move_bottom * outflow_factor);
						if (bottom_flag && (!d_chunk_updates[chunk_num + x_chunks])) // If on the bottom edge of the chunk, we need to flag the chunk in the next row for updating.
						{
							d_chunk_updates[chunk_num + x_chunks] = true;
						}
					}
					atomicAdd(&g_prime[adjusted_idx], -(g_move_left + g_move_right + g_move_top + g_move_bottom) * outflow_factor);
				}
			}
		}
	}
}

__global__ void process_transfer_pigment(int w, int h, int chunk_index, float* chunk_g, float* chunk_d, float* thickness, int padding, float gamma, float rho, float omega, int x_chunks, bool* M, bool* M_chunks)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int N = chunk_size * chunk_size;
	//if (M_chunks[chunk_index])  // If this chunk does not contain any wet cells, then no need for further calculation.
	//{
	if (idx < N)  // Confirm that this is within the chunk.  This should always be true, if the chunk_size to threads_per_block ratio is right.
	{
		int x0 = chunk_size * (chunk_index % x_chunks); // X offset for chunk.
		int y0 = chunk_size * (chunk_index / x_chunks); // Y offset for chunk.
		int img_x = x0 + idx % chunk_size; // X location in image.
		int img_y = y0 + idx / chunk_size; // Y location in image.

		if ((img_x < w) && (img_y < h))  // Don't proceed if location is outside the bounds of the image.
		{
			int img_idx = img_x + img_y * w; // Index for location in image.
			int thickness_idx = padding + img_x + (2 * padding + w) * (padding + img_y);
			if (M[img_idx]) // Only proceed if this cell is wet.
			{
				float delta_down = chunk_g[idx] * (1.0f - thickness[thickness_idx] * gamma) * rho;
				float delta_up = chunk_d[idx] * (1.0f + (thickness[thickness_idx] - 1.0f) * gamma) * rho / omega;
				chunk_d[idx] += (delta_down - delta_up);
				chunk_g[idx] += (delta_up - delta_down);
			}
		}
	}
	//}
}

__global__ void process_capillary_flow_step1(int w, int h, float* s, float* thickness, int padding, bool* M, bool* M_chunks)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int N = w * h;
	//if (M_chunks[chunk_index])  // If this chunk does not contain any wet cells, then no need for further calculation.
//{
	if (idx < N)  // Confirm that this is within the chunk.  This should always be true, if the chunk_size to threads_per_block ratio is right.
	{
		if (M[idx])
		{
			int x = idx % w;
			int y = idx / w;
			int thickness_idx = padding + x + (2 * padding + w) * (padding + y);
			float capacity = thickness[thickness_idx] * (capacity_max - capacity_min) + capacity_min;
			float value = capacity - s[idx];
			if (absorption_alpha < value)
			{
				value = absorption_alpha;
			}
			else if (value < 0.0) {
				value = 0.0;
			}
			s[idx] += value;
		}
	}
	//}
}

__global__ void process_capillary_flow_step2(int w, int h, float* s, float* s_prime, float* thickness, int padding, bool* M, bool* M_chunks)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int N = w * h;
	//if (M_chunks[chunk_index])  // If this chunk does not contain any wet cells, then no need for further calculation.
//{
	if (idx < N)  // Confirm that this is within the chunk.  This should always be true, if the chunk_size to threads_per_block ratio is right.
	{
		if (M[idx])
		{
			int x = idx % w;
			int y = idx / w;
			int thickness_idx = padding + x + (2 * padding + w) * (padding + y);
			float capacity = thickness[thickness_idx] * (capacity_max - capacity_min) + capacity_min;
			// k and l will give the four direct neighbors to i,j.
			int k = 0;
			int l = 0;
			for (int neighbor = 0; neighbor < 4; ++neighbor)  // Step through four adjacent neighbors.  Left, Up, Right, Down.
			{
				unsigned char bit0 = neighbor & 1;
				unsigned char bit1 = (neighbor & 2) >> 1;
				int dx = (1 - bit0) * (2 * bit1 - 1);
				int dy = bit0 * (2 * bit1 - 1);
				k = x + dx;
				l = y + dy;
				int neighbor_idx = idx + dx + dy * w;
				if ((k >= 0) && (k < w) && (l >= 0) && (l < h)) // Check to see if k or l is out of bounds.
				{
					if ((s[idx] > saturation_epsilon * capacity) && (s[idx] > s[neighbor_idx]))
					{
						float neighbor_capacity = thickness[thickness_idx + dx + dy * (2 * padding + w)] * (capacity_max - capacity_min) + capacity_min;
						float del_s = s[idx] - s[neighbor_idx];
						float value = neighbor_capacity - s[neighbor_idx];
						if (del_s > value)
						{
							del_s = value;
						}
						del_s = del_s / 4.0;
						if (saturation_max_diffusion < del_s)
						{
							del_s = saturation_max_diffusion;
						}
						if (del_s > 0.0)
						{
							atomicAdd(&s_prime[idx], -del_s);
							atomicAdd(&s_prime[neighbor_idx], del_s);
						}
					}
				}
			}
		}
	}
}

__global__ void process_capillary_flow_step3(int w, int h, float* s, float* s_prime, float* p, float* thickness, int padding, bool* M, bool* M_chunks)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int N = w * h;
	if (idx < N)  // Confirm that this is within the chunk.  This should always be true, if the chunk_size to threads_per_block ratio is right.
	{
		int x = idx % w;
		int y = idx / w;
		int thickness_idx = padding + x + (2 * padding + w) * (padding + y);
		float capacity = thickness[thickness_idx] * (capacity_max - capacity_min) + capacity_min;
		if ((!M[idx]) && ((s[idx] > saturation_sigma) || (s[idx] >= (capacity - 2.0 * absorption_alpha))))
		{
			M[idx] = true;
			p[idx] = 0.0f;
		}
	}
}

__global__ void process_dry(float* g_chunk, float* d_chunk)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int N = chunk_size * chunk_size;
	if (idx < N)  // Confirm that this is within the chunk.  This should always be true, if the chunk_size to threads_per_block ratio is right.
	{
		d_chunk[idx] += render_g_value * g_chunk[idx];
	}
}

__global__ void process_brush_location(cudaBrush* brush, FloatPointPair loc)
{
	brush->location.x = loc.x;
	brush->location.y = loc.y;
}

__global__ void process_brush_orientation(cudaBrush* brush, float orientation)
{
	brush->orientation = orientation;
	brush->orient_cosine = cos(orientation);
	brush->orient_sine = sin(orientation);
}

bool Convolve(float* input, int matrix_width, int matrix_height, int padding, float* output, float* in_kernel, float scale, int k_radius, float background)
{
	bool ret = true;
	cudaError_t cudaStatus;
	int N = matrix_width * matrix_height;
	int x = matrix_width + 2 * padding;
	int y = matrix_height + 2 * padding;
	int blocksPerGrid = (x * y + threads_per_block - 1) / threads_per_block;
	int k_size = k_radius + 1;

	if (ret)
	{
		process_convolution << <blocksPerGrid, threads_per_block >> > (input, output, in_kernel, scale, matrix_width, matrix_height, padding, k_radius);
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess)
		{
			throw std::runtime_error("Failed to process convolution in Convolve.\n");
			ret = false;
		}
	}
	return ret;
}

bool SlopeVelocities(float* thickness, float* u, float* v, int w, int h, int padding)
{
	bool ret = true;
	cudaError_t cudaStatus;
	int blocksPerGrid = (w * h + threads_per_block - 1) / threads_per_block;
	process_slope_velocities << <blocksPerGrid, threads_per_block >> > (thickness, u, v, w, h, padding);
	cudaDeviceSynchronize();
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed to update velocities in SlopeVelocities.\n");
		ret = false;
	}
	return ret;
}

float GetMaxVelocity(float* u, float* v, int w, int h, float* results1, float* results2, int results_length)
{
	float ret = 0.0;
	cudaError_t cudaStatus;
	int blocksPerGrid = (w * h + threads_per_block - 1) / threads_per_block;
	process_max_velocity_2 << <blocksPerGrid, threads_per_block >> > (w, h, u, v, results1);
	cudaDeviceSynchronize();
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed to calculate max velocities in GetMaxVelocity.\n");
		ret = -1.0;
	}
	results_length = blocksPerGrid;
	bool result_in_1 = true;
	while (results_length > 1)
	{
		result_in_1 = !result_in_1;
		blocksPerGrid = (results_length + threads_per_block - 1) / threads_per_block;
		if (result_in_1)
		{
			get_max_results_2 << <blocksPerGrid, threads_per_block >> > (results2, results_length, results1);
		}
		else {
			get_max_results_2 << <blocksPerGrid, threads_per_block >> > (results1, results_length, results2);
		}
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess)
		{
			throw std::runtime_error("Failed to calculate max velocities in GetMaxVelocity.\n");
			ret = -1.0;
		}
		results_length = blocksPerGrid;
	}
	if (result_in_1)
	{
		cudaStatus = cudaMemcpy(&ret, results1, sizeof(float), cudaMemcpyDeviceToHost);
	}
	else {
		cudaStatus = cudaMemcpy(&ret, results2, sizeof(float), cudaMemcpyDeviceToHost);
	}
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed to convert max velocity to host variable in GetMaxVelocity.\n");
		ret = -1.0;
	}
	return ret;
}

bool CalcVelocities(float dt, int w, int h, float* u, float* v, float* u_prime, float* v_prime, float* p, bool* M)
{
	bool ret = true;
	cudaError_t cudaStatus;
	int blocksPerGrid = (w * h + threads_per_block - 1) / threads_per_block;
	bool keep_going = true;
	float t = 0.0f; // To keep track of the progress through the full time step.
	float step_size = dt; // The step size, which is also the scale factor for u and v, so that they take into account the time step size.

	while (ret && keep_going)  // t is updated at the end of the step.
	{
		process_calc_velocities << <blocksPerGrid, threads_per_block >> > (step_size, w, h, u, v, u_prime, v_prime, p);
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess)
		{
			throw std::runtime_error("Failed to calculate velocities in CalcVelocities.\n");
			ret = false;
		}

		ret = ret && (CopyFloatArray(u_prime, u, w + 1, h) && CopyFloatArray(v_prime, v, w, h + 1));

		if (ret)
		{
			process_enforce_boundaries << <blocksPerGrid, threads_per_block >> > (w, h, u, v, M);
			cudaDeviceSynchronize();
			cudaStatus = cudaGetLastError();
			if (cudaStatus != cudaSuccess)
			{
				throw std::runtime_error("Failed to enforce boundaries in CalcVelocities.\n");
				ret = false;
			}
		}

		t += step_size;
		if (t >= 1.0f)
		{
			keep_going = false;
		}
		else {
			step_size = std::min(step_size, 1.0f - t);
			if (step_size < EFFECTIVE_ZERO)
			{
				keep_going = false;
			}
		}
	}
	return ret;
}

bool CalcRelaxDivergence(bool* M, float* u, float* v, float* p, float* delta_matrix, int w, int h)
{
	bool ret = true;
	cudaError_t cudaStatus;
	int blocksPerGrid = (w * h + threads_per_block - 1) / threads_per_block;
	int count = relaxation_steps;
	while (ret && (count > 0))
	{
		process_relax_divergence_step1 << <blocksPerGrid, threads_per_block >> > (M, u, v, delta_matrix, w, h);
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess)
		{
			throw std::runtime_error("Failed in first step of CalcRelaxDivergence.\n");
			ret = false;
		}
		if (ret)
		{
			process_relax_divergence_step2 << <blocksPerGrid, threads_per_block >> > (M, u, v, p, delta_matrix, w, h);
			cudaDeviceSynchronize();
			cudaStatus = cudaGetLastError();
			if (cudaStatus != cudaSuccess)
			{
				throw std::runtime_error("Failed in second step of CalcRelaxDivergence.\n");
				ret = false;
			}
		}
		count--;
	}
	return ret;
}

bool CalcRelaxDivergence_b(bool* M, float* u, float* v, float* p, float* delta_matrix, int w, int h)
{
	bool ret = true;
	cudaError_t cudaStatus;
	int blocksPerGrid = (w * h + threads_per_block - 1) / threads_per_block;
	int count = relaxation_steps;
	// Use 2D blocks and grid for call to second step (to simplify shared memory).
	dim3 threadsPerBlock(block_dimension, block_dimension);
	dim3 numBlocks((w + block_dimension - 1) / block_dimension, (h + block_dimension - 1) / block_dimension);
	while (ret && (count > 0))
	{
		process_relax_divergence_step1 << <blocksPerGrid, threads_per_block >> > (M, u, v, delta_matrix, w, h);
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess)
		{
			throw std::runtime_error("Failed in first step of CalcRelaxDivergence.\n");
			ret = false;
		}
		if (ret)
		{
			process_relax_divergence_step2_b << <numBlocks, threadsPerBlock >> > (M, u, v, p, delta_matrix, w, h);
			cudaDeviceSynchronize();
			cudaStatus = cudaGetLastError();
			if (cudaStatus != cudaSuccess)
			{
				throw std::runtime_error("Failed in second step of CalcRelaxDivergence.\n");
				ret = false;
			}
		}
		count--;
	}
	return ret;
}

bool CalcMprime(bool* M, float* Mprime, float* Mprime_Kernel, int k_radius, float MprimeKernelTotal, int w, int h, int x0, int y0, int x1, int y1)
{
	bool ret = true;
	cudaError_t cudaStatus;
	int blocksPerGrid = ((x1 - x0 + 1) * (y1 - y0 + 1) + threads_per_block - 1) / threads_per_block;
	process_calc_Mprime << <blocksPerGrid, threads_per_block >> > (M, Mprime, Mprime_Kernel, k_radius, MprimeKernelTotal, w, h, x0, y0, x1, y1);
	cudaDeviceSynchronize();
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed to process calculation of M_Prime.\n");
		ret = false;
	}
	return ret;
}

bool CalcFlowOutward(bool* M, float* Mprime, float* p, int w, int h)
{
	bool ret = true;
	cudaError_t cudaStatus;
	int blocksPerGrid = (w * h + threads_per_block - 1) / threads_per_block;
	process_flow_outward << <blocksPerGrid, threads_per_block >> > (M, Mprime, p, w, h);
	cudaDeviceSynchronize();
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed to process process_flow_outward in CalcFlowOutward.\n");
		ret = false;
	}
	return ret;
}

bool CalcActiveChunks(bool* M, bool* M_chunks, bool* host_M_chunks, int x_chunks, int y_chunks, int w, int h)
{
	bool ret = true;
	cudaError_t cudaStatus;
	int blocksPerGrid = (x_chunks * y_chunks + threads_per_block - 1) / threads_per_block;
	dim3 blockDim(block_dimension, block_dimension);
	dim3 gridDim(
		(w + blockDim.x - 1) / blockDim.x,
		(h + blockDim.y - 1) / blockDim.y
	);
	// First, we have to set M_chunks to all false.
	initialize_Bool_Array << <blocksPerGrid, threads_per_block >> > (M_chunks, false, x_chunks * y_chunks);
	cudaDeviceSynchronize();
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed to initialize M_chunks in CalcActiveChunks.\n");
		ret = false;
	}
	// Now, set M_chunks to true where any M is true.
	if (ret)
	{
		process_active_chunks << <gridDim, blockDim >> > (M, M_chunks, x_chunks, chunk_size / block_dimension, w, h);
		cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess)
		{
			throw std::runtime_error("Failed to calculate M_chunks in CalcActiveChunks.\n");
			ret = false;
		}
	}
	if (ret)
	{
		cudaStatus = cudaMemcpy(host_M_chunks, M_chunks, x_chunks * y_chunks * sizeof(bool), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess)
		{
			throw std::runtime_error("Failed to transfer device data to host_M_chunks.\n");
			ret = false;
		}
	}
	return ret;
}

bool CalcDab(bool* M, float* s, float* p, SparseFloatMatrix* g, int w, int h, int x, int y, int radius, float saturation, float concentration)
{
	bool ret = true;
	cudaError_t cudaStatus;
	int x0 = std::min(w - 2, std::max(1, x - radius));
	int y0 = std::min(h - 2, std::max(1, y - radius));
	int x1 = std::min(w - 2, std::max(1, x + radius));
	int y1 = std::min(h - 2, std::max(1, y + radius));
	int w_width = (x1 - x0 + 1);
	int w_height = (y1 - y0 + 1);
	int r2 = radius * radius;
	int blocksPerGrid = (w_width * w_height + threads_per_block - 1) / threads_per_block;

	bool* d_chunk_updates = g->GetDChunkUpdates();
	bool* h_chunk_updates = g->GetHChunkUpdates();
	ret = ResetBoolArray(d_chunk_updates, g->GetXChunks(), g->GetYChunks(), false);
	if (!ret)
	{
		throw std::runtime_error("Failed to reset d_chunk_updates in CalcDab.\n");
	}
	else {
		process_dab_step1 << <blocksPerGrid, threads_per_block >> > (M, d_chunk_updates, s, p, w, h, x, y, x0, y0, w_width, w_height, r2, saturation, dab_pressure);
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess)
		{
			throw std::runtime_error("Failed in process_dab_step1.\n");
			ret = false;
		}
		ret = ret && g->UpdateChunks();  // Updates chunk_set as well as h_chunk_updates.
		if (!ret)
		{
			throw std::runtime_error("Failed to update chunks in g for CalcDab.\n");

		}


		std::set<int> chunks = g->GetChunkSet();
		int x_chunks = g->GetXChunks();
		int y_chunks = g->GetYChunks();
		for (int chunk_num = 0; chunk_num < (x_chunks * y_chunks); ++chunk_num)
		{
			if (h_chunk_updates[chunk_num])
			{
				float* chunk = g->GetChunk(chunk_num);
				int chunk_x0 = (chunk_num % x_chunks) * chunk_size;
				int chunk_y0 = (chunk_num / x_chunks) * chunk_size;
				process_dab_step2 << <blocksPerGrid, threads_per_block >> > (chunk, concentration, w, h, x, y, x0, y0, w_width, w_height, chunk_x0, chunk_y0, r2);
			}
		}
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess)
		{
			throw std::runtime_error("Failed in process_dab_step2.\n");
			ret = false;
		}
	}
	return ret;
}

bool CSetVelocity(float* u, float* v, int x, int y, int w, int h, float u_vel, float v_vel, bool sum, bool wait)
{
	bool ret = true;
	cudaError_t cudaStatus;
	set_single_value << <1, 1 >> > (u, w + 1, h, x, y, u_vel, sum);
	set_single_value << <1, 1 >> > (u, w + 1, h, x + 1, y, u_vel, sum);
	set_single_value << <1, 1 >> > (v, w, h + 1, x, y, v_vel, sum);
	set_single_value << <1, 1 >> > (v, w, h + 1, x, y + 1, v_vel, sum);

	if (wait)
	{
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess)
		{
			throw std::runtime_error("Failed to set velocity in CSetVelocity.\n");
			ret = false;
		}
	}
	return true;
}

bool CRender(std::vector<Pigment*> pigments, unsigned char* data, int w, int h, float* substrate_color)
{
	// Performs the rendering to the provided data matrix (w*h*3 unsigned chars).
	// substrate_color is a 3-element vector of floats between 0.0 and 1.0.
	bool ret = true;
	cudaError_t cudaStatus;
	int blocksPerGrid = (chunk_size * chunk_size + threads_per_block - 1) / threads_per_block;
	unsigned char* d_image = NULL;
	cudaStatus = cudaMalloc(&d_image, sizeof(unsigned char) * w * h * 3);
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed to allocate memory for d_image in CRender.\n");
		ret = false;
	}
	float* reflected = FloatArray(chunk_size, chunk_size, true, 0.0f);
	float* upward_transmission = FloatArray(chunk_size, chunk_size, true, 1.0f);
	float* downward_transmission = FloatArray(chunk_size, chunk_size, true, 1.0f);
	if (ret)
	{
		int num_pgmnt = pigments.size();
		if (num_pgmnt > 0)
		{
			int x_chunks = pigments[0]->GetXChunks();
			int y_chunks = pigments[0]->GetYChunks();
			for (int chunk_num = 0; ret && (chunk_num < (x_chunks * y_chunks)); ++chunk_num)
			{
				int x0 = (chunk_num % x_chunks) * chunk_size;
				int y0 = (chunk_num / x_chunks) * chunk_size;
				for (int channel = 0; channel < 3; ++channel)
				{
					ResetFloatArray(reflected, chunk_size, chunk_size, 0.0f);
					ResetFloatArray(upward_transmission, chunk_size, chunk_size, 1.0f);
					ResetFloatArray(downward_transmission, chunk_size, chunk_size, 1.0f);
					for (int p_num = num_pgmnt - 1; ret && (p_num >= 0); --p_num)
					{
						Pigment* pgmnt = pigments[p_num];
						SparseFloatMatrix* g = pgmnt->Get_g();
						SparseFloatMatrix* d = pgmnt->Get_d();
						float* g_chunk = NULL;
						float* d_chunk = NULL;
						if (g->CheckChunkNumber(chunk_num))
						{
							g_chunk = g->GetChunk(chunk_num);
						}
						if (d->CheckChunkNumber(chunk_num))
						{
							d_chunk = d->GetChunk(chunk_num);
						}
						if ((NULL != g_chunk) && (NULL == d_chunk))
						{
							throw std::runtime_error("Lack of consistency between g and d in CRender.\n");
							ret = false;
						}
						if (ret && ((NULL != g_chunk) || (NULL != d_chunk)))
						{
							calc_render_step1 << <blocksPerGrid, threads_per_block >> > (reflected, upward_transmission, downward_transmission, g_chunk, d_chunk, pgmnt->GetK(channel), pgmnt->GetS(channel));
							cudaDeviceSynchronize();
							cudaStatus = cudaGetLastError();
							if (cudaStatus != cudaSuccess)
							{
								throw std::runtime_error("Failed in calc_render_step1.\n");
								ret = false;
							}
						}
					}
					if (ret)
					{
						calc_render_step2 << <blocksPerGrid, threads_per_block >> > (reflected, upward_transmission, downward_transmission, substrate_color[channel], d_image, w, h, x0, y0, channel);
						cudaDeviceSynchronize();
						cudaStatus = cudaGetLastError();
						if (cudaStatus != cudaSuccess)
						{
							throw std::runtime_error("Failed in calc_render_step2.\n");
							ret = false;
						}
					}
				}
			}
		}
	}
	FreeFloatArray(reflected);
	FreeFloatArray(upward_transmission);
	FreeFloatArray(downward_transmission);
	if (ret)
	{
		cudaStatus = cudaMemcpy(data, d_image, 3 * w * h * sizeof(unsigned char), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess)
		{
			throw std::runtime_error("Failed to copy memory from device in CRender.\n");
			ret = false;
		}
	}
	if (NULL != d_image)
	{
		cudaFree(d_image);
	}
	return true;
}

bool CMovePigment(int w, int h, float* u, float* v, bool* M, bool* host_M_chunks, std::vector<Pigment*> pigments, float dt, float* g, float* g_prime)
{
	// dt is the time-step, based on the maximum velocity.

	bool ret = true;
	cudaError_t cudaStatus;
	int blocksPerGrid = (chunk_size * chunk_size + threads_per_block - 1) / threads_per_block;
	float t = 0.0f; // To keep track of the progress through the full time step.

	if (pigments.size() > 0)
	{
		// Create a set of all chunks that have any true values in the M mask.
		std::set<int> active_M_chunks;
		active_M_chunks.clear();
		int x_chunks = pigments[0]->GetXChunks();
		int y_chunks = pigments[0]->GetYChunks();
		for (int i = 0; i < x_chunks * y_chunks; ++i)
		{
			if (host_M_chunks[i])
			{
				active_M_chunks.insert(i);
			}
		}

		std::vector<Pigment*>::iterator k;
		for (k = pigments.begin(); (ret && (k != pigments.end())); ++k)
		{
			Pigment* k_pntr = *k;
			SparseFloatMatrix* sparse_g = k_pntr->Get_g();
			// Calculate the intersection, if any, between the M mask and the chunks for which this pigment has g values.
			std::set<int> intersection;
			intersection.clear();
			//std::set_intersection(sparse_g->GetChunkSet().begin(), sparse_g->GetChunkSet().end(), active_M_chunks.begin(), active_M_chunks.end(), std::inserter(intersection, intersection.begin()));
			for (std::set<int>::iterator chunk_it = active_M_chunks.begin(); chunk_it != active_M_chunks.end(); ++chunk_it)
			{
				if (sparse_g->CheckChunkNumber(*chunk_it))
				{
					intersection.insert(*chunk_it);
				}
			}
			if (intersection.size() > 0)
			{
				ret = ret && sparse_g->ExpandToFloatArray(g);
				ret = ret && sparse_g->ExpandToFloatArray(g_prime);
				int x_chunks = sparse_g->GetXChunks();
				int y_chunks = sparse_g->GetYChunks();
				bool* d_chunk_updates = sparse_g->GetDChunkUpdates();
				ret = ret && ResetBoolArray(d_chunk_updates, x_chunks, y_chunks, false);
				if (!ret)
				{
					throw std::runtime_error("Failed to copy g or g_prime values in CMovePigment.\n");
				}
				float step_size = dt; // The step size, which is also the scale factor for u and v, so that they take into account the time step size.
				t = 0.0f;
				bool keep_going = true;
				while (ret && keep_going)  // t is updated at the end of the step.
				{
					//std::set<int> chunk_set = sparse_g->GetChunkSet();
					for (std::set<int>::iterator chunk_it = intersection.begin(); chunk_it != intersection.end(); ++chunk_it)  // Only process chunks in the intersection.
					{
						int chunk_index = *chunk_it;
						int wy = chunk_index / x_chunks;  // Y component of beginning of window represented by chunk.
						int wx = chunk_index - (wy * x_chunks); // X component of the beginning of windows represented by chunk. 
						wx = wx * chunk_size;
						wy = wy * chunk_size;
						calc_move_pigment << <blocksPerGrid, threads_per_block >> > (w, h, M, chunk_index, u, v, g, g_prime, d_chunk_updates, step_size);
					}
					cudaDeviceSynchronize();
					cudaStatus = cudaGetLastError();
					if (cudaStatus != cudaSuccess)
					{
						throw std::runtime_error("Failed to move pigements in calc_move_pigment.\n");
						ret = false;
					}
					ret = ret && CopyFloatArray(g_prime, g, w, h);
					if (!ret)
					{
						throw std::runtime_error("Failed to copy g values in CMovePigment.\n");
					}
					t += step_size;
					if (t >= 1.0f)
					{
						keep_going = false;
					}
					else {
						step_size = std::min(step_size, 1.0f - t);
						if (step_size < EFFECTIVE_ZERO)
						{
							keep_going = false;
						}
					}
				}
				ret = ret && sparse_g->CompressFromFloatArray(g);
				ret = ret && k_pntr->Get_d()->SyncChunks(sparse_g);
			}
		}
	}
	return true;
}

bool CTransferPigment(Pigment* pgmnt, float* thickness, int padding, bool* M, bool* M_chunks)
{
	bool ret = true;
	cudaError_t cudaStatus;
	int blocksPerGrid = (chunk_size * chunk_size + threads_per_block - 1) / threads_per_block;
	SparseFloatMatrix* pigment_g = pgmnt->Get_g();
	SparseFloatMatrix* pigment_d = pgmnt->Get_d();
	int w = pigment_g->GetWidth();
	int h = pigment_g->GetHeight();
	float gamma = pgmnt->Get_gamma();
	float rho = pgmnt->Get_rho();
	float omega = pgmnt->Get_omega();
	// *** Which is quicker?  To figure out which chunks are currently wet and only send those to the kernel, or send everything to the kernel and let it quickly use M_chunks to short-circuit those that are dry?
	// I am going with the first one.
	std::set<int> g_chunks = pigment_g->GetChunkSet();

	for (std::set<int>::iterator chunk_it = g_chunks.begin(); chunk_it != g_chunks.end(); ++chunk_it)
	{
		int chunk_index = *chunk_it;
		float* chunk_g = pigment_g->GetChunk(chunk_index);
		float* chunk_d = pigment_d->GetChunk(chunk_index);
		if ((NULL == chunk_g) || (NULL == chunk_d))
		{
			throw std::runtime_error("Chunk not allocated for g or d in CTransferPigment.\n");
			ret = false;
		}
		process_transfer_pigment << <blocksPerGrid, threads_per_block >> > (w, h, chunk_index, chunk_g, chunk_d, thickness, padding, gamma, rho, omega, pgmnt->GetXChunks(), M, M_chunks);
	}
	cudaDeviceSynchronize();  // Allow as many threads to run in parallel as possible.
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed to transfer pigements in process_transfer_pigment.\n");
		ret = false;
	}
	return ret;
}

bool CCapillaryFlow(int w, int h, float* s, float* s_prime, float* p, float* thickness, int gauss_radius, bool* M, bool* M_chunks)
{
	bool ret = true;
	cudaError_t cudaStatus;
	int N = w * h;
	int blocksPerGrid = (N + threads_per_block - 1) / threads_per_block;

	process_capillary_flow_step1 << <blocksPerGrid, threads_per_block >> > (w, h, s, thickness, gauss_radius, M, M_chunks);
	cudaDeviceSynchronize();
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed to update saturation in process_capillary_flow_step1.\n");
		ret = false;
	}
	ret = ret && CopyFloatArray(s, s_prime, w, h);
	if (!ret)
	{
		throw std::runtime_error("Failed to copy s to s_prime in CCapillaryFlow.\n");
	}
	if (ret)
	{
		process_capillary_flow_step2 << <blocksPerGrid, threads_per_block >> > (w, h, s, s_prime, thickness, gauss_radius, M, M_chunks);
		cudaDeviceSynchronize();

		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess)
		{
			throw std::runtime_error("Failed to move saturation values in process_capillary_flow_step2.\n");
			ret = false;
		}
	}
	ret = ret && CopyFloatArray(s_prime, s, w, h);
	if (!ret)
	{
		throw std::runtime_error("Failed to copy s_prime to s in CCapillaryFlow.\n");
	}
	if (ret)
	{
		process_capillary_flow_step3 << <blocksPerGrid, threads_per_block >> > (w, h, s, s_prime, p, thickness, gauss_radius, M, M_chunks);
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess)
		{
			throw std::runtime_error("Failed to update M and p values in process_capillary_flow_step3.\n");
			ret = false;
		}
	}
	return ret;
}

bool CDry(SparseFloatMatrix* pgmnt_g, SparseFloatMatrix* pgmnt_d)
{
	bool ret = true;
	cudaError_t cudaStatus;
	int N = chunk_size * chunk_size;
	int blocksPerGrid = (N + threads_per_block - 1) / threads_per_block;

	pgmnt_d->SyncChunks(pgmnt_g);
	if (0.0f != render_g_value)
	{
		int num_chunks = pgmnt_g->GetXChunks() * pgmnt_g->GetYChunks();
		std::set<int> chunk_set = pgmnt_g->GetChunkSet();
		for (std::set<int>::iterator chunk_it = chunk_set.begin(); chunk_it != chunk_set.end(); ++chunk_it)
		{
			int chunk_number = *chunk_it;
			process_dry << <blocksPerGrid, threads_per_block >> > (pgmnt_g->GetChunk(chunk_number), pgmnt_d->GetChunk(chunk_number));
		}
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess)
		{
			throw std::runtime_error("Failed to dry pigment in process_dry.\n");
			ret = false;
		}
	}
	pgmnt_g->Reset_Values();
	pgmnt_g->SyncChunks(pgmnt_d); // To allow for transfer of pigment from a dry layer if it is later painted over.
	return ret;
}

bool CPaintArea(bool* M, float* s, float* p, SparseFloatMatrix* g, int* SPdata, int w, int h, RectQuad window, int identifier, float saturation, float concentration)
{
	bool ret = true;
	cudaError_t cudaStatus;
	int x0 = window.x0;
	int y0 = window.y0;
	int x1 = window.x1;
	int y1 = window.y1;
	int w_width = (x1 - x0 + 1);
	int w_height = (y1 - y0 + 1);
	int blocksPerGrid = (w_width * w_height + threads_per_block - 1) / threads_per_block;

	bool* d_chunk_updates = g->GetDChunkUpdates();
	bool* h_chunk_updates = g->GetHChunkUpdates();
	ret = ResetBoolArray(d_chunk_updates, g->GetXChunks(), g->GetYChunks(), false);
	if (!ret)
	{
		throw std::runtime_error("Failed to reset d_chunk_updates in CPaintArea.\n");
	}
	else {
		process_paint_step1 << <blocksPerGrid, threads_per_block >> > (M, d_chunk_updates, s, p, w, h, SPdata, identifier, x0, y0, w_width, w_height, saturation, dab_pressure);
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess)
		{
			throw std::runtime_error("Failed in process_paint_step1.\n");
			ret = false;
		}
		ret = ret && g->UpdateChunks();  // Updates chunk_set as well as h_chunk_updates.
		if (!ret)
		{
			throw std::runtime_error("Failed to update chunks in g for CPaintArea.\n");
		}

		std::set<int> chunks = g->GetChunkSet();
		int x_chunks = g->GetXChunks();
		int y_chunks = g->GetYChunks();
		for (int chunk_num = 0; chunk_num < (x_chunks * y_chunks); ++chunk_num)
		{
			if (h_chunk_updates[chunk_num])
			{
				float* chunk = g->GetChunk(chunk_num);
				int chunk_x0 = (chunk_num % x_chunks) * chunk_size;
				int chunk_y0 = (chunk_num / x_chunks) * chunk_size;
				process_paint_step2 << <blocksPerGrid, threads_per_block >> > (chunk, concentration, w, h, SPdata, identifier, x0, y0, w_width, w_height, chunk_x0, chunk_y0);
			}
		}
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess)
		{
			throw std::runtime_error("Failed in process_paint_step2.\n");
			ret = false;
		}
	}
	return ret;
}

cudaBrush* CreateCudaBrush(cudaBrush* host_brush, std::vector<Bristle*> host_bristles)
{
	cudaBrush* device_brush_pointer = NULL;
	float* temp_bristle_kernel = host_brush->bristle_kernel;
	int bristle_kernel_side = 1 + BRISTLE_KERNEL / 2;

	cudaError_t cudaStatus;
	// First, allocate memory for bristles, and copy data.
	cudaStatus = cudaMalloc(&host_brush->bristles, sizeof(cudaBristle) * host_brush->num_bristles);
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed to allocate memory in CreateCudaBrush.\n");
		return NULL;
	}
	cudaBristle* temp_bristles = (cudaBristle*)malloc(sizeof(cudaBristle) * host_brush->num_bristles);
	if (NULL == temp_bristles)
	{
		throw std::runtime_error("Failed to allocate memory for temp_bristles in CreateCudaBrush.\n");
		return NULL;
	}
	int count = 0;
	for (std::vector<Bristle*>::iterator br_it = host_bristles.begin(); br_it != host_bristles.end(); ++br_it)
	{
		Bristle* bristle = *br_it;
		temp_bristles[count].down = bristle->GetBristleDown();
		temp_bristles[count].flow_difference = bristle->GetFlowDiff();
		temp_bristles[count].last_loc = bristle->GetLast();
		temp_bristles[count].offset = bristle->GetOffset();
		temp_bristles[count].wander = bristle->GetWander();
		++count;
	}
	cudaStatus = cudaMemcpy(host_brush->bristles, temp_bristles, sizeof(cudaBristle) * host_brush->num_bristles, cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed to copy bristle memory to device in CopyFromHost.\n");
		return NULL;
	}
	free(temp_bristles);

	// Next, allocate memory for bristle kernel, and copy data.
	cudaStatus = cudaMalloc(&host_brush->bristle_kernel, sizeof(float) * bristle_kernel_side * bristle_kernel_side);
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed to allocate bristle_kernel memory in CreateCudaBrush.\n");
		return NULL;
	}
	cudaStatus = cudaMemcpy(host_brush->bristle_kernel, temp_bristle_kernel, sizeof(float) * bristle_kernel_side * bristle_kernel_side, cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed to copy memory to device in CopyFromHost.\n");
		return NULL;
	}

	// And allocate memory for brush, and copy data (including pointers to device memory for bristles and bristle kernel).
	cudaStatus = cudaMalloc(&device_brush_pointer, sizeof(cudaBrush));
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed to allocate memory in CreateCudaBrush.\n");
		return NULL;
	}
	cudaStatus = cudaMemcpy(device_brush_pointer, host_brush, sizeof(cudaBrush), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed to copy memory to device in CopyFromHost.\n");
		return NULL;
	}
	return device_brush_pointer;
}

bool FreeCudaBrush(cudaBrush* brush, cudaBristle* bristles, float* bristle_kernel)
{
	bool ret = true;
	if (NULL != brush)
	{
		cudaFree(brush);
	}
	if (NULL != bristles)
	{
		cudaFree(bristles);
	}
	if (NULL != bristle_kernel)
	{
		cudaFree(bristle_kernel);
	}
	return ret;
}

bool SetCudaBrushLocation(cudaBrush* brush, FloatPointPair loc)
{
	bool ret = true;
	cudaError_t cudaStatus;
	process_brush_location << <1, 1 >> > (brush, loc);
	cudaDeviceSynchronize();
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed in process_brush_location.\n");
		ret = false;
	}
	return ret;
}

bool SetCudaBrushOrientation(cudaBrush* brush, float orientation)
{
	bool ret = true;
	cudaError_t cudaStatus;
	process_brush_orientation << <1, 1 >> > (brush, orientation);
	cudaDeviceSynchronize();
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		throw std::runtime_error("Failed in process_brush_orientation.\n");
		ret = false;
	}
	return ret;
}

bool CudaDab3(cudaBrush* brush, SparseFloatMatrix* g, float saturation, float concentration, FloatPointPair direction, int num_bristles, SPixelData* mask, int mask_value, float spot_radius, bool begin, float paint_scale)
{
	bool ret = true;
	cudaError_t cudaStatus;
	if (num_bristles > 0)
	{
		int blocksPerGrid = (num_bristles + threads_per_block - 1) / threads_per_block;  // One thread per bristle.

		// Use direction to measure whether the movement has been forward.
		if ((abs(direction.x) > EFFECTIVE_ZERO) || (abs(direction.y) > EFFECTIVE_ZERO))
		{
			float magnitude = sqrt(direction.x * direction.x + direction.y * direction.y);
			direction.x = direction.x / magnitude;
			direction.y = direction.y / magnitude;
		}

		float scaled_spot_radius_squared = spot_radius * spot_radius * paint_scale * paint_scale;
		bool new_stroke = false;
		int bristle_kernel_side = 1 + BRISTLE_KERNEL / 2;
		int mask_width = 0;
		int mask_height = 0;
		int* device_mask = NULL;
		if (NULL != mask)
		{
			mask_width = mask->GetWidth();
			mask_height = mask->GetHeight();
			device_mask = mask->GetDeviceData();
		}
		process_bristle_dab << < blocksPerGrid, threads_per_block >> > (brush, num_bristles, direction, device_mask, mask_width, mask_height, mask_value, scaled_spot_radius_squared, begin, paint_scale, saturation, dab_pressure, concentration, g->GetDChunkUpdates());
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess)
		{
			throw std::runtime_error("Failed in process_bristle_dab.\n");
			ret = false;
		}
	}
	return ret;
}

bool StartStroke(SparseFloatMatrix* sparse_g, float* full_g, int width, int height)
{
	bool ret = true;
	int x_chunks = sparse_g->GetXChunks();
	int y_chunks = sparse_g->GetYChunks();
	ret = ResetFloatArray(full_g, width, height, 0);
	ret = ret && ResetBoolArray(sparse_g->GetDChunkUpdates(), x_chunks, y_chunks, false);  // Reset the device version of the chunk updates.
	ret = ret && CopyToHost(sparse_g->GetDChunkUpdates(), x_chunks * y_chunks, sparse_g->GetHChunkUpdates()); // Copy cleared version to host.
	return ret;
}

bool EndStroke(SparseFloatMatrix* sparse_g, float* full_g, int width, int height)
{
	bool ret = true;
	cudaError_t cudaStatus;
	int blocksPerGrid = (chunk_size * chunk_size + threads_per_block - 1) / threads_per_block;
	ret = sparse_g->UpdateChunks(); // Copy d_chunk_updates to h_chunk_updates, and allocate any new chunks needed.
	if (!ret)
	{
		throw std::runtime_error("Failed to update chunks on EndStroke.\n");
	}
	else {
		int x_chunks = sparse_g->GetXChunks();
		int y_chunks = sparse_g->GetYChunks();
		bool* h_chunk_updates = sparse_g->GetHChunkUpdates();
		for (int chunk_num = 0; ret && (chunk_num < (x_chunks * y_chunks)); ++chunk_num)
		{
			if (h_chunk_updates[chunk_num]) // This chunk was identified as one with updates.
			{
				int x_offset = (chunk_num % x_chunks) * chunk_size;
				int y_offset = (chunk_num / x_chunks) * chunk_size;
				int portion_width = std::min(chunk_size, width - x_offset);
				int portion_height = std::min(chunk_size, height - y_offset);
				ret = ret && AddFloatArrayPortion(full_g, width, height, x_offset, y_offset, portion_width, portion_height, sparse_g->GetChunk(chunk_num), chunk_size, chunk_size, 0, 0);
			}
		}
	}
	return ret;
}


//bool cuda_Paper::cRelaxDivergence(bool** hM, float** hu, float** hv, float** hp, std::set<int> M_chunks, int x_chunks)
//{
//	// Arguments are the host versions (hence the "h") of M (mask), u (horizontal velocities), v (vertical velocities), and p (pressure).
//
//	bool ret = true;
//	int count = 0;
//	int blocksPerGrid = (w * h + threads_per_block - 1) / threads_per_block;
//
//	ret = ret && TransferFromBoolArray(M, hM, w, h, 0);
//	ret = ret && TransferFromFloatArray(u, hu, w + 1, h, 0);
//	ret = ret && TransferFromFloatArray(v, hv, w, h + 1, 0);
//	ret = ret && TransferFromFloatArray(p, hp, w, h, 0);
//
//	while (ret && (count < relaxation_steps))
//	{
//		if (ret)
//		{
//			process_relax_divergence_step1 <<<blocksPerGrid, threads_per_block >>> (M, u, v, delta_matrix, w, h);
//			cudaDeviceSynchronize();
//			cudaStatus = cudaGetLastError();
//			if (cudaStatus != cudaSuccess)
//			{
//				throw std::runtime_error("Failed to process first step in cuda_Paper::cRelaxDivergence.\n");
//				ret = false;
//			}
//		}
//		if (ret)
//		{
//			process_relax_divergence_step2 <<<blocksPerGrid, threads_per_block >>> (M, u, v, p, delta_matrix, w, h);
//			cudaDeviceSynchronize();
//			cudaStatus = cudaGetLastError();
//			if (cudaStatus != cudaSuccess)
//			{
//				throw std::runtime_error("Failed to process second step in cuda_Paper::cRelaxDivergence.\n");
//				ret = false;
//			}
//		}
//		count++;
//	}
//
//	ret = ret && TransferToFloatArray(hu, u, w + 1, h, 0);
//	ret = ret && TransferToFloatArray(hv, v, w, h + 1, 0);
//	ret = ret && TransferToFloatArray(hp, p, w, h, 0);
//	return ret;
//}