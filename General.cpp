#include "General.h"

float arccoth(float x) {
	if (std::abs(x) <= 1) {
		throw std::domain_error("Attempted to call arccoth on a value with an absolute value less than 1.\n");
		return 0.0f;
	}
	return 0.5 * std::log((x + 1.0f) / (x - 1.0f));
}

float** GaussKernel(int gauss_radius, float gauss_thin_factor)
{
	// Gaussian smoothing.  
	// i, j from 0 to +- gauss_radius;
	// gauss_kernel[i][j] = exp(-sqrt(i * i + j * j) * gauss_thin_factor);
	float** gauss_kernel = (float**)malloc(sizeof(float*) * (gauss_radius + 1));
	if (NULL == gauss_kernel)
	{
		throw std::runtime_error("Failed to allocate memory for gauss_kernel for GaussKernel.\n");
		return NULL;
	}
	for (int i = 0; i <= gauss_radius; ++i)
	{
#ifdef _WIN32
		gauss_kernel[i] = static_cast<float*>(_aligned_malloc(sizeof(float) * (gauss_radius + 1), align_bits));
#else
		gauss_kernel[i] = static_cast<float*>(aligned_alloc(align_bits, sizeof(float) * (gauss_radius + 1)));
#endif
		if (NULL == gauss_kernel[i])
		{
			throw std::runtime_error("Failed to allocate memory for gauss_kernel[i] for GaussKernel.\n");
			return NULL;
		}
		for (int j = 0; j <= gauss_radius; ++j)
		{
			gauss_kernel[i][j] = exp(-sqrt(i * i + j * j) * gauss_thin_factor);
		}
	}
	return gauss_kernel;
}

float GaussTotal(float** gauss_kernel, int gauss_radius)
{
	float gauss_total = 0.0;
	for (int i = 0; i <= gauss_radius; ++i)
	{
		for (int j = 0; j <= gauss_radius; ++j)
		{
			if ((i == 0) || (j == 0))
			{
				if ((i == 0) && (j == 0))
				{
					gauss_total += gauss_kernel[i][j];  // Special case at the origin.  There is only one instance of this.
				}
				else {
					gauss_total += 2.0 * gauss_kernel[i][j]; // Case where either i or j (but not both) are on the axis.  Account for the instance directly across from the origin.
				}
			}
			else {
				gauss_total += 4.0 * gauss_kernel[i][j]; // General case.  Account for all four quadrants.
			}
		}
	}
	if (gauss_total <= 0.0f)
	{
		gauss_total = 1.0f;
	}
	return gauss_total;
}

bool FreeGaussKernel(float** kernel, int gauss_radius)
{
	bool ret = true;
	for (int i = 0; i <= gauss_radius; ++i)
	{
#ifdef _WIN32
		_aligned_free(kernel[i]);
#else
		free(kernel[i]);
#endif
	}
	free(kernel);
	return ret;
}

float** FloatArray(int x, int y, bool initialize_value, float value)
{
	float** ret = NULL;
	y = ((y + AVX2_stride - 1) / AVX2_stride) * AVX2_stride; // Increase y as needed to ensure each column is a multiple of align_bits.
	ret = (float**)malloc(sizeof(float*) * x);
	if (NULL == ret)
	{
		throw std::runtime_error("Failed to allocate memory in FloatArray.\n");
		return NULL;
	}
	for (int i = 0; i < x; ++i)
	{
#ifdef _WIN32
		ret[i] = static_cast<float*>(_aligned_malloc(sizeof(float) * y, align_bits));
#else
		ret[i] = static_cast<float*>(aligned_alloc(align_bits, sizeof(float) * y));
#endif
		if (NULL == ret[i])
		{
			throw std::runtime_error("Failed to allocate memory in FloatArray.\n");
			return NULL;
		}
		if (initialize_value)
		{
			for (int j = 0; j < y; ++j)
			{
				ret[i][j] = value;
			}
		}
	}
	return ret;
}

bool RenormalizeFloatArray(float** matrix, int x, int y, float min_value, float max_value)
{
	bool ret = true;
	float low_val = matrix[0][0];
	float high_val = matrix[0][0];
	for (int i = 0; i < x; ++i)
	{
		for (int j = 0; j < y; ++j)
		{
			float val = matrix[i][j];
			if (val > high_val)
			{
				high_val = val;
			}
			if (val < low_val)
			{
				low_val = val;
			}
		}
	}
	float old_range = high_val - low_val;
	float new_range = max_value - min_value;
	for (int i = 0; i < x; ++i)
	{
		for (int j = 0; j < y; ++j)
		{
			matrix[i][j] = min_value + (matrix[i][j] - low_val) * new_range / old_range;
		}
	}
	return ret;
}

bool** BoolArray(int x, int y, bool value)
{
	bool** ret = NULL;
	y = ((y + AVX2_stride - 1) / AVX2_stride) * AVX2_stride; // Increase y as needed to ensure each column is a multiple of align_bits.

	ret = (bool**)malloc(sizeof(bool*) * x);
	if (NULL == ret)
	{
		throw std::runtime_error("Failed to allocate memory in BoolArray.\n");
		return NULL;
	}
	for (int i = 0; i < x; ++i)
	{
		//ret[i] = (bool*)malloc(sizeof(bool) * y);
#ifdef _WIN32
		ret[i] = static_cast<bool*>(_aligned_malloc(sizeof(bool) * y, align_bits));
#else
		ret[i] = static_cast<bool*>(aligned_alloc(align_bits, sizeof(bool) * y));
#endif
		if (NULL == ret[i])
		{
			throw std::runtime_error("Failed to allocate memory in BoolArray.\n");
			return NULL;
		}
		for (int j = 0; j < y; ++j)
		{
			ret[i][j] = value;
		}
	}
	return ret;
}

bool FreeFloatArray(float** a, int x)
{
	if (NULL == a)
	{
		return false;
	}
	for (int i = 0; i < x; ++i)
	{
		if (NULL != a[i])
		{
#ifdef _WIN32
			_aligned_free(a[i]);
#else
			free(a[i]);
#endif
		}
	}
	free(a);
	return true;
}

bool FreeBoolArray(bool** a, int x)
{
	if (NULL == a)
	{
		return false;
	}
	for (int i = 0; i < x; ++i)
	{
		if (NULL != a[i])
		{
			_aligned_free(a[i]);
			//#ifdef _WIN32
			//			_aligned_free(a[i]);
			//#else
			//			free(a[i]);
			//#endif
		}
	}
	free(a);
	return true;
}

bool CopyFloatArray(float** source, float** target, int x, int y)
{
	if ((NULL == source) || (NULL == target))
	{
		throw std::runtime_error("NULL source or target passed to CopyFloatArray.\n");
		return false;
	}
	y = ((y + AVX2_stride - 1) / AVX2_stride) * AVX2_stride; // Increase y as needed to ensure each column is a multiple of align_bits.
	if (use_AVX2)
	{
		for (int i = 0; i < x; ++i)
		{
			for (int j = 0; j < y; j += AVX2_stride)
			{
				__m256 data = _mm256_load_ps(&source[i][j]);
				_mm256_store_ps(&target[i][j], data);
			}
		}
	}
	else {
		for (int i = 0; i < x; ++i)
		{
			for (int j = 0; j < y; ++j)
			{
				target[i][j] = source[i][j];
			}
		}
	}
	return true;
}

bool AddPartialFloatArray(float** source, int source_width, int source_height, float** target, int x_offset, int y_offset)
{
	bool ret = true;
	int remainder_y = source_height % AVX2_stride;
	if (use_AVX2)
	{
		for (int i = 0; i < source_width; ++i)
		{
			for (int j = 0; j < source_height - remainder_y; j += AVX2_stride)
			{
				__m256 source_data = _mm256_load_ps(&source[i][j]);
				__m256 target_data = _mm256_load_ps(&target[i + x_offset][j + y_offset]);
				target_data = _mm256_add_ps(source_data, target_data);
				_mm256_store_ps(&target[i + x_offset][j + y_offset], target_data);
			}
			for (int j = source_height - remainder_y; j < source_height; ++j)
			{
				target[x_offset + i][y_offset + j] += source[i][j];
			}
		}
	}
	else {
		for (int i = 0; i < source_width; ++i)
		{
			for (int j = 0; j < source_height; ++j)
			{
				target[x_offset + i][y_offset + j] += source[i][j];
			}
		}
	}
	return true;
}

bool CopyPartialFloatArray(float** source, int source_x_offset, int source_y_offset, int source_width, int source_height, float** target, int target_x_offset, int target_y_offset)
{
	bool ret = true;
	int remainder_y = source_height % AVX2_stride;
	if (use_AVX2)
	{
		int source_i = source_x_offset;
		for (int i = 0; i < source_width; ++i)
		{
			int source_j = source_y_offset;
			for (int j = 0; j < source_height - remainder_y; j += AVX2_stride)
			{
				__m256 source_data = _mm256_load_ps(&source[source_i][source_j]);
				_mm256_store_ps(&target[i + target_x_offset][j + target_y_offset], source_data);
				source_j += AVX2_stride;
			}
			source_j = source_y_offset + source_height - remainder_y;
			for (int j = source_height - remainder_y; j < source_height; ++j)
			{
				target[target_x_offset + i][target_y_offset + j] = source[source_i][source_j];
				++source_j;
			}
			++source_i;
		}
	}
	else {
		int source_i = source_x_offset;
		for (int i = 0; i < source_width; ++i)
		{
			int source_j = source_y_offset;
			for (int j = 0; j < source_height; ++j)
			{
				target[target_x_offset + i][target_y_offset + j] = source[source_i][source_j];
				++source_j;
			}
			++source_i;
		}
	}
	return true;
}

bool ResetFloatArray(float** matrix, int x, int y, float value)
{
	bool ret = true;
	y = ((y + AVX2_stride - 1) / AVX2_stride) * AVX2_stride; // Increase y as needed to ensure each column is a multiple of align_bits.
	if (use_AVX2)
	{
		__m256 value_vector = _mm256_set1_ps(value);
		for (int i = 0; i < x; ++i)
		{
			for (int j = 0; j < y; j += AVX2_stride)
			{
				_mm256_store_ps(&matrix[i][j], value_vector);
			}
		}
	}
	else {
		for (int i = 0; i < x; ++i)
		{
			for (int j = 0; j < y; ++j)
			{
				matrix[i][j] = value;
			}
		}
	}
	return true;
}

void ResetFloatArrayAVX2(float** matrix, int x, int y, __m256* value_vector)
{
	//y = ((y + AVX2_stride - 1) / AVX2_stride) * AVX2_stride; // Increase y as needed to ensure each column is a multiple of align_bits.   *** Is this needed? ***

	for (int i = 0; i < x; ++i)
	{
		for (int j = 0; j < y; j += AVX2_stride)
		{
			_mm256_store_ps(&matrix[i][j], *value_vector);
		}
	}
	return;
}

bool WriteOutFloatArray(float** source, int x, int y, std::string name, float min, float max)
{
	bool ret = true;
	unsigned char* data;
	float delta = std::abs(max - min);
	data = (unsigned char*)malloc(sizeof(unsigned char) * x * y * 3);
	if (NULL == data)
	{
		throw std::runtime_error("Failed to allocate data for image.\n");
		ret = false;
	}
	if (ret)
	{
		for (int j = 0; j < y; ++j)  // Using y direction on outside to help with caching on writes to data.
		{
			for (int i = 0; i < x; ++i)
			{
				long pos = 3 * (j * x + i);
				unsigned char value;
				if (source[i][j] > max)
				{
					value = 255;
				}
				else if (source[i][j] < min)
				{
					value = 0;
				}
				else {
					value = 255 * ((source[i][j] - min) / delta);
				}
				//if (source[i][j] < 0.0f)
				//{
				//	data[pos] = value;
				//	data[pos + 1] = 0;
				//}else{
				//	data[pos] = 0;
				//	data[pos + 1] = value;
				//}		
				data[pos] = value;
				data[pos + 1] = value;
				data[pos + 2] = value;
			}
		}
		if (0 == stbi_write_png(name.c_str(), x, y, 3, data, x * 3))
		{
			throw std::runtime_error("Unable to write out debug image.\n");
			ret = false;
		}
	}
	free(data);
	return ret;
}

bool WriteOutSparseFloatMatrix(SparseFloatMatrix* source, int x, int y, std::string name, float min, float max)
{
	bool ret = true;
	float** matrix = FloatArray(x, y, false);
	if (NULL == matrix)
	{
		throw std::runtime_error("Failed to allocate matrix for WriteOutSparseFloatMatrix.\n");
		ret = false;
	}
	if (ret)
	{
		for (int i = 0; i < x; ++i)
		{
			for (int j = 0; j < y; ++j)
			{
				matrix[i][j] = source->Get_Value(i, j);
			}
		}
		ret = WriteOutFloatArray(matrix, x, y, name, min, max);
	}
	FreeFloatArray(matrix, x);
	return ret;
}

float* PadMatrix(float** input, int x, int y, int min_pad, int& pad_width, int& pad_height, float background)
{
	float* output = NULL;  // Data is stored by columns, so the entire height of the first column comes before the second column.
	pad_width = x + 2 * min_pad;
	pad_height = (y + 2 * min_pad + align_bits - 1) & ~(align_bits - 1);
#ifdef _WIN32
	output = static_cast<float*>(_aligned_malloc(pad_width * pad_height * sizeof(float), align_bits));
#else
	output = static_cast<float*>(aligned_malloc(align_bits, pad_width * pad_height * sizeof(float)));
#endif
	if (NULL == output)
	{
		throw std::runtime_error("Failed to allocate memory in PadMatrix.\n");
		return NULL;
	}
	if (0.0f == background)
	{
		memset(output, 0, pad_width * pad_height * sizeof(float)); // Set values to zero, so padded regions are zero.
	}
	else {
		for (int i = 0; i < pad_width; ++i)
		{
			if ((i < min_pad) || (i >= (min_pad + x)))
			{
				for (int j = 0; j < pad_height; ++j)
				{
					output[i * pad_height + j] = background;
				}
			}
			else {
				for (int j = 0; j < min_pad; ++j)
				{
					output[i * pad_height + j] = background;
				}
				for (int j = min_pad + y; j < pad_height; ++j)
				{
					output[i * pad_height + j] = background;
				}
			}
		}
	}
	for (int i = 0; i < x; ++i)
	{
		for (int j = 0; j < y; ++j)  // *** Todo: Use SIMD to do this more efficiently.
		{
			output[(i + min_pad) * pad_height + j + min_pad] = input[i][j];
		}
	}
	return output;
}

bool FreePadMatrix(float* matrix)
{
	bool ret = true;
#ifdef _WIN32
	_aligned_free(matrix);
#else
	free(matrix);
#endif
	return ret;
}

bool Convolve(const float* input, float* output, float** kernel, const int pad_width, const int pad_height, const int kernel_radius)
{
	// The input and output matrices need to be padded using the PadMatrix function.  input and output are stored by columns.
	// pad_width and pad_height apply to input and output matrices.
	// kernel_size is twice the kernel_radius plus 1 (so, always odd).  The kernel_radius is the radius beyond the origin point.

	bool ret = true;
	int kernel_size = 2 * kernel_radius + 1;
	float* alignedKernel = NULL;
#ifdef _WIN32
	alignedKernel = static_cast<float*>(_aligned_malloc(sizeof(float) * kernel_size * kernel_size, align_bits));
#else
	alignedKernel = static_cast<float*>(aligned_alloc(align_bits, sizeof(float) * kernel_size * kernel_size));
#endif
	for (int k = 0; k < kernel_size * kernel_size; ++k)
	{
		int kx = abs(k / kernel_size - kernel_radius);
		int ky = abs(k % kernel_size - kernel_radius);
		alignedKernel[k] = kernel[kx][ky];
	}
	for (int i = kernel_radius; i < (pad_width - kernel_radius); ++i)
	{
		for (int j = kernel_radius; j < (pad_height - kernel_radius - AVX2_stride); j += AVX2_stride)
		{
			__m256 sum = _mm256_setzero_ps(); // Zero out sum, which is the accumulator.
			for (int kx = 0; kx < kernel_size; ++kx)
			{
				for (int ky = 0; ky < kernel_size; ++ky)
				{
					__m256 inputVals = _mm256_loadu_ps(
						&input[(j + ky - kernel_radius) + (i + kx - kernel_radius) * pad_height]
					);
					__m256 kernelVal = _mm256_broadcast_ss(&alignedKernel[ky + kx * kernel_size]);
					sum = _mm256_add_ps(sum, _mm256_mul_ps(inputVals, kernelVal));
				}
			}
			_mm256_storeu_ps(&output[j + i * pad_height], sum);  // *** Can't I use the _mm256_store_ps version?
		}
		// If there are any more elements to the columns, handle those.
		for (int j = pad_height - kernel_radius - (pad_height - 2 * kernel_radius) % AVX2_stride;
			j < (pad_height - kernel_radius); ++j)
		{
			float sum = 0.0f;
			for (int kx = 0; kx < kernel_size; ++kx)
			{
				for (int ky = 0; ky < kernel_size; ++ky)
				{
					sum += input[(j + ky - kernel_radius) + (i + kx - kernel_radius) * pad_height]
						* alignedKernel[ky + kx * kernel_size];
				}
			}
			output[j + i * pad_height] = sum;
		}
	}
	return ret;
}


SparseFloatMatrix::SparseFloatMatrix(int width, int height, float value)
{
	w = width;
	h = height;
	background_value = value;

	x_chunks = w / chunk_size;
	if (x_chunks * chunk_size < w)
	{
		x_chunks++;
	}

	y_chunks = h / chunk_size;
	if (y_chunks * chunk_size < h)
	{
		y_chunks++;
	}

	chunk_set.clear(); // Start with no allocated chunks.

	chunk = (float***)malloc(sizeof(float**) * x_chunks * y_chunks);
	if (NULL == chunk)
	{
		throw std::runtime_error("Failed to allocate chunk in SparseFloatMatrix.\n");
		return;
	}
	for (int chunk_count = 0; chunk_count < (x_chunks * y_chunks); ++chunk_count)
	{
		chunk[chunk_count] = NULL;
	}
}

SparseFloatMatrix::~SparseFloatMatrix()
{
	std::set<int>::iterator chunk_iterator;
	for (chunk_iterator = chunk_set.begin(); chunk_iterator != chunk_set.end(); ++chunk_iterator)
	{
		int chunk_index = *chunk_iterator;
		FreeFloatArray(chunk[chunk_index], chunk_size);
	}
	free(chunk);
}

bool SparseFloatMatrix::Set_Value(int x, int y, float value, bool add)
{
	bool ret = true;
	int chunk_index = GetChunkNumber(x, y);
	if (chunk_index >= 0)
	{
		int x_offset, y_offset;
		if (false == CheckChunkNumber(chunk_index))
		{
			chunk[chunk_index] = FloatArray(chunk_size, chunk_size, true, background_value);
			chunk_set.insert(chunk_index);
		}
		x_offset = x - (x / chunk_size) * chunk_size;  // *** Calculate offsets once at creation of SparseFloatMatrix, store in vector indexed by chunk number. ***
		y_offset = y - (y / chunk_size) * chunk_size;
		if (add)
		{
			chunk[chunk_index][x_offset][y_offset] += value;
		}
		else {
			chunk[chunk_index][x_offset][y_offset] = value;
		}
	}
	else {
		throw std::runtime_error("Set_Value called with x or y out of bounds.\n");
		ret = false;
	}
	return ret;
}

bool SparseFloatMatrix::Reset_Values()
{
	bool ret = true;
	std::set<int>::iterator chunk_it;
	for (chunk_it = chunk_set.begin(); chunk_it != chunk_set.end(); ++chunk_it)
	{
		int chunk_index = *chunk_it;
		FreeFloatArray(chunk[chunk_index], chunk_size);
		chunk[chunk_index] = NULL;
	}
	chunk_set.clear();
	return ret;
}

float SparseFloatMatrix::Get_Value(int x, int y)
{
	float ret = 0.0f;
	int chunk_index = GetChunkNumber(x, y);
	if ((false == CheckChunkNumber(chunk_index)) || (x < 0) || (x >= w) || (y < 0) || (y >= h))
	{
		ret = background_value;
	}
	else {
		int x_offset = x - (x / chunk_size) * chunk_size;
		int y_offset = y - (y / chunk_size) * chunk_size;
		ret = chunk[chunk_index][x_offset][y_offset];
	}
	return ret;
}

bool SparseFloatMatrix::Copy(SparseFloatMatrix* src)
{
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

				//for (int i = 0; i < chunk_size; ++i)
				//{
				//	for (int j = 0; j < chunk_size; ++j)
				//	{
				//		chunk[chunk_index][i][j] = src->chunk[chunk_index][i][j];  // *** This can be speeded up with AVX2. ***
				//	}
				//}
			}
			else { // The chunk is not present in the source matrix.
				if (dest_chunk)
				{
					FreeFloatArray(chunk[chunk_index], chunk_size);
					chunk_set.erase(chunk_index);
				}
			}
		}
	}
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
	return CheckChunkNumber(GetChunkNumber(x, y));
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

float** SparseFloatMatrix::GetChunk(int chunk_index)
{
	return chunk[chunk_index];
}

bool SparseFloatMatrix::TestMask(bool** M, int chunk_index)
{
	bool ret = false;
	int wy = chunk_index / x_chunks;  // Y component of beginning of window represented by chunk.
	int wx = chunk_index - (wy * x_chunks); // X component of the beginning of windows represented by chunk. 
	wx = wx * chunk_size;
	wy = wy * chunk_size;
	for (int i = wx; (i < (wx + chunk_size)) && (i < w); ++i)
	{
		for (int j = wy; (j < (wy + chunk_size)) && (j < h); ++j)
		{
			if (M[i][j])
			{
				return true;
			}
		}
	}
	return ret;
}

