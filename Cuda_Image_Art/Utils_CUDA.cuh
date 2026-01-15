#pragma once

// Copyright (c) 2023-2025 Steve Young
// Licensed under the MIT License

#include <iostream>
#include <vector>
#include <set>
#define _USE_MATH_DEFINES
#include <math.h>
#include <array>

#include "..\stb_image_write.h"
#include "..\stb_image.h"

#include <cuda_runtime.h>
#include "device_launch_parameters.h"
#define threads_per_block 256
#define block_dimension 16  // Must be square root of threads_per_block
#define chunk_size 96  // This is the size for chunking out the matrices for efficiency.  * MUST BE A MULTIPLE OF SQUARE ROOT OF 'threads_per_block' * .

__global__ void initialize_Float_Array(float* array, float value, int N);
__global__ void initialize_Int_Array(int* array, int value, int N);
__global__ void initialize_UChar_Array(unsigned char* array, unsigned char value, int N);
__global__ void renormalize_Float_Array(float* array, int N, float min_value, float max_value);
__global__ void initialize_Bool_Array(bool* array, bool value, int N);
__global__ void test_Bool_Array(bool* array, int w, int h, int wx, int wy, bool* result);
__global__ void process_Copy_Unsigned_Char_Array(unsigned char* source, int N, unsigned char* target);
__global__ void process_Copy_Int_Array(int* source, int N, int* target);
__global__ void process_add_partial_Float_Array(float* source, int source_width, int source_height, float* target, int target_width, int target_height, int target_x_offset, int target_y_offset);
__global__ void process_copy_partial_Float_Array(float* source, int source_width, int source_height, float* target, int target_width, int target_height, int target_x_offset, int target_y_offset);
__global__ void process_copy_Float_Array_portion(float* source, int source_image_width, int source_image_height, int source_x_offset, int source_y_offset, int portion_width, int portion_height, float* target, int target_width, int target_height, int target_x_offset, int target_y_offset);
__global__ void process_add_Float_Array_portion(float* source, int source_image_width, int source_image_height, int source_x_offset, int source_y_offset, int portion_width, int portion_height, float* target, int target_width, int target_height, int target_x_offset, int target_y_offset);
__global__ void initialize_gauss_kernel(float* matrix, int g_radius, float thin_factor);
__global__ void set_single_value(float* matrix, int w, int h, int x, int y, float value, bool sum);



class SparseFloatMatrix
{
private:
	int w, h; // Overall width and height of the matrix.
	float background_value; // Value for any element that has not otherwise been initialized.
	int x_chunks, y_chunks; // The number of chunks in the x and y direction (each of size chunk_size).
	std::set<int> chunk_set; // The list of all chunks that have been allocated.
	float** chunk; // Pointers to each 2D float matrix.
	// Temp variables.
	bool* b_result; // Device single-value array for returning a single bool result.
	bool* d_chunk_updates = NULL; // Device vector of chunks to be updated.
	bool* h_chunk_updates = NULL; // Host vector of chunks to be updated.
public:
	SparseFloatMatrix(int width, int height, float value = 0.0f);
	~SparseFloatMatrix();
	bool Set_Value(int x, int y, float value, bool add = false); // *** Try to eliminate this. ***
	bool Reset_Values(bool free_memory = false);  // Argument determines whether memory is actually freed and chunk_set is cleared, or values are just reset to background_value.
	float Get_Value(int x, int y); // *** Try to eliminate this. ***
	bool Copy(SparseFloatMatrix* src);
	int GetChunkNumber(int x, int y);
	bool CheckChunkNumber(int chunk_index);
	bool Check(int x, int y);
	std::set<int> GetChunkSet();
	int GetXChunks();
	int GetYChunks();
	int GetWidth();
	int GetHeight();
	float* GetChunk(int chunk_index);
	bool TestMask(bool* M, int chunk_index);
	float GetBackround();
	bool* GetDChunkUpdates();
	bool* GetHChunkUpdates();
	bool UpdateChunks(); // Use the d_chunk_updates to allocate any additional chunks that are needed.
	bool ExpandToFloatArray(float* d_dst); // Write the data to a FloatArray that has already been allocated on the device.
	bool CompressFromFloatArray(float* d_src); // Uses the FloatArray to update values, and uses d_chunk_updates to determine what new chunks are needed.
	bool SyncChunks(SparseFloatMatrix* src); // Allocate any chunks that are not present but are present in the src matrix.
};

bool WriteOutSparseFloatMatrix(SparseFloatMatrix* source, int x, int y, std::string name, float min, float max);

// Matrix operations for gaussian kernels.
float* GaussKernel(int gauss_radius, float gauss_thin_factor);
float GaussTotal(float* gauss_kernel, int gauss_radius);
bool FreeGaussKernel(float* kernel);

// Matrix operations for floating point arrays.
float* FloatArray(int x, int y, bool initialize_value, float value = 0.0f);
bool RenormalizeFloatArray(float* matrix, int x, int y, float min_value, float max_value);
bool FreeFloatArray(float* a);
bool CopyFloatArray(float* source, float* target, int x, int y);
bool AddPartialFloatArray(float* source, int source_width, int source_height, float* target, int target_width, int target_height, int target_x_offset, int target_y_offset);
bool CopyPartialFloatArray(float* source, int source_width, int source_height, float* target, int target_width, int target_height, int target_x_offset, int target_y_offset);
bool CopyFloatArrayPortion(float* source, int source_image_width, int source_image_height, int source_x_offset, int source_y_offset, int portion_width, int portion_height, float* target, int target_width, int target_height, int target_x_offset, int target_y_offset);
bool AddFloatArrayPortion(float* source, int source_image_width, int source_image_height, int source_x_offset, int source_y_offset, int portion_width, int portion_height, float* target, int target_width, int target_height, int target_x_offset, int target_y_offset);
bool ResetFloatArray(float* matrix, int x, int y, float value);
bool WriteOutFloatArray(float* source, int x, int y, std::string name, float min, float max);
bool CopyFromHost(float* source, int N, float* dest);
bool CopyToHost(float* source, int N, float* dest);
bool PartialCopyToHost(float* source, int source_width, int source_height, int source_x_offset, int source_y_offset, float* dest, int dest_width, int dest_height, int dest_x_offset, int dest_y_offset, int copy_width, int copy_height);

// Matrix operations for boolean arrays.
bool* BoolArray(int x, int y, bool value);
bool FreeBoolArray(bool* a);
bool ResetBoolArray(bool* matrix, int x, int y, bool value);
bool CopyFromHost(bool* source, int N, bool* dest);
bool CopyToHost(bool* source, int N, bool* dest);

// Matrix operations for integer arrays.
int* IntArray(int x, int y, bool initialize_value, int value = 0);
bool FreeIntArray(int* a);
bool ResetIntArray(int* a, int w, int h, int value = 0);
bool CopyFromHost(int* source, int N, int* dest);
bool CopyToHost(int* source, int N, int* dest);
bool OnDeviceCopy(int* source, int N, int* dest);

// Matrix operations for unsigned char arrays.
unsigned char* UCharArray(int x, int y, bool initialize_value, unsigned char value = 0);
bool FreeUCharArray(unsigned char* a);
bool ResetUCharArray(unsigned char* a, int w, int h, unsigned char value = 0);
bool CopyFromHost(unsigned char* source, int N, unsigned char* dest);
bool CopyToHost(unsigned char* source, int N, unsigned char* dest);
bool OnDeviceCopy(unsigned char* source, int N, unsigned char* dest);