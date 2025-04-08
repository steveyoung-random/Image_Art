#pragma once

#include <cuda_runtime.h>
#include <vector>
#include "brush.h"
#include "device_launch_parameters.h"
#define threads_per_block 256
#define block_dimension 16  // Must be square root of threads_per_block


__global__ void initialize_Float_Array(float* array, float value, int N);
__global__ void initialize_Int_Array(int* array, int value, int N);
__global__ void renormalize_Float_Array(float* array, int N, float min_value, float max_value);
__global__ void initialize_Bool_Array(bool* array, bool value, int N);
__global__ void test_Bool_Array(bool* array, int w, int h, int wx, int wy, bool* result);
__global__ void process_add_partial_Float_Array(float* source, int source_width, int source_height, float* target, int target_width, int target_height, int target_x_offset, int target_y_offset);
__global__ void process_copy_partial_Float_Array(float* source, int source_width, int source_height, float* target, int target_width, int target_height, int target_x_offset, int target_y_offset);
__global__ void process_copy_Float_Array_portion(float* source, int source_image_width, int source_image_height, int source_x_offset, int source_y_offset, int portion_width, int portion_height, float* target, int target_width, int target_height, int target_x_offset, int target_y_offset);
__global__ void process_add_Float_Array_portion(float* source, int source_image_width, int source_image_height, int source_x_offset, int source_y_offset, int portion_width, int portion_height, float* target, int target_width, int target_height, int target_x_offset, int target_y_offset);
__global__ void process_relax_divergence_step1(bool* M, float* u, float* v, float* delta_matrix, int w, int h);
__global__ void process_relax_divergence_step2(bool* M, float* u, float* v, float* p, float* delta_matrix, int w, int h);
__global__ void process_relax_divergence_step2_b(bool* M, float* u, float* v, float* p, float* delta_matrix, int w, int h);
__global__ void process_convolution(float* C_input, float* C_output, float* C_kernel, float scale, int width, int height, int padding, int kernel_radius);
__global__ void initialize_gauss_kernel(float* matrix, int g_radius, float thin_factor);
__global__ void process_slope_velocities(float* thickness, float* u, float* v, int w, int h, int padding);
__global__ void process_max_velocity(float* u, float* v, float* results, int w, int h);
__global__ void process_max_velocity_2(int w, int h, float* u, float* v, float* result);
__global__ void get_max_results(float* results);
__global__ void get_max_results_2(float* input, int input_length, float* output);
__global__ void process_calc_velocities(float dt, int w, int h, float* u, float* v, float* u_prime, float* v_prime, float* p);
__global__ void process_enforce_boundaries(int w, int h, float* u, float* v, bool* M);
__global__ void process_calc_Mprime(bool* M, float* Mprime, float* Mprime_Kernel, int k_radius, float MprimeKernelTotal, int w, int h, int x0, int y0, int x1, int y1);
__global__ void process_flow_outward(bool* M, float* Mprime, float* p, int w, int h);
__global__ void process_active_chunks(bool* M, bool* M_chunks, int x_chunks, int blocks_per_chunk, int w, int h);
__global__ void process_dab_step1(bool* M, bool* chunk_updates, float* s, float* p, int w, int h, int x, int y, int wx, int wy, int w_width, int w_height, int radius_squared, float saturation, float pressure);
__global__ void process_dab_step2(float* chunk, float concentration, int w, int h, int x, int y, int wx, int wy, int w_width, int w_height, int chunk_x0, int chunk_y0, int radius_squared);
__global__ void process_paint_step1(bool* M, bool* chunk_updates, float* s, float* p, int w, int h, int* SPdata, int identifier, int wx, int wy, int w_width, int w_height, float saturation, float pressure);
__global__ void process_paint_step2(float* chunk, float concentration, int w, int h, int* SPdata, int identifier, int wx, int wy, int w_width, int w_height, int chunk_x0, int chunk_y0);
__global__ void process_bristle_dab(cudaBrush* brush, int num_bristles, FloatPointPair direction, int* mask, int mask_width, int mask_height, int mask_value, float scaled_spot_radius_squared, bool begin, float paint_scale, float saturation, float pressure, float concentration, bool* chunk_updates);
__global__ void set_single_value(float* matrix, int w, int h, int x, int y, float value, bool sum);
__global__ void calc_render_step1(float* reflected, float* upward_transmission, float* downward_transmission, float* g_chunk, float* d_chunk, float K, float S);
__global__ void calc_render_step2(float* reflected, float* upward_transmission, float* downward_transmission, float substrate_color, unsigned char* d_image, int w, int h, int x0, int y0, int channel);
__global__ void calc_move_pigment(int w, int h, bool* M, int chunk_num, float* u, float* v, float* g, float* g_prime, bool* d_chunk_updates, float step_size);
__global__ void process_transfer_pigment(int w, int h, int chunk_index, float* chunk_g, float* chunk_d, float* thickness, int padding, float gamma, float rho, float omega, int x_chunks, bool* M, bool* M_chunks);
__global__ void process_capillary_flow_step1(int w, int h, float* s, float* thickness, int padding, bool* M, bool* M_chunks);
__global__ void process_capillary_flow_step2(int w, int h, float* s, float* s_prime, float* thickness, int padding, bool* M, bool* M_chunks);
__global__ void process_capillary_flow_step3(int w, int h, float* s, float* s_prime, float* p, float* thickness, int padding, bool* M, bool* M_chunks);
__global__ void process_dry(float* g_chunk, float* d_chunk);
__global__ void process_brush_location(cudaBrush* brush, FloatPointPair loc);
__global__ void process_brush_orientation(cudaBrush* brush, float orientation);

class Pigment;

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

struct cudaBristle
{
	FloatPointPair offset;
	FloatPointPair wander;
	FloatPointPair last_loc;
	float flow_difference;
	bool down;
};

struct cudaBrush
{
	brush_shape shape;
	float orientation; // Radians. Indicates direction of movement of brush (default is moving towards x=1.0, y=0.0).  Positive moves a point at (1.0, 0.0) in the diretion of positive Y.
	float orient_cosine; // Cosine value for orientation.
	float orient_sine;  // Sine value for orientation.
	float brush_width;
	float brush_depth;
	int num_bristles;
	FloatPointPair location;
	Color color;
	Color second;
	cudaBristle* bristles;
	float* bristle_kernel;
	Paint_Properties paint_prop;
	bool watercolor;
	Paper* watercolor_paper;  // Is this needed?
	int watercolor_pigment_index;
	// Paper or pigment values are stored below, to avoid having to send these repeatedly to the device.
	bool* M;
	float* s;
	float* p;
	SparseFloatMatrix* sparse_g;
	float* full_g;
	int width;  // Paper width.
	int height; // Paper height.
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

bool Convolve(float* input, int matrix_width, int matrix_height, int padding, float* output, float* in_kernel, float scale, int k_radius, float background);   
bool SlopeVelocities(float* thickness, float* u, float* v, int w, int h, int padding);
float GetMaxVelocity(float* u, float* v, int w, int h, float* results1, float* results2, int results_length);
bool CalcVelocities(float dt, int w, int h, float* u, float* v, float* u_prime, float* v_prime, float* p, bool* M);
bool CalcRelaxDivergence(bool* M, float* u, float* v, float* p, float* delta_matrix, int w, int h);
bool CalcRelaxDivergence_b(bool* M, float* u, float* v, float* p, float* delta_matrix, int w, int h);
bool CalcMprime(bool* M, float* Mprime, float* Mprime_Kernel, int k_radius, float MprimeKernelTotal, int w, int h, int x0, int y0, int x1, int y1);
bool CalcFlowOutward(bool* M, float* Mprime, float* p, int w, int h);
bool CalcActiveChunks(bool* M, bool* M_chunks, bool* host_M_chunks, int x_chunks, int y_chunks, int w, int h);
bool CalcDab(bool* M, float* s, float* p, SparseFloatMatrix* g, int w, int h, int x, int y, int radius, float saturation, float concentration);
bool CSetVelocity(float* u, float* v, int x, int y, int w, int h, float u_vel, float v_vel, bool sum, bool wait);
bool CRender(std::vector<Pigment*> pigments, unsigned char* data, int w, int h, float* substrate_color);
bool CMovePigment(int w, int h, float* u, float* v, bool* M, bool* host_M_chunks, std::vector<Pigment*> pigments, float dt, float* g, float* g_prime);
bool CTransferPigment(Pigment* pgmnt, float* thickness, int padding, bool* M, bool* M_chunks);
bool CCapillaryFlow(int w, int h, float* s, float* s_prime, float* p, float* thickness, int gauss_radius, bool* M, bool* M_chunks);
bool CDry(SparseFloatMatrix* pgmnt_g, SparseFloatMatrix* pgmnt_d);
bool CPaintArea(bool* M, float* s, float* p, SparseFloatMatrix* g, int* SPdata, int w, int h, RectQuad window, int identifier, float saturation, float concentration);

// Brush related functions
cudaBrush* CreateCudaBrush(cudaBrush* host_brush, std::vector<Bristle*> host_bristles);
bool FreeCudaBrush(cudaBrush* brush, cudaBristle* bristles, float* bristle_kernel);
bool SetCudaBrushLocation(cudaBrush* brush, FloatPointPair loc);
bool SetCudaBrushOrientation(cudaBrush* brush, float orientation);
bool CudaDab3(cudaBrush* brush, SparseFloatMatrix* g, float saturation, float concentration, FloatPointPair direction, int num_bristles, SPixelData* mask, int mask_value, float spot_radius, bool begin, float paint_scale);
bool StartStroke(SparseFloatMatrix* sparse_g, float* full_g, int width, int height);
bool EndStroke(SparseFloatMatrix* sparse_g, float* full_g, int width, int height);
