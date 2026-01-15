#pragma once

// Copyright (c) 2023-2025 Steve Young
// Licensed under the MIT License

#include "Utils_CUDA.cuh"
#include "..\brush.h"

__global__ void process_relax_divergence_step1(bool* M, float* u, float* v, float* delta_matrix, int w, int h);
__global__ void process_relax_divergence_step2(bool* M, float* u, float* v, float* p, float* delta_matrix, int w, int h);
__global__ void process_relax_divergence_step2_b(bool* M, float* u, float* v, float* p, float* delta_matrix, int w, int h);
__global__ void process_convolution(float* C_input, float* C_output, float* C_kernel, float scale, int width, int height, int padding, int kernel_radius);
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
