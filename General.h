#pragma once

// Copyright (c) 2023-2024 Steve Young
// Licensed under the MIT License

#include <stdexcept>
#include <set>
#include <thread>
#include <algorithm>
#include <execution>
#include <immintrin.h>
#include <cmath>
#include "stb_image_write.h"
#include "stb_image.h"

#define EFFECTIVE_ZERO 0.01f // 0.001f
#define CURVE_ADJUSTMENT 4
#define MAX_CURVE_FACTOR 350
#define OUTLINE false
#define BACKGROUND 0
#define GLITCH1 false
#define GLITCH2 false
#define GLITCH3 false
#define PAINT_SCALE 1
#define SCALE_FLOW true
#define BRISTLES 40
#define FLOW 20.0
#define FLOW_VARIATION 10.0
#define SUBPIXEL false
#define BRUSH_WIDTH_FACTOR (float)1.8
#define BRISTLE_THIN_FACTOR 2.0
#define BRUSH_SHAPE shape_straight
#define BRUSH_WIDTH_OVERRIDE false
#define BRUSH_WIDTH 30
#define MIX_PAINTS true
#define RADIUS_VARIATION true
#define DATA_FILE_MAJOR 0
#define DATA_FILE_MINOR 2
#define FINGERPRINT "IMRT"
#define CONTRAST_RADIUS 0
#define CONTRAST_BOX_MARGIN 2
#define PAINT_MASK true
#define WATERCOLOR true

// Watercolor defines

#define mu 0.03f // 0.005 - 0.05f
#define kappa 0.001f // 0.001f - 0.040f
#define absorption_alpha 0.1f // 0.1f
#define saturation_sigma 0.1f // 0.1f
#define saturation_dry_value 0.005f
#define saturation_epsilon 0.48f  // 0.48 0.53 0.55 0.57 
#define saturation_delta -0.1f // -0.1f
#define capacity_max 0.9f // 0.9f
#define capacity_min 0.3f // 0.3f
#define relaxation_steps 70 // 50, 120
#define tau 0.005f // 0.01f, 0.005f
#define xi 0.15f // 0.1f
#define K_radius 5 // 5
#define eta 0.001f // 0.03f
#define slope_factor 1.0f // 1.0f
#define slope_velocity_steps 400
#define max_velocity 2.5f // 2.5f
#define render_g_value 1.0f // render_g_value 0.2 -0.035f
#define process_steps 400 // 1400
// Set pressure_delta to -1 for speckled image with more water dispertion
#define pressure_delta 1
#define dab_pressure 0.00005f
#define chunk_size 96  // This is the size for chunking out the matrices for efficiency.  * MUST BE A MULTIPLE OF AVX2_STRIDE * .

// Values for AVX2
#define AVX2_stride 8
#define align_bits 32
#define use_AVX2 true

struct RectQuad {
	int x0, y0, x1, y1;
};

struct PointPair {
	int x, y;
};

struct FloatPointPair {
	float x, y;
};

struct Color {
	float channel[3];
};

struct Corner
{
	// p0 is the first midpoint between vertices.
	// p1 is the second midpoint.
	// c0 is the first control point, unless this is a corner, in which case it is the corner point.
	// c1 is the second control point, if there is one.
	// smooth indicates whether this is a curve (using control points) or a corner.
	bool smooth;
	FloatPointPair p0; // Midpoints
	FloatPointPair p1;
	FloatPointPair c0; // Control Points
	FloatPointPair c1;
	float radius_p0; // Radius needed to touch edge of region.  Can be calculated different ways.
	float radius_p1;
	float radius_c0;
};

enum brush_shape { shape_round, shape_straight, shape_test };

struct Paint_Properties
{
	bool outline = OUTLINE;
	int background = BACKGROUND;  // 0 - white, 1 - black
	bool glitch1 = GLITCH1;
	bool glitch2 = GLITCH2;
	bool glitch3 = GLITCH3;
	float paint_scale = PAINT_SCALE;
	float bristles = BRISTLES;
	float flow = FLOW;
	float flow_variation = FLOW_VARIATION;
	bool sub_pixel = SUBPIXEL;
	float max_curve_factor = MAX_CURVE_FACTOR;
	float brush_width_factor = BRUSH_WIDTH_FACTOR;
	float bristle_thin_factor = BRISTLE_THIN_FACTOR;
	brush_shape shape = BRUSH_SHAPE;
	bool brush_width_override = BRUSH_WIDTH_OVERRIDE;
	float brush_width = BRUSH_WIDTH;
	bool mix_paints = MIX_PAINTS;
	bool radius_variation = RADIUS_VARIATION;
	bool paint_mask = PAINT_MASK;
	bool extinguish_paint = false;
	bool watercolor = WATERCOLOR;
};

float arccoth(float x);

class SparseFloatMatrix
{
private:
	int w, h; // Overall width and height of the matrix.
	float background_value; // Value for any element that has not otherwise been initialized.
	int x_chunks, y_chunks; // The number of chunks in the x and y direction (each of size chunk_size).
	std::set<int> chunk_set; // The list of all chunks that have been allocated.
	float*** chunk; // Pointers to each 2D float matrix.
public:
	SparseFloatMatrix(int width, int height, float value = 0.0f);
	~SparseFloatMatrix();
	bool Set_Value(int x, int y, float value, bool add = false); // Sets value in matrix, unless add is true in which case the value is added to existing value.
	bool Reset_Values(); // Reset to background, by removing allocated chunks.
	float Get_Value(int x, int y);
	bool Copy(SparseFloatMatrix* src);
	int GetChunkNumber(int x, int y);
	bool CheckChunkNumber(int chunk_index);
	bool Check(int x, int y);
	std::set<int> GetChunkSet();
	int GetXChunks();
	int GetYChunks();
	float** GetChunk(int chunk_index);
	bool TestMask(bool** M, int chunk_index); // Return true if any element of M[x][y] is true within the identified chunk.
};
bool WriteOutSparseFloatMatrix(SparseFloatMatrix* source, int x, int y, std::string name, float min, float max);

// Matrix operations for gaussian kernels.
float** GaussKernel(int gauss_radius, float gauss_thin_factor);
float GaussTotal(float** gauss_kernel, int gauss_radius);
bool FreeGaussKernel(float** kernel, int gauss_radius);

// Matrix operations for floating point arrays that are 32-bit aligned in the y direction for AVX2 optimization.
float** FloatArray(int x, int y, bool initialize_value, float value = 0.0f);
bool RenormalizeFloatArray(float** matrix, int x, int y, float min_value, float max_value);
bool FreeFloatArray(float** a, int x);
bool CopyFloatArray(float** source, float** target, int x, int y);
bool AddPartialFloatArray(float** source, int source_width, int source_height, float** target, int x_offset, int y_offset);  // Add source array to sub-section of target array, at offset.
bool CopyPartialFloatArray(float** source, int source_x_offset, int source_y_offset, int source_width, int source_height, float** target, int target_x_offset, int target_y_offset);  // Copy source array to sub-section of target array, at offset.
bool ResetFloatArray(float** matrix, int x, int y, float value);
void ResetFloatArrayAVX2(float** matrix, int x, int y, __m256* value_vector);
bool WriteOutFloatArray(float** source, int x, int y, std::string name, float min, float max);

// Matrix operations for boolean arrays that are 32-bit aligned in the y direction for AVX2 optimization.
bool** BoolArray(int x, int y, bool value);
bool FreeBoolArray(bool** a, int x);

bool Convolve(const float* input, float* output, float** kernel, const int pad_width, const int pad_height, const int kernel_radius);  // AVX2 optimized convolution function.
float* PadMatrix(float** input, int x, int y, int min_pad, int& pad_width, int& pad_height, float background = 0.0f);  // Function to create a matrix for Convolve function.
bool FreePadMatrix(float* matrix);  // Free the padded matrix used for the Convolve function.