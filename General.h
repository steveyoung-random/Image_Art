#pragma once

// Copyright (c) 2023-2025 Steve Young
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


// CUDA
//#define USE_CUDA

#define EFFECTIVE_ZERO 0.001f // 0.001f
#define CURVE_ADJUSTMENT 4
#define MAX_CURVE_FACTOR 350
#define OUTLINE false
#define BACKGROUND 0
#define GLITCH1 false
#define GLITCH2 false
#define GLITCH3 false
#define PAINT_SCALE 1
#define SCALE_FLOW true
#define BRISTLES 80
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
#define BRISTLE_KERNEL 5
#define DEFAULT_TO_WATERCOLOR true
#ifdef USE_CUDA
#define WATERCOLOR DEFAULT_TO_WATERCOLOR
#else
#define WATERCOLOR false
#endif

// Watercolor defines

#define mu 0.1f // wet: 0.020f dry: 0.05f 0.005 - 0.05f
#define kappa 0.04f // wet: 0.001f dry: 0.001f 0.001f - 0.040f
#define absorption_alpha 0.001f // wet: 0.1f dry: 0.001f 0.1f
#define saturation_sigma 0.1f // wet: 0.1 f dry: 0.1f 0.1f
#define saturation_dry_value 0.005f // wet: 0.005f dry: 0.005f
#define saturation_epsilon 0.999f // wet: 0.9f dry: 0.999f 0.9
#define saturation_max_diffusion 0.1f // wet: 0.1f dry: 0.1f 0.1
#define capacity_max 0.9f // wet: 0.9f dry: 0.9f 0.9f
#define capacity_min 0.1f // wet: 0.1f dry: 0.1f 0.3f
#define relaxation_steps 70 // wet: 70 dry: 70 50, 120
#define tau 0.005f // wet: 0.005f dry: 0.005f 0.01f, 0.005f
#define xi 0.15f // wet: 0.15f dry: 0.15f 0.1f
#define K_radius 5 // wet: 5 dry: 5 5
#define eta 0.05f // wet: 0.001f dry: 0.05f 0.03f
#define pigment_lag 0.85f // wet: 1.0f dry: 0.85f
#define slope_factor 1.0f // wet: 1.0f dry: 1.0f 1.0f
#define slope_velocity_steps 400 // wet: 400 dry: 400
#define max_velocity 2.5f // wet: 2.5f dry: 2.5f 2.5f
#define render_g_value 0.0f // wet: 0.8f dry: 0.45f 0.2 -0.035f
#define process_steps 100 // wet: 400 dry: 400 1400
// Set pressure_delta to -1 for speckled image with more water dispertion
#define pressure_delta 1 // wet: 1 dry: 1
#define dab_pressure 0.01f // wet: 0.01f dry: 0.01f
#define dab_concentration 0.0075f // wet: 0.01f dry:
#define dab_saturation 0.01f // wet: 0.05f dry:
#define chunk_size 96  // This is the size for chunking out the matrices for efficiency.  * MUST BE A MULTIPLE OF SQUARE ROOT OF 'threads_per_block' * .
#define brush_velocity 2.0f // The induced velocity during painting.

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
