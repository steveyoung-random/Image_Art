#pragma once

// Copyright (c) 2023-2024 Steve Young
// Licensed under the MIT License

#define EFFECTIVE_ZERO 0.001
#define CURVE_ADJUSTMENT 4
#define MAX_CURVE_FACTOR 80
#define OUTLINE false;
#define BACKGROUND 0;
#define GLITCH1 false;
#define GLITCH2 false;
#define GLITCH3 false;
#define PAINT_SCALE 1;
#define SCALE_FLOW true;
#define BRISTLES 40;
#define FLOW 20.0;
#define FLOW_VARIATION 10.0;
#define SUBPIXEL false;
#define BRUSH_WIDTH_FACTOR 1.8;
#define BRISTLE_THIN_FACTOR 2.0;
#define BRUSH_SHAPE shape_straight;
#define BRUSH_WIDTH_OVERRIDE false;
#define BRUSH_WIDTH 30;
#define MIX_PAINTS true;
#define RADIUS_VARIATION true;

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
};