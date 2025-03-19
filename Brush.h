#pragma once

// Copyright (c) 2023-2025 Steve Young
// Licensed under the MIT License

#include <iostream>
#include <vector>
#include "General.h"
#include "Paper.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <array>

class SPixelData;

class Bristle {
private:
	FloatPointPair offset;
	FloatPointPair wander;
	FloatPointPair last_loc;
	float flow_difference;
	bool down;
public:
	Bristle(FloatPointPair o, float fd);
	FloatPointPair GetOffset();
	FloatPointPair GetUnadjustedOffset();
	float GetFlowDiff();
	bool AdjustOffset(FloatPointPair o);
	bool AdjustWander(FloatPointPair w);
	bool SetLast(FloatPointPair loc);
	FloatPointPair GetLast();
	bool GetBristleDown();
	bool SetBristleDown(bool d);
};


class Brush
{
private:
	brush_shape shape;
	float orientation; // Radians. Indicates direction of movement of brush (default is moving towards x=1.0, y=0.0).  Positive moves a point at (1.0, 0.0) in the diretion of positive Y.
	float brush_width;
	float brush_depth;
	int num_bristles;
	FloatPointPair location;
	Color color;
	Color second;
	std::vector<Bristle*> bristles;
	float** bristle_kernel;
	Paint_Properties paint_prop;
	bool watercolor;
	Paper* watercolor_paper;
	int watercolor_pigment_index;

public:
	Brush(FloatPointPair start, Color c, Color sec, float w, float d, Paint_Properties prop, Paper* paper=NULL, int pigment_index=-1);
	~Brush();

	bool Dab3(FloatPointPair direction, float* data, int width, int height, SPixelData* mask = NULL, int mask_value = 0, float spot_radius = -1, bool begin=false, SPixelData* extinguish_mask = NULL);
	bool MoveTo(FloatPointPair loc);
	bool PaintTo2(FloatPointPair loc2, FloatPointPair o2, float* data, int width, int height, SPixelData* mask, int mask_value, bool use_mask, float rad1 = -1, float rad2 = -1, bool begin=false, SPixelData* extinguish_mask = NULL);
	bool ChangeColor(Color c, Color sec);
	bool PaintCorner2(Corner corner, float* data, int width, int height, SPixelData* mask, int mask_value, bool use_mask=false, SPixelData* extinguish_mask=NULL);
	bool ExtinguishCorner(Corner corner, float* data, int width, int height, SPixelData* mask, int mask_value);
	bool SetOrientation(FloatPointPair o);
	bool SetOrientation(float o);
	float GetOrientation();
	float CalculateOrientation(FloatPointPair direction);
	float OrientationDifference(float o1, float o2);
	float KernelAdjustment(int i, int j, float x, float y);
	bool ExtinguishQuadrilateral(std::array<FloatPointPair, 4> points, int width, int height, SPixelData* extinguish_mask, int value);
	bool ExtinguishLineSegment(std::array<FloatPointPair, 2> points, int width, int height, SPixelData* extinguish_mask, int value);
};
