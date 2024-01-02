#pragma once

// Copyright (c) 2023-2024 Steve Young
// Licensed under the MIT License

#include <iostream>
#include <vector>
#include "General.h"
#define _USE_MATH_DEFINES
#include <math.h>

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
	float orientation; // Radians. Indicates direction of movement of brush (default is moving towards x=1.0, y=0.0).
	float brush_width;
	float brush_depth;
	int num_bristles;
	FloatPointPair location;
	Color color;
	Color second;
	std::vector<Bristle*> bristles;
	float** bristle_kernel;
	Paint_Properties paint_prop;

public:
	Brush(FloatPointPair start, Color c, Color sec, float w, float d, Paint_Properties prop);
	~Brush();

	bool Dab(unsigned char* data, int width, int height);
	bool Dab2(float* data, int width, int height, SPixelData* mask = NULL, int mask_value = 0, float spot_radius = -1, float flow_adjustment = 1.0);
	bool Dab3(FloatPointPair direction, float* data, int width, int height, SPixelData* mask = NULL, int mask_value = 0, float spot_radius = -1, float flow_adjustment = 1.0, bool begin=false);
	bool MoveTo(FloatPointPair loc);
	bool PaintTo(FloatPointPair loc, float* data, int width, int height, SPixelData* mask = NULL, int mask_value = 0, float rad1 = -1, float rad2 = -1);
	bool PaintTo2(FloatPointPair loc2, FloatPointPair o2, float* data, int width, int height, SPixelData* mask = NULL, int mask_value = 0, float rad1 = -1, float rad2 = -1, bool begin=false);
	bool ChangeColor(Color c, Color sec);
	bool PaintCorner(Corner corner, float* data, int width, int height, SPixelData* mask = NULL, int mask_value = 0, bool variable_radius = false);
	bool PaintCorner2(Corner corner, float* data, int width, int height, SPixelData* mask = NULL, int mask_value = 0, bool variable_radius = false);
	bool SetOrientation(FloatPointPair o);
	bool SetOrientation(float o);
	float GetOrientation();
	float KernelAdjustment(int i, int j, float x, float y);
};
