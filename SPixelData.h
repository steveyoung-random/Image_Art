#pragma once

// Copyright (c) 2023-2024 Steve Young
// Licensed under the MIT License

#include <iostream>
#include <set>
#include "ImageData.h"

class SuperPixel;

class SPixelData {
private:
	int* data = NULL;
	int width, height;
	std::set<int> meeting_points;

public:
	SPixelData(int w, int h);
	SPixelData(const SPixelData&);
	~SPixelData();
	int GetPixel(int x, int y);
	bool GetMeetingPoint(PointPair point);
	int GetNeighbors(int identifier, PointPair point);
	bool SetMeetingPoint(PointPair point);
	bool SetPixel(int x, int y, int value);
	bool Reset();
	bool CopyData(SPixelData* sp);
	bool erode(int mode = 0, int struct_size = 3);
	bool dilate(SuperPixel* sp, int mode = 0, int struct_size = 3);
	bool dilate_erode(SuperPixel* sp, bool dilate, int mode = 0, int struct_size = 3);
	ImageData* GenerateImage(SuperPixel* sp, Color background);
	int GetWidth();
	int GetHeight();
	float CalculateRadius(std::vector<Corner>::iterator curve_begin, std::vector<Corner>::iterator curve_end, int mask_value);
	float RadiusTransverse(FloatPointPair p, FloatPointPair c, int mask_value);
	FloatPointPair AxialExtent(FloatPointPair p, FloatPointPair c, int mask_value);
	bool FloodReplace(int p, int orig, int updated);
};
