#pragma once

// Copyright (c) 2023-2024 Steve Young
// Licensed under the MIT License

#include <iostream>
#include "stb_image_write.h"
#include "General.h"
#include "Brush.h"
#include <sstream>
#include <string>

class GradData;
class SPixelData;

void write_png_to_mem(void* context, void* data, int size);

class ImageData {
private:
	unsigned char* data;
	float* data_wide;
	int width, height, colorchannels;
	Brush* brush;
public:
	ImageData(unsigned char* data_in, int w, int h, int n, bool frac_values=false);
	~ImageData();
	Color GetPixel(int x, int y);
	bool CollapseWideData();
	bool CreateBrush(FloatPointPair start, Color c, Color sec, int r, Paint_Properties prop);
	bool PaintCurve(std::vector<Corner> curve, SPixelData* mask = NULL, int mask_value = 0);
	unsigned char* GetData();
	int GetWidth();
	int GetHeight();
	int GetColorChannels();
	bool SetPixel(int x, int y, Color c);
	bool Reset();
	bool SetBackground(Color c);
	bool write_file(std::string filename);

	GradData* gen_gray(int channel = 0, int nchannel = 0);
	GradData* gen_diff(ImageData* image2);

	Color CIELABconvert(unsigned char r, unsigned char g, unsigned char b);
	Color CIELABconvert(Color input);
	Color RGBconvert(float L, float a, float b);
	Color RGBconvert(Color input);
};

