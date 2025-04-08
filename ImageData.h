#pragma once

// Copyright (c) 2023-2025 Steve Young
// Licensed under the MIT License

#include <iostream>
#include "stb_image_write.h"
#include "General.h"
#include "Brush.h"
#ifdef USE_CUDA
#include "Paper.h"
#endif
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
#ifdef USE_CUDA
	Paper* paper;
#endif
	Color background_color;
public:
#ifdef USE_CUDA
	ImageData(unsigned char* data_in, int w, int h, int n, bool frac_values, bool watercolor=false);
#else
	ImageData(unsigned char* data_in, int w, int h, int n, bool frac_values);
#endif
	~ImageData();
	Color GetPixel(int x, int y);
	bool CollapseWideData(bool dither = false);
	bool CreateBrush(FloatPointPair start, Color c, Color sec, int r, Paint_Properties prop, int pigment_index = -1);
	bool PaintCurve(std::vector<Corner> curve, SPixelData* mask, int mask_value, bool use_mask=false, SPixelData* extinguish_mask = NULL);
	unsigned char* GetData();
	int GetWidth();
	int GetHeight();
	int GetColorChannels();
	bool SetPixel(int x, int y, Color c);
	bool Reset();
	bool SetBackground(Color c);
	bool write_file(std::string filename);
#ifdef USE_CUDA
	bool ProcessWatercolor();
	bool RenderWatercolor();
	Paper* GetPaper();
#endif

	GradData* gen_gray(int channel = 0, int nchannel = 0);
	GradData* gen_diff(ImageData* image2);

	Color CIELABconvert(float r, float g, float b);
	Color CIELABconvert(Color input);
	Color RGBconvert(float L, float a, float b);
	Color RGBconvert(Color input);
};

