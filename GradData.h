#pragma once

// Copyright (c) 2023-2025 Steve Young
// Licensed under the MIT License

#include <iostream>
#include "stb_image_write.h"
#include "stb_image.h"
#include <sstream>
#include <string>

class ImageData;

class GradData {
private:
	unsigned char* data = NULL;
	int width, height;
public:
	GradData(unsigned char* gradient, int w, int h);
	GradData(const GradData&);
	GradData(std::string filename);
	~GradData();
	unsigned char GetPixel(int x, int y);
	int GetWidth();
	int GetHeight();
	bool write_file(std::string filename);
	GradData* gen_dilate_erode(bool dilate, int mode, int struct_size);
	GradData* gen_edge(GradData* dilate, GradData* erode, int xdiv, int ydiv, float konst);
	ImageData* Gradient2Image(int mode);
	GradData* Preprocess_Gray(int num_steps, unsigned char steps, unsigned char modes, int structsize);
	GradData* Generate_Gradient(int mode = 0, int struct_size = 3, int xd = 100, int yd = 100, float konst = 0);
};

