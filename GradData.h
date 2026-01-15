#pragma once

// Copyright (c) 2023-2025 Steve Young
// Licensed under the MIT License

#include <iostream>
#include "stb_image_write.h"
#include "stb_image.h"
#include <sstream>
#include <string>
#ifdef USE_CUDA
#include "Cuda_Image_Art\Image_CUDA.cuh"
#endif

class ImageData;

class GradData {
	// Object for storing one-channel image data.  Name refers to gradiant data, which is
	// stored in this format.  Effectively, though, this is used for monochromatic images.
private:
	unsigned char* data = NULL; // The image data, stored in row major form.
#ifdef USE_CUDA
	unsigned char* c_device_data = NULL;  // This is the CUDA on-device version of data.
#endif
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

