#pragma once

// Copyright (c) 2023-2025 Steve Young
// Licensed under the MIT License

#include <iostream>
#include "SPixelData.h"
#include "Cell.h"
#include "SuperPixel.h"
#include "GradData.h"
#include <fstream>
#include <vector>
#include <algorithm>

#include "stb_image_write.h"
#include "stb_image.h"

class ColorPalette {
private:
	int size;
	std::vector<Color> centroid_color;

public:
	ColorPalette(int num, SuperPixel* sp);
	~ColorPalette();
	Color GetColorByIndex(int index);
	Color TranslateColor(Color c);
	std::string TranslateColorToString(Color c);
};

struct MCColorBucket { // Color bucket used for Median Cut method.
	unsigned char max[3]; // Use indexes for the three color channels.
	unsigned char min[3];
	unsigned char centroid[3];
};

struct SuperPixelIndex
{
	int identifier;
	SuperPixel* sp;
	bool operator< (const SuperPixelIndex& rhs) const
	{
		return identifier < rhs.identifier;
	}
};

class WorkSpace {
private:
	int width, height, colorchannels;
	int xdiv = 100;
	int ydiv = 100;
	GradData* gray = NULL;
	GradData* temp_gray = NULL;
	GradData* processed_gray = NULL;
	GradData* dilate = NULL;
	GradData* erode = NULL;
	GradData* edge = NULL;
	unsigned char* data = NULL; // Image data loaded by stbi_load and released by stbi_image_free.
	ImageData* data_revised = NULL;
	GradData* data_diff = NULL;
	ImageData* image = NULL;
	SPixelData* pixeldata = NULL;
	Cell* cell = NULL;
	SuperPixel* list_head = NULL;
	SuperPixel* list_tail = NULL;
	SuperPixel* skeleton_head = NULL;
	SuperPixel* skeleton_tail = NULL;
	SPixelData* skeleton_pixeldata = NULL;
	SuperPixel* processed_head = NULL;
	SuperPixel* processed_tail = NULL;
	SPixelData* processed_pixeldata = NULL;
	SuperPixel* current = NULL;
	RectQuad bbox;
	bool diagonals;
	std::set<SuperPixelIndex> superpixel_set;
	ColorPalette* palette = NULL;

public:

	WorkSpace(std::string filename, int channel = 0, int nchannel = 0, bool diagonals = false);
	WorkSpace(std::string sp_filename, bool diagonals = false);
	WorkSpace(const WorkSpace& t);
	~WorkSpace();
	bool Reset();
	bool write_file(std::string filename);
	bool Generate_Gradient(int mode = 0, int struct_size = 3, int xd = 100, int yd = 100, float konst = 0);
	bool InitialSuperPixels(std::string seeds = "");
	bool Watershed();
	bool SetAveColors();
	bool CombineSuperPixels(float colormatch);
	bool Postprocess_SuperPixels(int num_steps, unsigned char steps, unsigned char modes, int structsize);
	bool WriteSeeds(std::string filename);
	bool WriteSuperPixelCSV(std::string filename);
	bool WriteSuperPixels(std::string filename);
	bool WriteSuperPixelsSVG(std::string filename, int mode = 0, bool polygon = false, bool fine = false, int palette_size=0, int contrast_radius=0);
	bool WritePaintCurvesSVG(std::string filename);
	bool SplitSuperPixels(float num_sigmas);
	bool ThinSuperPixels(bool only_one=false, bool glitch3 = false, SuperPixel* sp = NULL, SPixelData* sd = NULL);
	bool FindPaths(int mode = 0, bool polygon = true, bool fine = true);
	bool FindNeighbors(int mode = 0);
	bool InsertSPIntoIndex(SuperPixel* sp);
	SuperPixel* GetByIdentifier(int key);
	bool RecalcSuperPixelSet();
	bool CalculateSizes();
	bool FindPalette(int mode, int palette_size);
	bool ReduceToPalette(int mode, int palette_size);
	bool SeparateDiscontinuousSuperPixel(SuperPixel* sp, int point);  // If SuperPixel has disconnected parts, split each contiguous part into its own SuperPixel.
	bool SimplifySuperPixels(SuperPixel* sp);
	bool DeleteSkeleton();
	bool DeleteSuperPixelSet(SuperPixel* sp);

	ImageData* GenerateImage(int mode, Paint_Properties prop);
	ImageData* Gradient2Image(int mode = 0);

	bool GenerateContrasts(int mode = 0, int contrast_radius = 0);

	int GetXdiv();
	int GetYdiv();
	SuperPixel* GetHead();
	ImageData* GetImage();
	GradData* GetGray();
	SPixelData* GetPixeldata();

	bool SetXdiv(int xd);
	bool SetYdiv(int yd);

	// Pass through to ImageData class.
	GradData* gen_gray(ImageData* image, int channel = 0, int nchannel = 0);
	GradData* gen_diff(ImageData* image1, ImageData* image2);

	// Pass through to GradData class.
	GradData* gen_dilate_erode(GradData* gray, bool dilate, int mode, int struct_size);
	GradData* gen_edge(GradData* dilate, GradData* erode, int xdiv, int ydiv, float konst);
	bool Preprocess_Gray(int num_steps, unsigned char steps, unsigned char modes, int structsize);
};
