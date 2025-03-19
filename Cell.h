#pragma once

// Copyright (c) 2023-2025 Steve Young
// Licensed under the MIT License

#include "SPixelData.h"
#include "GradData.h"
#include "SuperPixel.h"


class Cell {
private:

	SPixelData* pixeldata;
	GradData* gradientdata;
	int numPixels;
	RectQuad box, window;
	int center_x, center_y, buffer;
	SuperPixel* head;
	SuperPixel* tail;
	PointPair seed;

public:

	Cell(GradData* graddat, SPixelData* pixdat, int cx, int cy, int w, int h, int buf);
	~Cell();
	bool FindSeed(SPixelData* mask=NULL, int mask_value=0, bool diagonals=false);
	PointPair GetSeed();

};
