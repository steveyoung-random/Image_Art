#pragma once
#include <iostream>

int accumulate(unsigned char* data, int w, int h, int x, int y, float c);

int ipart(float x);

float wround(float x);

float fpart(float x);

float rfpart(float x);

int WDist(unsigned char* data, int w, int h, int x0, int y0, int x1, int y1);
