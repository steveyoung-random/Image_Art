#pragma once

// Copyright (c) 2023-2026 Steve Young
// Licensed under the MIT License

#include "Utils_CUDA.cuh"
#include "..\SuperPixel.h"

__global__ void update_values(int width, int height, int* pixel_data, bool* ResultsArray);
__global__ void watershed_step(int level, int width, int height, int* pixel_data, unsigned char* gradient_data, bool* ResultsArray);
__global__ void calc_bbox(int width, int height, int* pixel_data, int* d_bbox_array, int sp_count);
__global__ void find_initial_cell_regions(int width, int height, int xdiv, int ydiv, int buffer, int* pixdata, unsigned char* edge_data);
__global__ void cell_grow(int width, int height, int xdiv, int ydiv, unsigned char* edge_data, int* temp_spdata, int* delta_spdata, bool* collision_data, int* stage_sizes);
__global__ void resolve_seeds(int width, int height, int xdiv, int ydiv, int buffer, int* spdata, int* temp_spdata, int* d_seeds);

bool c_watershed_grow(int x, int y, SuperPixel* list_head); // CUDA version of watershed growth algorithm.
bool c_find_bbox(int x, int y, SuperPixel* list_head); // CUDA function to populate bounding boxes in SuperPixels based on pixelmap.
std::vector <PointPair> c_find_seeds(int x, int y, GradData* edge, SPixelData* pixeldata, SuperPixel* list_head, int xdiv, int ydiv, int buffer); // CUDA function to find seeds for 