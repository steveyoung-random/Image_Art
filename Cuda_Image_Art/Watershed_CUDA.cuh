#pragma once

// Copyright (c) 2023-2026 Steve Young
// Licensed under the MIT License

#include "Utils_CUDA.cuh"
#include "..\SuperPixel.h"

__global__ void update_values(int width, int height, int* pixel_data, bool* ResultsArray);
__global__ void watershed_step(int level, int width, int height, int* pixel_data, unsigned char* gradient_data, bool* ResultsArray);
__global__ void calc_bbox(int width, int height, int* pixel_data, int* d_bbox_array, int sp_count);

bool c_watershed_grow(int x, int y, SuperPixel* list_head); // CUDA version of watershed growth algorithm.
bool c_find_bbox(int x, int y, SuperPixel* list_head); // CUDA function to populate bounding boxes in SuperPixels based on pixelmap.

/* Better strategy for finding bounding boxes (instead of slow testing of each row and column for each superpixel):
* Create a shared memory array for a single row (or column) at a time.  If this is more than can go in shared memory, then handle in strips.
* Set up a global memory array of the bounding boxes and sizes for each superpixel.
* Spin up threads which each map to a superpixel (use striding to address all superpixels).
* Each thread goes through all pixels in the row (or column), which are in shared memory.
*    - Use an offset mechanism for them to read from different memory locations at each step (avoid all reading from the same location at the same time).
* As the thread goes through the pixels, it is keeping track of whether its target superpixel shows up in the row, and uses this information to keep track of
* the information for the superpixel (bounding box and cumulative size).
* Synchronize between rows (or columns) and then read in the next row (or column).
* Only write out to global memory once at the end, from the local values (think about this -- it seems to require more local memory for the striding across
* superpixels).
* 
* The key to this is that each thread is mapped 1-to-1 to each superpixel, so there is no worry about another thread overwriting data for the superpixel that
* it is working on.
* 
* Is there an issue with multiple blocks trying to read the same row (or column) into memory at the same time?  Offsetting these seems a little more complicated,
* but it may be doable.

*/