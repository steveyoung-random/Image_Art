#pragma once

// Copyright (c) 2023-2026 Steve Young
// Licensed under the MIT License

#include "Utils_CUDA.cuh"
#include "..\SuperPixel.h"

__global__ void watershed_step_1(int level, int width, int height, int* pixel_data, unsigned char* gradient_data, bool* ResultsArray);
__global__ void watershed_step_2(int width, int height, int* pixel_data, bool* ResultsArray);
__global__ void calc_bbox(int width, int height, int* pixel_data, int* d_bbox_array, int sp_count);
__global__ void find_initial_cell_regions(int width, int height, int xdiv, int ydiv, int buffer, int* pixdata, unsigned char* edge_data);
__global__ void find_seeds_cell_grow(int width, int height, int xdiv, int ydiv, unsigned char* edge_data, int* temp_spdata, int* delta_spdata, bool* collision_data, int* stage_sizes);
__global__ void resolve_seeds_initial(int width, int height, int xdiv, int ydiv, int buffer, int* region_spdata, int* temp_spdata, int* d_seeds);

__device__ void find_seeds_cell_grow_main_old(int width, int height, int x_offset, int y_offset, int cell_width, int cell_height, unsigned char* edge_data, int* temp_spdata, int* delta_spdata, bool* collision_data, int* stage_sizes, int* mask, int mask_identifier);
__device__ void calc_stage_sizes_old(int width, int N_cell_elements, int cell_width, int x_offset, int y_offset, int* temp_spdata, int* stage_sizes);
__device__ void detect_collision_old(int width, int cell_width, int cell_height, int x_offset, int y_offset, unsigned char level, unsigned char* edge_data, int* temp_spdata, bool* collision_data);
__device__ void detect_touching_regions(int width, int cell_width, int cell_height, int x_offset, int y_offset, int* temp_spdata, bool* collision_data, int* mask, int mask_identifier);
__device__ void handle_collisions_old(int width, int cell_width, int cell_height, int x_offset, int y_offset, int* temp_spdata, bool* collision_data, int* stage_sizes);
__device__ int find_target_identifier_old(int width, int cell_width, int cell_height, int x_offset, int y_offset, int iter_pos, int* contig_identifier, int* contig_size, int* temp_spdata, bool* collision_data, int* stage_sizes);
__device__ void zero_out_cell_old(int width, int cell_width, int cell_height, int x_offset, int y_offset, int* matrix);
__device__ void growth_step_1(int width, int cell_width, int cell_height, int x_offset, int y_offset, unsigned char level, unsigned char* edge_data, int* temp_spdata, int* delta_spdata, int* mask, int mask_identifier);
__device__ void growth_step_2(int width, int cell_width, int cell_height, int x_offset, int y_offset, int* temp_spdata, int* delta_spdata, bool* continue_step_ptr, int* mask, int mask_identifier);
__device__ void detect_multiple_regions(int width, int cell_width, int cell_height, int x_offset, int y_offset, int* numbers_seen, int* temp_spdata, bool* all_done_ptr, int* mask, int mask_identifier);
__device__ void region_flood_fill(int width, int height, int x_offset, int y_offset, int region_width, int region_height, unsigned char min, int seed_identifier, int* pixdata, unsigned char* edge_data, int* mask, int mask_identifier);
__device__ void resolve_seeds_main_old(int width, int height, int x_region_offset, int y_region_offset, int dx, int dy, int n_x_cells, int buffer, int* region_spdata, int* temp_spdata, int* d_seeds);


bool c_watershed_grow(int x, int y, SuperPixel* list_head); // CUDA version of watershed growth algorithm.
bool c_find_bbox(int x, int y, SuperPixel* list_head); // CUDA function to populate bounding boxes in SuperPixels based on pixelmap.
std::vector <PointPair> c_find_seeds(int width, int height, GradData* edge, SPixelData* pixeldata, int xdiv, int ydiv, int buffer); // CUDA function to find seeds for an image.

// Split-specific (new)
__global__ void resolve_seeds_split(int width, int height, int* d_geometry, int buffer, int* region_spdata, int* temp_spdata, int* d_seeds, int max_nx, int max_ny, int* mask, int* d_identifiers);
__global__ void find_split_cell_regions(int width, int height, int* d_geometry, int buffer, int* mask, int* d_identifiers, int* local_pixdata, unsigned char* edge_data);
__global__ void split_cell_grow(int width, int height, int* d_geometry, int buffer, unsigned char* edge_data, int* temp_spdata, int* delta_spdata, bool* collision_data, int max_identifiers, int* mask, int* d_identifiers);
__global__ void find_split_cell_regions(int width, int height, int* d_geometry, int buffer, int* mask, int* d_identifiers, int* local_pixdata, unsigned char* edge_data);
__global__ void split_grow(int width, int height, int* d_geometry, unsigned char* image_data, int color_channels, int* mask, unsigned char* edge, int* temp_pixdata, int* d_seeds, int max_nx, int max_ny, int* d_identifiers);

__device__ void resolve_seeds_main(int width, int height, int x_region_offset, int y_region_offset, int dx, int dy, int nx, int ny, int buffer, int* region_spdata, int* temp_spdata, int* d_seeds, int* mask, int mask_identifier);
__device__ void zero_out_cell(int width, int cell_width, int cell_height, int x_offset, int y_offset, int* matrix, int* mask, int mask_identifier);
__device__ void detect_collision(int width, int cell_width, int cell_height, int x_offset, int y_offset, unsigned char level, unsigned char* edge_data, int* temp_spdata, bool* collision_data, int* mask, int mask_identifier);
__device__ int find_target_identifier(int width, int cell_width, int cell_height, int x_offset, int y_offset, int iter_pos, int* contig_identifier, int* contig_size, int* temp_spdata, bool* collision_data, int* stage_sizes, int* mask, int mask_identifier);
__device__ void handle_collisions(int width, int cell_width, int cell_height, int x_offset, int y_offset, int* temp_spdata, bool* collision_data, int* stage_sizes, int* mask, int mask_identifier);
__device__ void calc_stage_sizes(int width, int N_cell_elements, int cell_width, int x_offset, int y_offset, int* temp_spdata, int* stage_sizes, int max_identifiers, int* mask, int mask_identifier);
__device__ void find_seeds_cell_grow_main(int width, int height, int x_offset, int y_offset, int cell_width, int cell_height, unsigned char* edge_data, int* temp_spdata, int* delta_spdata, bool* collision_data, int* stage_sizes, int max_identifiers, int* mask, int mask_identifier);


std::vector<PointPair> c_split_seeds(SuperPixel* head, std::vector<int>split_identifiers, ImageData* image);