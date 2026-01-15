#pragma once

// Copyright (c) 2023-2025 Steve Young
// Licensed under the MIT License

#include "Utils_CUDA.cuh"

__global__ void process_dilate_erode_rect_1(int x, int y, unsigned char* input, unsigned char* output, bool isdilate, int struct_size);
__global__ void process_dilate_erode_rect_2(int x, int y, unsigned char* input, unsigned char* output, bool isdilate, int struct_size);
__global__ void process_dilate_erode_disc(int x, int y, unsigned char* input, unsigned char* output, bool isdilate, int struct_size);
__global__ void process_gen_grayscale(unsigned char* c_device_data_input, int N, unsigned char* c_device_data_output);
__global__ void process_gen_grayscale_c(unsigned char* c_device_data_input, int N, int c, unsigned char* c_device_data_output);
__global__ void process_gen_grayscale_c_nc(unsigned char* c_device_data_input, int N, int c, int nc, unsigned char* c_device_data_output);

unsigned char* c_gen_dilate_erode(int x, int y, unsigned char* input, bool isdilate, int mode, int struct_size);
bool c_gen_gray(unsigned char* c_device_data_input, int x, int y, int channels, int c, int n, unsigned char* c_device_data_output);
