#pragma once
#define _USE_MATH_DEFINES
#include <vector>
#include <set>
#include "General.h"
#include "Cuda_Image_Art\Brush_CUDA.cuh"
#include <iostream>
#include <random>
#include <sstream>
#include <iomanip>


class Pigment // This is an instance of a particular pigment, and keeps track of its presence in a Paper object (including in the water layer on it).
{
private:
	int w, h;
	SparseFloatMatrix* g = NULL;  // Pigment concentration in the shallow water layer.  0 to 1.
	SparseFloatMatrix* d = NULL;  // Pigment concentration in the deposition layer.  0 to 1.
	float K[3], S[3]; // K and S values for Kulbelka-Munk method, for each color channel (R G B).
	float rho, omega, gamma;
	int x_chunks, y_chunks; // The number of chunks in the x and y direction (each of size chunk_size).
public:
	Pigment(int width, int height, float Kr_in, float Kg_in, float Kb_in, float Sr_in, float Sg_in, float Sb_in, float rho_in, float omega_in, float gamma_in);
	Pigment(int width, int height, Color* c, float r, float rho_in, float omega_in, float gamma_in);
	~Pigment();
	SparseFloatMatrix* Get_g();
	SparseFloatMatrix* Get_d();
	float Get_rho();
	float Get_omega();
	float Get_gamma();
	int GetXChunks();
	int GetYChunks();
	bool UpdateColor(Color* c, float r);
	float GetK(int channel);
	float GetS(int channel);
};

class Paper  // This is the object representing the paper on which painting will be done.  It includes the variables to calculate pigment movement, and also the details of each pigment.
{
private:
	int w, h;
	Color color;
	bool* M; // Mask for the wet-area of the paper.
	float* thickness;  // Thickness of paper at point x,y.  Values are between 0 and 1.0.  Padded with gauss_radius on left, right, top, and bottom.
	float* u;  // Horizontal velocity on a staggered grid.  Size is [w+1][h].
	float* v;  // Vertical velocity on a staggered grid.  Size is [w][h+1].
	float* p;  // Pressure at point x,y.
	float* s;  // Water saturation.  Values are between 0 and 1.0.
	bool* M_chunks = NULL; // Device-side vector of length x_chunks*y_chunks, to track which chunks have any true M values.
	bool* host_M_chunks = NULL; // Host-side vector corresponding to M_chunks.
	std::vector<Pigment*> pigments;  // Set of pointers to Pigments.
	int x_chunks, y_chunks; // The number of chunks in the x and y direction (each of size chunk_size).
	int gauss_radius = 50; // Size of gauss_radius is added as padding to thickness, so it needs to be consistent.

	// Temporary working variables
	float* g = NULL; // Pigment concentration used for MovePigment function (pre-move).
	float* g_prime = NULL; // Pigment concentration used for MovePigment function (post-move).
	float* results1 = NULL; // Vector of length threads_per_block for storing intermediate results.
	float* results2 = NULL; // Vector of length threads_per_block for storing intermediate results.
	int results_length = 0; // Length of results vector.
	float* s_prime = NULL; // Saturation used for CapillaryFlow function.
	float* u_prime = NULL; // Horizonal velocity used in UpdateVelocities.
	float* v_prime = NULL; // Vertical velocity used in UpdateVelocities.
	float* M_prime = NULL; // Gaussian blurred version of mask M, used in EnforceBoundaries.  Values saved are actually 1 - M', to save on later calculations.
	float* M_prime_kernel = NULL; // Kernel for M_prime.
	float M_prime_kernel_total = 1.0f; // Kernal total value for M_prime.
	std::set<int> dried; // Set of pigment layers that have been dried.
	float* delta_matrix = NULL; // Holds delta values calculated in RelaxDivergence function.

	bool MoveWater(bool slope_velocity, bool debug_output);  // Move water in the shallow water layer.
	bool MovePigment(); // Move pigments in the shallow water layer.
	bool TransferPigment(); // Calculate pigment adsorption and desorption.
	bool CapillaryFlow(); // Simulate backruns by similating the capillary flows.
	bool UpdateVelocities(bool slope_velocity);  // Update water velocities.
	bool RelaxDivergence(); // Relax the divergence of the velocity field.
	bool FlowOutward(); // Removes some water near the boundary of the wet-area mask (M), to induce flow outward.
	bool PaperSlope(); // Update u,v values based on the local slope of the paper.
	float MaxVelocity(); // Return the absolute value of the maximum velocity from the u and v vector fields.

public:
	Paper(int width, int height, Color c, float saturation = 0.0f);
	~Paper();
	int GetWidth();
	int GetHeight();
	bool* GetM();
	float* GetS();
	float* GetP();
	float* GetFullG();
	int GetPadding();
	bool GetThickness(float* out_matrix);
	bool Tick(bool slope_velocity, bool debug_output=false); // Move main loop one clock tick forward in time.
	bool UpdateM_chunks(); // Update M_chunks as well as host_M_chunks.
	bool Process(bool out=false); // Run the full set of ticks.
	std::vector<Pigment*> GetPigments(); // Return the set of Pigment objects for this Paper.
	int SetPigment(Pigment* pgmnt); // Sets a new pigment, and returns the zero-based index for it in the pigments vector.
	bool Dab(int i, int j, int radius, float saturation, float concentration, int pgmnt_num);
	bool PaintArea(int* SPdata, int w, int h, RectQuad window, int identifier, float saturation, float concentration, int pgmnt_num);
	bool Calc_M_prime(int x = -1, int y = -1); // Function to calculate M_prime.  If x and y are given, then the calcuation only updates the area around x,y.  Otherwise, full recalculation.
	bool OutputWaterConcentration(int layer);
	bool OutputDeposition(int layer);
	bool OutputPressure(int frame=0);
	bool OutputSaturation(int frame = 0);
	bool Render(unsigned char* data);  // Render pigment layers onto a background defined by the color variable, under white light.
	bool Dry();
	bool CheckDry(int layer); // True if the indicated layer has already been dried.
	bool SetVelocity(int x, int y, float u_vel, float v_vel, bool sum = false, bool wait = false); // Set velocity at x and y.  Set wait to true if this is the last velocity to set before other calculations.  
	bool SetBackgroundColor(Color c);
};