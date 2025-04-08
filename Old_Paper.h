#pragma once
#define _USE_MATH_DEFINES
#include <vector>
#include <set>
#include "General.h"
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
	bool Add_g(int x, int y, float concentration);
	bool Add_d(int x, int y, float concentration);
	float GetR(int channel, float thickness); // Return R value (Kulbelka-Munk method) for a given thickness value and color channel.
	float GetT(int channel, float thickness); // Return T value (Kulbelka-Munk method) for a given thickness value and color channel.
	bool UpdateColor(Color* c, float r);
};

class Paper  // This is the object representing the paper on which painting will be done.  It includes the variables to calculate pigment movement, and also the details of each pigment.
{
private:
	int w, h;
	Color color;
	bool** M; // Mask for the wet-area of the paper.  Index order is [x][y].
	int** M_column; // Tracks values in each column of each chunk.  Addressed as M_column[chunk_number][column number in chunk].  0 = all false, 1 = some true, 2 = all true.
	std::set<int> M_chunks;  // The list of chunks (based on pigments[0]) for which there is at least one true value for M.
	bool TestMask(int chunk_index); // Return true if any element of M[x][y] is true within the identified chunk.  Also updates M_column.
	float** thickness;  // Thickness of paper at point x,y.  Values are between 0 and 1.0.  Index order is [x][y].
	float** u;  // Horizontal velocity on a staggered grid.  Index order is [x][y].  Size is [w+1][h].
	float** v;  // Vertical velocity on a staggered grid.  Index order is [x][y].  Size is [w][h+1].
	float** p;  // Pressure at point x,y.  Index order is [x][y].
	float** s;  // Water saturation.  Values are between 0 and 1.0.  Index order is [x][y].
	std::vector<Pigment*> pigments;  // Set of pointers to Pigments.
	int x_chunks, y_chunks; // The number of chunks in the x and y direction (each of size chunk_size).
	// Temporary working variables
	SparseFloatMatrix* g = NULL; // Pigment concentration used for MovePigment function (pre-move).
	SparseFloatMatrix* g_prime = NULL; // Pigment concentration used for MovePigment function (post-move).
	float** s_prime = NULL; // Saturation used for CapillaryFlow function.
	float** u_prime = NULL; // Horizonal velocity used in UpdateVelocities.
	float** v_prime = NULL; // Vertical velocity used in UpdateVelocities.
	float** M_prime = NULL; // Gaussian blurred version of mask M, used in EnforceBoundaries.  Values saved are actually 1 - M', to save on later calculations.
	float** M_prime_kernel = NULL; // Kernel for M_prime.
	float M_prime_kernel_total = 1.0f; // Kernal total value for M_prime.
	std::set<int> dried; // Set of pigment layers that have been dried.
	float*** u_delta;  // These store the chunk_size+1 by chunk_size+1 FloatArray's for each chunk that will be used for temporary calculation of velocity deltas.
	float*** v_delta;  // These store the chunk_size+1 by chunk_size+1 FloatArray's for each chunk that will be used for temporary calculation of velocity deltas.

	bool MoveWater(bool slope_velocity, bool debug_output);  // Move water in the shallow water layer.
	bool MovePigment(); // Move pigments in the shallow water layer.
	bool TransferPigment(); // Calculate pigment adsorption and desorption.
	bool CapillaryFlow(); // Simulate backruns by similating the capillary flows.
	bool UpdateVelocities(bool slope_velocity);  // Update water velocities.
	bool UpdateVelocitiesThreads(bool slope_velocity); // Update water velocities.  This version uses multi-threading.
	bool UpVelThreadworker(int chunk_id, float step_size); // Worker function for the multi-threading version of the UpdateVelocities function.
	bool RelaxDivergence(); // Relax the divergence of the velocity field.
	bool RelaxDivergenceThreads(); // Relax the divergence of the velocity field.  This version uses multi-threading.
	bool RelDivThreadWorker(int chunk_id); // Worker function for the multi-threading version of the RelaxDivergence function.  Returns the max delta value for the chunk.
	bool FlowOutward(); // Removes some water near the boundary of the wet-area mask (M), to induce flow outward.
	bool EnforceBoundaries(); // Set the velocities at each point beyond the wet-area mask (M) to zero.
	bool PaperSlope(); // Update u,v values based on the local slope of the paper.

	float MaxVelocity(); // Return the absolute value of the maximum velocity from the u and v vector fields.
	float capacity(int i, int j); // Return the calculation of the capacity c based on the thickness at i,j, and the maximum and minimum values of capacity.
	float u_vel(int i, int j); // Return the interpolated value of u at the center point of i,j (not staggered).
	float v_vel(int i, int j); // Return the interpolated value of v at the center of point i,j (not staggered).
	float uv_corner(int i, int j, int x, int y); // Return the value of u*v at the corner defined by x,y where x=0 means i-.5, x=1 means i+.5, y=0 means j-.5, y=1 means j+.5.

public:
	Paper(int width, int height, Color c, float saturation = 0.0f);
	~Paper();
	int GetWidth();
	int GetHeight();
	float** GetThickness();
	bool Tick(bool slope_velocity, bool debug_output=false); // Move main loop one clock tick forward in time.
	bool UpdateM_chunks(); // Update M_chunks and the values of M_column.
	bool Process(bool out=false); // Run the full set of ticks.
	std::vector<Pigment*> GetPigments(); // Return the set of Pigment objects for this Paper.
	int SetPigment(Pigment* pgmnt); // Sets a new pigment, and returns the zero-based index for it in the pigments vector.
	bool Dab(int i, int j, int radius, float saturation, float concentration, int pgmnt_num);
	bool Calc_M_prime(int x = -1, int y = -1); // Function to calculate M_prime.  If x and y are given, then the calcuation only updates the area around x,y.  Otherwise, full recalculation.
	bool OutputWaterConcentration(int layer);
	bool OutputDeposition(int layer);
	bool OutputPressure(int frame=0);
	bool OutputSaturation(int frame = 0);
	bool Render(unsigned char* data);  // Render pigment layers onto a background defined by the color variable, under white light.
	bool Dry();
	bool CheckDry(int layer); // True if the indicated layer has already been dried.
	bool SetVelocity(int x, int y, float u_vel, float v_vel, bool sum=false); // Set velocity at x and y.  If last term is true, add values rather than replace.
	bool SetBackgroundColor(Color c);
};