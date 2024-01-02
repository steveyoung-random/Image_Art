#pragma once
#include <stdlib.h>
#include <math.h>
#include "WDist.h"

struct VectorPair {
	float x, y;
};

class Seed {
public:
	unsigned char* gradientdata;
	int width, height;
	VectorPair location, velocity, force;
	bool stationary;

	Seed(unsigned char* graddat, VectorPair initial, int w, int h, bool stat);
	~Seed();

	VectorPair GetPairForce(Seed* otherseed);
	void ResetForce();
	void AccumulateForce(VectorPair f);
	VectorPair GetLocation();
	void Move();
	bool GetStationary();
};
