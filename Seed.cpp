#include "Seed.h"
#define WDTHRESH 5000
#define DISTTHRESH 500
#define FORCECONST 250
#define FRICTION 0.5

Seed::Seed(unsigned char* graddat, VectorPair initial, int w, int h, bool stat)
{
	gradientdata = graddat;
	width = w;
	height = h;
	location.x = initial.x;
	location.y = initial.y;
	velocity.x = 0.0;
	velocity.y = 0.0;
	stationary = stat;
	ResetForce();
}

Seed::~Seed()
{
}

VectorPair Seed::GetPairForce(Seed* otherseed)
{
	VectorPair ret;
	int wd;
	float dist;
	VectorPair other;
	other = otherseed->GetLocation();
	dist = sqrt((location.x - other.x) * (location.x - other.x) + (location.y - other.y) * (location.y - other.y));
	if (dist > DISTTHRESH)
	{
		ret.x = 0;
		ret.y = 0;
		return ret;
	}
	if (dist < 1.0)
	{
		ret.x = rand() * 10.0;
		ret.y = rand() * 10.0;
		return ret;
	}
	wd = WDist(gradientdata, width, height, location.x, location.y, other.x, other.y);
	if (wd > WDTHRESH)
	{
		ret.x = 0;
		ret.y = 0;
		return ret;
	}
	ret.x = FORCECONST * (location.x - other.x) / (dist * wd);
	ret.y = FORCECONST * (location.y - other.y) / (dist * wd);
	if (false == stationary)
	{
		AccumulateForce(ret);
	}
	ret.x = -ret.x;
	ret.y = -ret.y;
	otherseed->AccumulateForce(ret);
	return ret;
}

void Seed::ResetForce()
{
	force.x = 0.0;
	force.y = 0.0;
}

void Seed::AccumulateForce(VectorPair f)
{
	force.x += f.x;
	force.y += f.y;
}

VectorPair Seed::GetLocation()
{
	return location;
}

void Seed::Move()
{
	if (false == stationary)
	{
		velocity.x = FRICTION * velocity.x + force.x;
		velocity.y = FRICTION * velocity.y + force.y;
		if ((abs(velocity.x) < 0.1) && (abs(velocity.y) < 0.1))
		{
			velocity.x = 0;
			velocity.y = 0;
		}
		location.x += velocity.x;
		if (location.x < 0)
		{
			location.x = 0;
		}
		else if (location.x >= width)
		{
			location.x = width - 1;
		}
		location.y += velocity.y;
		if (location.y < 0)
		{
			location.y = 0;
		}
		else if (location.y >= height)
		{
			location.y = height - 1;
		}
	}
}

bool Seed::GetStationary()
{
	return stationary;
}
