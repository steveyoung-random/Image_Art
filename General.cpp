#include "General.h"

float arccoth(float x) {
	if (std::abs(x) <= 1) {
		throw std::domain_error("Attempted to call arccoth on a value with an absolute value less than 1.\n");
		return 0.0f;
	}
	return 0.5 * std::log((x + 1.0f) / (x - 1.0f));
}



