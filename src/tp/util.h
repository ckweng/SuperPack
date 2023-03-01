#ifndef SP_UTILS_H__
#define SP_UTILS_H__

#include <chrono>

using std::string;
using std::chrono::time_point;
using std::chrono::high_resolution_clock;

inline time_point<high_resolution_clock> clock_start() {
	return high_resolution_clock::now();
}

inline double time_from(const time_point<high_resolution_clock>& s) {
	return std::chrono::duration_cast<std::chrono::microseconds>(high_resolution_clock::now() - s).count();
}

#endif
