#ifndef TIMER_H
#define TIMER_H

#include <ctime>
#include <unordered_map>
#include <iostream>

class Timer {
public:
	Timer(const std::string &name_p0 = "t0") { 
		AddTimePoint(name_p0);
	}

	void AddTimePoint(const std::string &name) {
		timePoints.emplace(name, clock());	
	}

	double SecondsElapsed(const std::string &name = "t0") {
		return (clock() - timePoints[name]) / (double)CLOCKS_PER_SEC;
	}

private:
	std::unordered_map<std::string, clock_t> timePoints;
};

#endif
