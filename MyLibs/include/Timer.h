#ifndef TIMER_H
#define TIMER_H

#include <ctime>
#include <unordered_map>
#include <iostream>
#include <chrono>

class Timer {
	using WClock = std::conditional< std::chrono::high_resolution_clock::is_steady,
												std::chrono::high_resolution_clock,
												std::chrono::steady_clock >::type;	// the Wall clock
public:
	Timer() { UpdateT0(); }

	void UpdateT0() { t0CPU = std::clock(); t0Wall = WClock::now(); }

	double CPUSec(const std::string &name = "t0") {
		return (std::clock() - t0CPU) / (double)CLOCKS_PER_SEC;
	}

	double WallSec(const std::string &name = "t0") {
		return std::chrono::duration_cast<std::chrono::duration<double>>(WClock::now() - t0Wall).count();
	}

private:
	std::clock_t t0CPU;
	WClock::time_point t0Wall;
};


class Timer2 {
	using WClock = std::conditional< std::chrono::high_resolution_clock::is_steady,
												std::chrono::high_resolution_clock,
												std::chrono::steady_clock >::type;	// the Wall clock
public:
	Timer2(const std::string &name_p0 = "t0") { 
		AddTimePoint(name_p0);
	}

	void AddTimePoint(const std::string &name) {
		timePoints.emplace(name, std::make_pair(std::clock(),WClock::now()));	
	}

	double CPUSec(const std::string &name = "t0") {
		return (std::clock() - timePoints[name].first) / (double)CLOCKS_PER_SEC;
	}

	double WallSec(const std::string &name = "t0") {
		return std::chrono::duration_cast<std::chrono::duration<double>>(WClock::now() - timePoints[name].second).count();
	}

private:
	std::unordered_map<std::string, std::pair<std::clock_t, WClock::time_point>> timePoints;
};

#endif
