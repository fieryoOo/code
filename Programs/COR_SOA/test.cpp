#include <iostream>
#include <thread>
#include <cmath>
#include "MyOMP.h"
#include "Threads.h"
#include "Timer.h"

template<class Func, class ...Args>
void Caller( Func &&func, Args &&...args ) { func(args...); }

//int ncalls1 = 0, ncalls2 = 0;
int sum(0); std::mutex msum;
void DoSum(const int& val) {
//void DoSum(int &sum, const int val) {
	//for(int i=0; i<10000; i++) std::this_thread::yield();
	std::this_thread::sleep_for(std::chrono::milliseconds(1));
	auto adder = log(1+std::sqrt(abs(log(val))));
	adder = pow(adder-1, 100);
	msum.lock(); sum += val; msum.unlock();
//	#pragma omp critical
//	sum += adder;
}

int main() {
	int npts = 10, npts2 = 1000; float data[npts];
	for(int i=0; i<npts; i++) data[i] = i;

	Timer timer;

	Threads threads(12);
	for(int i=1; i<npts; i++) {
		for(int j=0; j<npts2; j++) threads.Call(DoSum, j+data[i]);  
	}
	threads.join();

/*
	#pragma omp parallel for schedule (auto)
	for(int i=1; i<npts; i++) {
		for(int j=0; j<npts2; j++) Caller(DoSum, j+data[i]);
	}
*/

//std::this_thread::sleep_for(std::chrono::seconds(1));
	std::cerr<<timer.CPUSec()<<" "<<timer.WallSec()<<" "<<sum<<" "<<std::endl;
	timer.UpdateT0();
	return 0;
}
