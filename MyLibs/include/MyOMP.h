#ifndef MYOMP_H
#define MYOMP_H

#ifdef _OPENMP
#include <omp.h>
#else
typedef int omp_int_t;
inline omp_int_t omp_get_num_threads() { return 1; }
inline omp_int_t omp_get_max_threads() { return 1; }
inline omp_int_t omp_get_thread_num() { return 0; }
inline void omp_set_nested(bool) {}
typedef char omp_lock_t;
inline void omp_init_lock(omp_lock_t*) {}
inline void omp_destroy_lock(omp_lock_t*) {}
inline void omp_set_lock(omp_lock_t*) {}
inline void omp_unset_lock(omp_lock_t*) {}
#endif

#include <vector>
class RWLock {
public:
	// FCFS: any lock func (lockRead or lockWrite) that's called first locks first
	// when FCFS==false, lockWrite may appear to have a lower priorty (because it
	//	sets more locks and is hence slower) and be constantly blocked by lockRead
	RWLock(bool FCFS = true) 
	: nt(omp_get_max_threads()), lockV(nt), FCFS(FCFS) {
		if( FCFS ) omp_init_lock(&plock);
		for(auto &lock : lockV) omp_init_lock(&lock);
	}
	~RWLock() {
		if( FCFS ) omp_destroy_lock(&plock);
		for(auto &lock : lockV) omp_destroy_lock(&lock);
	}
	void lock4Read() {
		if( FCFS ) omp_set_lock(&plock);
		omp_set_lock(&( lockV[omp_get_thread_num()]) ); 
		if( FCFS ) omp_unset_lock(&plock);
	}
	void unlock4Read() { 
		omp_unset_lock(&( lockV[omp_get_thread_num()]) ); 
	}
	void lock4Write() { 
		if( FCFS ) omp_set_lock(&plock);
		for(auto &lock : lockV) omp_set_lock(&lock);
		if( FCFS ) omp_unset_lock(&plock);
	}
	void unlock4Write() { 
		for(auto &lock : lockV) omp_unset_lock(&lock);
	}
private:
	int nt;
	bool FCFS;
	omp_lock_t plock;
	std::vector<omp_lock_t> lockV;
};

#endif
