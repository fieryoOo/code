#ifndef THREADS_H
#define THREADS_H

#include <deque>
#include <thread>
#include <future>
#include <mutex>
#include <condition_variable>

class Threads {
public:
	const int nthreads;
	Threads(const int nthreads = std::thread::hardware_concurrency()) : nthreads(nthreads) {}

	~Threads() { join();	}

	template< typename Func, class ...Args >
	void Call( Func &&func, Args &&...args ) {
		{ // wait for any thread to finish
			std::unique_lock<std::mutex> Clocker(mCount);
			if( nworking == nthreads ) 
				Cthreaddone.wait(Clocker, [&]{ return nworking<nthreads; });
			nworking++;
		}
		//std::thread(&Threads::CountedCall<Func>, this, std::forward<Func>(func), wrap<Args>(args)...).detach();	//doesn't work for mysterious reasons
		std::thread([=, &func](){CountedCall(func, args...);}).detach();
	}

	void join() {
		// wait for all threads to finish
		std::unique_lock<std::mutex> Clocker(mCount);
		Cthreaddone.wait(Clocker, [&]{ return nworking==0; });
	}


private:
	int nworking = 0;
	std::mutex mCount;
	std::condition_variable Cthreaddone;

	template< class Func, class ...Args >
	void CountedCall( Func &&func, Args &&...args ) {
		func(args...);
		{ // notify the Caller when the function call is done
			std::unique_lock<std::mutex> Clocker(mCount);
			nworking--; 
		}
		Cthreaddone.notify_all();
	}
};


class Threads0 {
public:
	const int nthreads;
	Threads0(const int nthreads = std::thread::hardware_concurrency()) 
		: nthreads(nthreads), _futureQ(nthreads) {}

	~Threads0() { join(); }

	template< class... Args >
	void Call( Args&&... args ) {
		bool isWaiting = true;
		while( isWaiting ) {
			for(auto &future : _futureQ) {
				if( future.valid() && future.wait_for(std::chrono::milliseconds(0)) != std::future_status::ready ) continue;
				future = std::async(std::launch::async, std::forward<Args>(args)...); 
				isWaiting = false; break;
			}
		}
	}

	void join() {
		for(auto &future : _futureQ)
			if( future.valid() ) future.wait();
	}

private:
	std::deque<std::future<void>> _futureQ;
};

#endif
