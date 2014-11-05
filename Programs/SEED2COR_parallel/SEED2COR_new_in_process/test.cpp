#include "MyOMP.h"
#include <cmath>
#include <iostream>

class ITER {
public:
	ITER(int imaxin) 
		: i(0), imax(imaxin) {}
	void Next() { i++; }
	bool InRange() { return i<imax; }
	void Process() { 
		data = sqrt(i*imax*0.75); 
		std::cout<<"from thread "<<omp_get_thread_num()<<": "<<data<<std::endl;
	}
	int index() { return i; }
private:
	int i, imax;
	float data;
};


int main() {
	ITER iter(100);
#pragma omp parallel shared(iter)
{
	int i;
	bool inrange = true;
	while( inrange ) {
		#pragma omp critical 
		{
			inrange = iter.InRange();
			i = iter.index();
			iter.Next();
		}
	
		if( inrange ) {
			int ithread=omp_get_thread_num();
			#pragma omp critical
			std::cerr<<ithread<<" "<<i<<std::endl;
		}
	}
}
	return 0;
}
