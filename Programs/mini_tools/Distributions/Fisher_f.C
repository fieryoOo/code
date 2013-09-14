#include <iostream>
#include "boost/random.hpp"
#include "boost/generator_iterator.hpp"
//#include "boost/math/distributions/fisher_f.hpp"
//#include "boost/random/fisher_f_distribution.hpp"
#include "fisher_f_distribution.hpp"
#include <time.h>
#include <sys/time.h>
using namespace std;

uint64_t ClockGetTime()
{
    timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    return (uint64_t)ts.tv_sec * 1000000LL + (uint64_t)ts.tv_nsec / 1000LL;
}

int main(int argc, char *argv[]) {
   if(argc!=2){
      cout<<"usage: "<<argv[0]<<" [array size]"<<endl;
      return 0;
   }

   typedef boost::mt19937 RNGType;
   RNGType gener( ClockGetTime() );
   boost::fisher_f_distribution<> fisher(1.,100.);
   boost::variate_generator< RNGType, boost::fisher_f_distribution<> > F(gener, fisher);

   cout<<F();

}
