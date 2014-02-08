#include <cstdio>
#include <iostream>
#include <algorithm>
#include <chrono>
#include <random>
#include <cstring>

int main() {
   unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
   std::default_random_engine generator (seed);
   std::uniform_real_distribution<float> distribution(-1.,1.);
   auto dice = std::bind ( distribution, generator );
   for(int i=0; i<10; i++) std::cerr<<dice()<<std::endl;
   return 0;
}
