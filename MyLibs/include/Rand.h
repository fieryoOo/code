#ifndef RAND_H
#define RAND_H

#include <random>
#include <chrono>
//#include <thread>

/* -------------------- the RNG class-------------------- */
class Rand {
   std::default_random_engine generator1;
   std::uniform_real_distribution<float> d_uniform;
   std::normal_distribution<float> d_normal;
public:
   Rand() /* add a true random number from std::random_device to the time seed to ensure thread safety */
      : generator1( std::chrono::system_clock::now().time_since_epoch().count() + std::random_device{}() )
      , d_uniform(0., 1.)
      , d_normal(0., 1.) {}

   //~Rand() {} // should not be defined!!!

   inline float Uniform() { return d_uniform(generator1); }
   inline float Normal() { return d_normal(generator1); }

};

#endif
