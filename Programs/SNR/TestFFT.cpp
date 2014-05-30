#include "SacRec.h"
#include <cstring>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <random>

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

/* --- main --- */
inline int nint( float x ) { return (int)floor(x+0.5); }
int main() {
   SacRec sacin("SACs/COR_M12A_P21A.SAC");
   sacin.Load();
   int npts = sacin.shd.npts;
   float delta = sacin.shd.delta;
   static constexpr float pi = 3.1415927;

   /* ----- 1. sin * sin ----- */
   // sin(2 pi t / per1) * sin( 2 pi t / per2) for half per2 period
   float per1 = 10., per2 = 100.;
   for(int i=0; i<npts; i++) sacin.sig[i] = 0.;
   for(int i=npts/2; i<npts/2+nint(per2*0.5/delta); i++) {
      float t = sacin.shd.b + i * delta;
      sacin.sig[i] = 10. * sin(2.*pi*t/per1) * sin(2.*pi*t/per2);
   }
   // FFT
   SacRec sacam;
   sacin.ToAm(sacam);
   sacin.Write("SinSin.SAC");
   sacam.Write("SinSin.am");

   /* ----- 2. sin ----- */
   for(int i=0; i<npts; i++) {
      float t = sacin.shd.b + i * delta;
      sacin.sig[i] = 10. * sin(2.*pi*t/per1);
   }
   // FFT
   sacin.ToAm(sacam);
   sacin.Write("Sin.SAC");
   sacam.Write("Sin.am");

   /* ----- 3. Gaussian noise ----- */
   Rand randO;
   for(int i=0; i<npts; i++) {
      float t = sacin.shd.b + i * delta;
      sacin.sig[i] = randO.Normal();
   }
   // FFT
   sacin.ToAm(sacam);
   // cut the noise shorter to see what happens
   SacRec sacam_2;
   sacin.cut(0, 50, sacam_2);
   sacam_2.ToAm(sacam_2);
   sacin.Write("Gaussian.SAC");
   sacam.Write("Gaussian.am");
   sacam_2.Write("Gaussian.am2");


   return 0;
}
