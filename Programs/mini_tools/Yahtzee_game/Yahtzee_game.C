#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <limits.h>
#include "boost/random.hpp"
#include "boost/generator_iterator.hpp"
#include <sys/time.h>
using namespace std;

uint64_t ClockGetTime()
{
    timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    return (uint64_t)ts.tv_sec * 1000000LL + (uint64_t)ts.tv_nsec / 1000LL;
}

void Roll(int *dices, int nd) {
   typedef boost::mt19937 RNGType;
   RNGType gener( ClockGetTime() );
   boost::uniform_int<> uni(1,6);
   boost::variate_generator< RNGType, boost::uniform_int<> > dice(gener, uni);
   for(int i=0;i<nd;i++) dices[i] = dice();
}

void ReRoll(int *dices, int nd) {
   typedef boost::mt19937 RNGType;
   RNGType gener( ClockGetTime() );
   boost::uniform_int<> uni(1,6);
   boost::variate_generator< RNGType, boost::uniform_int<> > dice(gener, uni);
   int i, counter[6]={0}, cmax=0, tar;
   for(i=0;i<nd;i++) counter[dices[i]-1]++;
   for(i=0;i<6;i++) if(cmax<counter[i]) { cmax = counter[i]; tar = i+1; }
   if(cmax==1) tar = 0;
   for(i=0;i<nd;i++) if(dices[i]!=tar) dices[i] = dice();
}

int Yahtzee(int nd, int nr) {
    int i, dices[nd];
    Roll(&dices[0], nd);
    for(i=1;i<nr;i++) ReRoll(&dices[0], nd);
    for(i=1;i<nd;i++) if(dices[i]!=dices[0]) return 0;
    return 1;
}

int main (int argc, char *arg[]) {
   int ndice, niter=100000;
   int i, nyahtzee;
   double P;
   FILE *fout = fopen("Dice.num_Probability.txt", "w");
   for(ndice=1; ndice<100; ndice++) {
      cout<<"dice# "<<ndice;
      for(nyahtzee=0,i=0;i<niter;i++)
         nyahtzee += Yahtzee(ndice, ndice-2);
      P = (double)nyahtzee/i;
      cout<<": "<<P<<endl;
      fprintf(fout, "%d %lf\n", ndice, P);
   }
  fclose(fout);
}
