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

int **Initial(int n, int *n1) {
   typedef boost::mt19937 RNGType;
   RNGType gener( ClockGetTime() );
   boost::uniform_int<> uni(0,1);
   boost::variate_generator< RNGType, boost::uniform_int<> > coin(gener, uni);

   int i, j;
   int **S = (int **) malloc (n * sizeof(int *));
   for(i=0;i<n;i++) S[i] = (int *) malloc (n * sizeof(int));
   for(*n1=0,i=0;i<n;i++)
      for(j=0;j<n;j++) {
         S[i][j] = 2*coin()-1;
         if(S[i][j]==1) *n1 = *n1 + 1;
      }
   return S;
}

int mod (int a, int b) {
   int c = a%b;
   if(c<0) c+=b;
   return c;
}

int ComputeE(int **S, int n) {
   if(n==2) return -S[0][0]*S[0][1]-S[0][1]*S[1][1]-S[1][1]*S[1][0]-S[1][0]*S[0][0];
   int i, j, E = 0;
   for(i=0;i<n;i++)
      for(j=0;j<n;j++) {
         E -= S[i][j]*S[mod(i-1,n)][j] + S[i][j]*S[i][mod(j-1,n)];
      }
   return E;
}

int Neighbor(int **S, int *in, int *jn, int E, int n) {
   typedef boost::mt19937 RNGType;
   RNGType gener( ClockGetTime() );
   boost::uniform_int<> uni(0,n-1);
   boost::variate_generator< RNGType, boost::uniform_int<> > rand(gener, uni);
   int i=rand(), j=rand();
   if(n==2) E += 2*(S[i][j]*S[mod(i-1,n)][j] + S[i][j]*S[i][mod(j-1,n)]);
   else E += 2*(S[i][j]*S[mod(i-1,n)][j] + S[i][j]*S[mod(i+1,n)][j] + S[i][j]*S[i][mod(j-1,n)] + S[i][j]*S[i][mod(j+1,n)]);
   *in = i; *jn = j;
   return E;
}

int accept(int E, int Enew, long double T) {
   typedef boost::mt19937 RNGType;
   RNGType gener( ClockGetTime() );
   boost::uniform_real<> uni(0,1);
   boost::variate_generator< RNGType, boost::uniform_real<> > rand(gener, uni);
   if(Enew<E) return 1;
   if(rand()<exp((E-Enew)/T)) return 1;
   return 0;
}

void Print(int **S, int n, int E, long double T, float perc) {
   int i, j;
   char name[100];
cout<<T<<" "<<perc<<" "<<E<<endl;
   sprintf(name, "snapshot_T%.2g_p%.2f_E%d.txt", double(T), perc, E);
cout<<name<<endl;
   FILE *fout = fopen(name, "w");
   for(i=0;i<n;i++) for(j=0;j<n;j++)
      fprintf(fout, "%2d %2d %d\n", j, i, S[i][j] );
   fclose(fout);
}

int SimAnnealing(int **S, float alpha, long double T, int n, int n1) {
   int i, in, jn, E, Enew;
   int flag = 1, nrej = 0;
   float perc;
 //for(i=0;i<n;i++) {for(int j=0;j<n;j++) cout<<S[i][j]<<" "; cout<<endl;}
   E = ComputeE(S, n);
cout<<"initial: "<<T<<" "<<n1<<" "<<E<<endl;
   Print(S, n, E, T, float(n1)/n/n*100.);
   while (nrej < n*n) {
      Enew = Neighbor(S, &in, &jn, E, n);
      if(accept(E, Enew, T)) {
         S[in][jn] *= -1;
         n1 += S[in][jn];
         E = Enew;
         nrej = 0;
         perc = float(n1)/n/n;
         if(flag && (perc>=0.8||perc<=0.2)) {
            Print(S, n, E, T, perc*100.);
            flag = 0;
         }
      }
      else nrej++;
      T *= alpha;
   }
   Print(S, n, E, T, float(n1)/n/n*100.);
   return E;
}

int main (int argc, char *arg[]) {
   if(argc!=4) {
      cout<<"usage: "<<arg[0]<<" [lattice size] [alpha] [Temperature]"<<endl;
      return -1;
   }
   int n1, n = atoi(arg[1]);
   if(n==1) {
      cout<<"lattice size has to be >1"<<endl;
      return -1;
   }
   int **S;
   float alpha = atof(arg[2]);
   long double T = atof(arg[3]);
   S = Initial(n, &n1);
   int E = SimAnnealing(S, alpha, T, n, n1);
   cout<<E<<endl;
   return 1;
}
