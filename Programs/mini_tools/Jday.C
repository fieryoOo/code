#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <time.h>
using namespace std;


int main ( int argc, char *argv[] )
{
   if(argc!=4){
      cout<<"Usage: calc_jday [year] [month] [day]"<<endl;
      return 0;
   }

   int y=atoi(argv[1]), m=atoi(argv[2]), d=atoi(argv[3]);
   int jd = 0;
   int i;

   for ( i = 1; i < m; i++ ) {
      if ( (i==1) || (i==3) || (i==5) || (i==7) || (i==8) || (i==10) ) jd += 31;
      else if (i==2) {
         if ( (y%400 == 0) || (y%4 == 0 && y%100 != 0 ) ) jd += 29;
         else jd += 28;
      }
      else jd += 30;
   }
   jd += d;

   printf("%03d\n",jd);

   return 1;
}

