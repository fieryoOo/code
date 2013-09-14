#include <string.h>
#include <iostream>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
using namespace std;

main(int argc, char *argv[])
{
   if (argc!=2) {
      cout<<"usage: Calc_time [in_file (hour1 min1 sec1 hour2 min2 sec2)]"<<endl;
      return 0;
   }

   FILE *ff;
   char buff[300];
   double t0, t1, t2, tt=0;
   double h1, m1, s1, h2, m2, s2;
   int hour, min, sec;
   if((ff=fopen(argv[1], "r")) == NULL) {
      cout<<"Cannot open file "<<argv[1]<<endl;
      return 0;
   }
   for(;;) {
      if((fgets(buff, 300, ff)) == NULL) break;
      sscanf(buff, "%lf %lf %lf %lf %lf %lf", &h1, &m1, &s1, &h2, &m2, &s2);
      t1 = h1*60*60+m1*60+s1;
      t2 = h2*60*60+m2*60+s2;
      t0 = t2 - t1;
      if(t0<0) t0+=86400.;
      tt += t0;
   }
   fclose(ff);
   hour = (int)(tt/3600.);
   min = (int)((tt-hour*3600.)/60.);
   sec = (int)(tt-hour*3600.-min*60.+0.5);
   cout<<"Total working time: "<<hour<<"hours "<<min<<"minutes "<<sec<<"seconds."<<endl;
   if(hour>19) cout<<"Your warning is on its way!!!"<<endl;

  return 1;
}
