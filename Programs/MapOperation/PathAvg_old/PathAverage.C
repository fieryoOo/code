#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <errno.h>
#include "PathAverage.h"
using namespace std;

double Path::PathAverage(double lamda) {
   double Nmin = 3.; //(2 ~ 20?) Don't know much about sw kernel
   // define ellipse
   double f = dist/2.; // known values
   double bsqrmax = pow(f+lamda/2./Nmin, 2)-f*f; //maximum affective b from Nmin
   //double asqr = pow(f+lamda/2./N, 2), bsqr = asqr-f*f;

   char buff[300];
   float wmax = sqrt(bsqrmax); // project using wmax
   char ftmpname[20];
   sprintf(ftmpname, "./projecttmp.txt");
   sprintf(buff, "project %s -Dg -C%lf/%lf -E%lf/%lf -Lw -Q -S -W-%f/%f -Fpqz > %s", fname, P1.Lon(), P1.Lat(), P2.Lon(), P2.Lat(), wmax, wmax, ftmpname);
   system(buff);

   // read in project data
   //double N = 12., asqr = pow(f+lamda/2./N, 2), bsqr = asqr-f*f;
   double Nhaf = 12., Ncur, weight, weit=0., zsum=0.;
   FILE *fproj = NULL;
   if( (fproj = fopen(ftmpname, "r")) == NULL ) {
      perror("open file");
      return -1;
   }
   float a, xdat, ydat, zdat;
   while( fgets(buff, 300, fproj)!=NULL ) {
      sscanf(buff, "%f %f %f", &xdat, &ydat, &zdat);
      if( zdat != zdat ) continue;
      xdat -= f;
      // compute N from xdat, ydat, f, lamda
      a = 0.5 * (sqrt(pow((xdat-f),2)+ydat*ydat) + sqrt(pow((xdat+f),2)+ydat*ydat));
      Ncur = lamda/(2*(a-f)); // the closer to the great circle path the larger
      Ncur = Nhaf/Ncur;
      weight = exp(-0.5*Ncur*Ncur);
      weit += weight;
      zsum += (zdat * weight);
//cerr<<xdat<<" "<<ydat<<" "<<zdat<<"   "<<weight<<" "<<zsum<<" "<<Ncur<<" "<<Nhaf<<endl;
   }
   fclose(fproj);
   remove(ftmpname);
   if(weit==0.) zsum = -12345.;
   else zsum /= weit;
   return zsum;
 
}


int main(int argc, char *argv[])
{
   if (argc!=7) {
      cout<<"Usage: "<<argv[0]<<" [Input_Map (lon lat value)] [lon1] [lat1] [lon2] [lat2] [wavelength (km)]"<<endl;
      return -1;
   }

   /* get input parameters */
   char fname[100]; sprintf(fname, "%s", argv[1]);
   double lamda = atof(argv[6]);
   double lon1 = atof(argv[2]), lon2 = atof(argv[4]);
   if( lon1 < 0. ) lon1 += 360.;
   if( lon2 < 0. ) lon2 += 360.;
   Point<double> P1(lon1, atof(argv[3])), P2(lon2, atof(argv[5]));
   //cerr<<P1<<" "<<P2<<endl;  

   /* read in the map */
   Path pathcur(fname, P1, P2);
   cout<<pathcur.PathAverage(lamda)<<endl;

   return 0;
}
