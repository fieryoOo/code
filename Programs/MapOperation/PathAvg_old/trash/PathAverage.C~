#include <stdio.h>
#include <math.h>
#include <iostream>
#include <errno.h>
using namespace std;

int calc_dist(double lati1, double long1, double lati2, double long2, double *dist);

int main(int argc, char *argv[])
{
   if (argc!=7) {
      cout<<"Usage: "<<argv[0]<<" [Input_Map (lon lat value)] [lon1] [lat1] [lon2] [lat2] [wavelength (km)]"<<endl;
      return -1;
   }

   double dist, Nmin = 3.; //(2 ~ 20?) Don't know much about sw kernel
   double lon1 = atof(argv[2]), lat1 = atof(argv[3]), lon2 = atof(argv[4]), lat2 = atof(argv[5]);
   calc_dist(lat1, lon1, lat2, lon2, &dist);
   // define ellipse
   double f = dist/2., lamda = atof(argv[6]); // known values
   double bsqrmax = pow(f+lamda/2./Nmin, 2)-f*f; //maximum affective b from Nmin
   //double asqr = pow(f+lamda/2./N, 2), bsqr = asqr-f*f;

   char buff[300];
   float wmax = sqrt(bsqrmax); // project using wmax
   char fname[100], ftmpname[20];
   sprintf(fname, "%s", argv[1]);
   sprintf(ftmpname, "./projecttmp.txt");
   sprintf(buff, "project %s -Dg -C%lf/%lf -E%lf/%lf -Lw -Q -S -W-%f/%f -Fpqz > %s", fname, lon1, lat1, lon2, lat2, wmax, wmax, ftmpname);
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
      Ncur = lamda/2./(a-f)/Nhaf;
      weight = exp(-0.5*Ncur*Ncur);
      weit += weight;
      zsum += (zdat * weight);
//cerr<<xdat<<" "<<ydat<<" "<<zdat<<"   "<<weight<<" "<<zsum<<endl;
   }
   fclose(fproj);
   remove(ftmpname);
   if(weit==0.) zsum = -12345.;
   else zsum /= weit;
   cout<<zsum<<endl;
   
   return 0;
}
