#include <iostream>
#include "boost/random.hpp"
#include "boost/generator_iterator.hpp"
using namespace std;

#define NSTA 3000

int main (int argc, char *argv[])
{
   if(argc!=4){
      cout<<"usage: add_noise [infile] [outf_name] [Gaussian_half_width]"<<endl;
      return 0;
   }

   typedef boost::mt19937 RNGType;
   RNGType gener( time(0) );
   boost::normal_distribution<> normal(0,atof(argv[3]));
   boost::variate_generator< RNGType, boost::normal_distribution<> > Gauss(gener, normal);

   FILE *ff;
   char buff[300];
   int i, nsta;
   double lon[NSTA], lat[NSTA], dat[NSTA];

   if((ff=fopen(argv[1],"r"))==NULL) {
      cout<<"Cannot open file: "<<argv[1]<<endl;
      exit(0);
   }
   for(i=0;;) {
      if(fgets(buff, 300, ff)==NULL) break;
      if((sscanf(buff, "%lf %lf %lf", &lon[i], &lat[i], &dat[i]))!=3) {
         cout<<"Wrong format: "<<buff<<"  Skipped."<<endl;
         continue;
      }
      dat[i] += Gauss();
      i++;
   }
   fclose(ff);
   nsta=i;

   ff=fopen(argv[2],"w");
   for(i=0;i<nsta;i++) fprintf(ff, "%lf %lf %lf\n", lon[i], lat[i], dat[i]);
   fclose(ff);

   return 1;
}

