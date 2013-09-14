#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
using namespace std;

class Temp1 {
public:
   int i, ndat;
   int *xdat, *ydat;
   int Nint(float a) {
      return (int)floor(a+0.5);
   }
   char buff[300];
   bool Producer(char *fname) {
      if( !Read(fname, &xdat, &ydat, &ndat) ) return false;
      for(i=0;i<ndat;i++) {cout<<i<<": "; cout<<xdat[i]<<" "<<ydat[i]<<endl; }
      delete [] xdat;
      delete [] ydat;
      return true;
   }

private:
   FILE *fin;
   int Read(char *fname, int **xdat, int **ydat, int *ndat) {
      if((fin=fopen(fname, "r"))==NULL) {
         cout<<"Cannot access file "<<fname<<endl;
         return 0;
      }
      float xtmp, ytmp;
      *xdat = new int[200];
      *ydat = new int[200];
      for(i=0; i<100 && fgets(buff, 300, fin)!=NULL; i++) {
         sscanf(buff, "%f %f", &xtmp, &ytmp);
         (*xdat)[i] = Nint(xtmp);
         (*ydat)[i] = Nint(ytmp);
      }
      fclose(fin);
      *ndat = i;
      return 1;
   }
};

int main(int argc, char **argv)
{
   if(argc!=2) {
      cout<<"Usage: "<<argv[0]<<" [input file]"<<endl;
      exit(-1);
   }
   double (*a)[200];
   cout<<sizeof(*a)<<endl;;
   a = new double[2][200]();
   cout<<sizeof(*a)<<endl;;
   int i, j;
   for(i=0;i<2;i++) for(j=0;j<200;j++) a[i][j] = i;
//   for(i=0;i<2;i++) for(j=0;j<200;j++) cout<<a[i][j]<<endl;
   delete [] a;
   cout<<sizeof(*a)<<endl;
   cout<<"!!"<<endl;
   cout<<a[0][100]<<endl;
   cout<<a[1][100]<<endl;
   //for(i=0;i<2;i++) for(j=0;j<200;j++) cout<<a[i][j]<<endl;
exit(0);
   Temp1 tmp;
   cout<<sizeof(Temp1)<<" "<<sizeof(tmp)<<endl;
   if( !tmp.Producer(argv[1]) ) return 0;
   return 1;
}
