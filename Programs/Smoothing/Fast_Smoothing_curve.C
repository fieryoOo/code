#define MAIN
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
using namespace std;

#define NMOD 5000

int Smoothing (char *infile, double hwidth) {
  FILE *ff;
  int i, j, jstart, idata, ism, nmod=0;
  char buff[300];
  double *datax = NULL, *datay = NULL, xtmp, ytmp, midx;
  double weight, weit;
  double dx=hwidth/2.,  alpha=0.5/hwidth/hwidth;

  if((ff=fopen(infile,"r"))==NULL) {
     cout<<"Can't open file "<<infile<<endl;
     return 0;
  }

  for(idata=0;;idata++) {
     if( nmod*NMOD <= idata ) {
	datax = (double *) realloc (datax, (++nmod)*NMOD * sizeof(double));
	datay = (double *) realloc (datay, nmod*NMOD * sizeof(double));
     }
     if((fgets(buff, 300, ff))==NULL) break;
     sscanf(buff,"%lf %lf", &xtmp, &ytmp);
     for(i=idata-1;i>=0&&xtmp<datax[i];i--) { 
        datax[i+1] = datax[i];
        datay[i+1] = datay[i];
     }
     datax[i+1] = xtmp; datay[i+1] = ytmp;
  }
  fclose(ff);

  i = (int)ceil((datax[0]+hwidth)/dx);
  sprintf(buff,"%s_sm\0",infile);
  ff=fopen(buff,"w");
  for(j=0, jstart=0;;i++) {
     midx = i*dx;
     if(midx>datax[idata-1]-hwidth) break;
     ytmp = 0;
     for(j=jstart;datax[j]<midx-hwidth;j++);
     jstart = j;
     for(weit=0.;datax[j]<=midx+hwidth&&j<idata;j++) {
        weight=exp(-alpha*pow((datax[j]-midx),2));
        ytmp += datay[j]*weight;
        weit += weight;
     }
     if(weit==0) continue;
     ytmp /= weit;
     fprintf(ff, "%lf %lf\n", midx, ytmp);
  }
  fclose(ff);
  free(datax); free(datay);

  return 1;
}

int main (int argc, char *argv[])
{ 
  if (argc != 3) {
  cout<<argv[0]<<" [input_file] [half_width]"<<endl;
  return 0;
  }

  Smoothing(argv[1], atof(argv[2]));

  return 1;
}
