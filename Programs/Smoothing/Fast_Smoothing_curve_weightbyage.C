#define MAIN
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
using namespace std;

int main (int argc, char *argv[])
{
  if (argc != 4) {
  cout<<argv[0]<<" [input_file] [half_width] [age]"<<endl;
  return 0;
  }

  FILE *ff;
  int i, j, jstart, idata, ism;
  char buff[300];
  double datax[5000], datay[5000], age[5000], xtmp, ytmp, atmp, midx;
  double weight, weitl, weitr, agel, ager, dage, vl, vr;
  double hwidth=atof(argv[2]), dx=hwidth/2.,  alpha=0.5/hwidth/hwidth, cage = atof(argv[3]);

  if((ff=fopen(argv[1],"r"))==NULL) {
     cout<<"Can't open file "<<argv[1]<<endl;
     return 0;
  }

  for(idata=0;;idata++) {
     if((fgets(buff, 300, ff))==NULL) break;
     sscanf(buff,"%lf %lf %lf", &xtmp, &ytmp, &atmp);
     for(i=idata-1;i>=0,xtmp<datax[i];i--) { 
        datax[i+1] = datax[i];
        datay[i+1] = datay[i];
        age[i+1] = age[i];
     }
     datax[i+1] = xtmp; datay[i+1] = ytmp; age[i+1] = atmp;
  }
  fclose(ff);

  i = (int)ceil((datax[0]+hwidth)/dx);
  sprintf(buff,"%s_sm\0",argv[1]);
  ff=fopen(buff,"w");
  for(j=0, jstart=0;;i++) {
     midx = i*dx;
     if(midx>datax[idata-1]-hwidth) break;
     vl = 0.; vr = 0.;
     agel = 0.; ager = 0.; 
     for(j=jstart;datax[j]+hwidth<midx;j++);
     jstart = j;
     for(weitr=0., weitl=0.;datax[j]-hwidth<=midx;j++) {
        dage = age[j]-cage;
        weight=exp(-alpha*pow((datax[j]-midx),2))/(fabs(age[j]-cage)+0.3);
        if(dage>0) {
           vr += datay[j]*weight;
           ager += age[j]*weight;
           weitr += weight;
        }
        else {
           vl += datay[j]*weight;
           agel += age[j]*weight;
           weitl += weight;
        }
     }
     if(weitr==0||weitl==0) continue;
     vr /= weitr; ager /= weitr;
     vl /= weitl; agel /= weitl;
    //cout<<agel<<": "<<vl<<"  "<<ager<<": "<<vr<<"  cage: "<<cage<<endl;
     ytmp = vl + (vr-vl)*(cage-agel)/(ager-agel);
    //cout<<ytmp<<endl;
     fprintf(ff, "%lf %lf %lf\n", midx, ytmp, ager);
  }
  fclose(ff);

  return 1;
}
