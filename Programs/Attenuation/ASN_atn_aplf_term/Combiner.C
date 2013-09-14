#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
using namespace std;

//main(int na, char *arg[])
int Combiner (char *csta, char *dir, float *perlst, int nper, double *slong, double *slati, char **stalst, double **trvt, int *nsta)
{
/*
  if (na!=4)
    {
      cout<<"usage: File_Combiner [in_file_lst] [column # (0 for entire row)] [out_file]"<<endl;
      return 0;
    }
*/
  FILE *f1, *f2, *fout;
  char filename[100], buff[100];
  int i, j, ista=0;
  double dtmp;

  *nsta = 0;
  for(i=0;i<nper;i++){
     sprintf(filename,"%s/Ph_Amp_Map_%.0fsec/%s_center_ph_amp_map_v2",dir,perlst[i],csta);
     if(!(f2=fopen(filename,"r"))){
        cout<<"Can't open file "<<filename<<" to read!"<<endl;
        exit(0);
     }
     for(;;) {
        if(fgets(buff,100,f2) == NULL)break;
        sscanf(buff,"%lf %lf %lf %lf %lf %s",&slong[ista], &slati[ista], &trvt[i][ista], &dtmp, &dtmp, &stalst[ista][0]);
        if(slong[ista]<0) slong[ista] += 360;
        for(j=0;j<*nsta;j++)if(strcmp(stalst[ista],stalst[j]) == 0) {
           trvt[i][j] = trvt[i][ista];
           ista--;
           break;
        }
        ista++;
     }
     fclose(f2);
     *nsta = ista;
  }
  for(i=0;i<*nsta;i++) {
     if(strcmp(stalst[i], csta)) continue;
     if(i==0) break;
     sprintf(buff,"%s",stalst[0]);
     sprintf(stalst[0],"%s",stalst[i]);
     sprintf(stalst[i],"%s",buff);
     dtmp=slong[0]; slong[0]=slong[i]; slong[i]=dtmp;
     dtmp=slati[0]; slati[0]=slati[i]; slati[i]=dtmp;
     for(j=0;j<nper;j++) {
        dtmp=trvt[j][0]; trvt[j][0]=trvt[j][i]; trvt[j][i]=dtmp;
     }
     break;
  }
  if(i==*nsta) {
     cout<<"Can't find center station "<<csta<<" in any of the ph_amp_map_v2 files"<<endl;
     exit(0);
  }

}
