#define MAIN
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
using namespace std;

main (int na, char *argv[])
{
   if(na!=2)
     {
      cout<<"usage: extract_prediction [input_file]"<<endl;
      return 0;
     }
   FILE *f1,*f2;
   int i,j;
   int iev,ista,pn,pn2;
   char evn[10],evn1[10],stan[10],out_fname[20],str[30];
   float evlong,evlati,stalong,stalati,per[100],vel[100];

   f1=fopen(argv[1],"r");
   i=0;
   for(;;)
     {
      if(fscanf(f1,"%d %d %d %s %s %f %f %f %f",&iev,&ista,&pn,&evn[0],&stan[0],&evlati,&evlong,&stalati,&stalong)!=9) break;
      if(strcmp(evn,evn1)!=0){
	 printf("Extracting event: %s\n",evn);
         sprintf(str,"mkdir -p %s",evn);
         system(str);
        }
      strcpy(evn1,evn);
      for(i=0;i<pn;i++)
         if(fscanf(f1,"%f %f",&per[i],&vel[i])!=2) { printf("Error: Wrong period_num!"); return 0;}
      sprintf(out_fname,"%s/%s_%s.dat",evn,evn,stan);
      f2=fopen(out_fname,"w");
//      fprintf(f2,"%d  %d  %d  %s  %s  %.5f %.5f  %.5f %.5f\n",iev,ista,pn,evn,stan,evlati,evlong,stalati,stalong);
      for(i=0;i<pn;i++) fprintf(f2,"%9.5f  %8.5f\n",per[i],vel[i]);
      fclose(f2);
     }
}

