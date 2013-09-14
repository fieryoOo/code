#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

#include <iostream>

//#include <sys/types.h>
#include <sys/stat.h>
//#include <fcntl.h>
#include "/home/tianye/code/Programs/head/mysac.h"

using namespace std;

void bin2hex( short* bin, char* hex )
{
   short i, j, hv;
   for(i=0;i<8;i++){
      hv = 0;
      for(j=0;j<4;j++){
         hv += bin[i*4+3-j]*(int)pow(2,j);
      }
      if( hv >= 0 && hv <=9 ) hex[i] = hv + '0';
      else hex[i] = hv - 10 + 'A';
   }
   hex[i]='\0';
}



SAC_HD *read_sac (char *fname, float *sig, SAC_HD *SHD, int nmax)
/*--------------------------------------------------------------------------*/
/* function to read sac files given the name, fname. The function outputs the time signal to the pointer sig
, fills the header SHD, if the signal has fewer than nmax points */
{
  FILE *fsac;

  if((fsac = fopen(fname, "rb")) == NULL) {
    fprintf(stderr,"read_sac: Could not open %s\n", fname);
    return NULL;
  }

  if ( !fsac ) {
    /*fprintf(stderr,"file %s not find\n", fname);*/
    return NULL;
  }

//  if ( !SHD ) SHD = &SAC_HEADER;

  fread(SHD,sizeof(SAC_HD),1,fsac);

  if ( SHD->npts > nmax ) {
    fprintf(stderr,"ATTENTION !!! %s npts is limited to %d.\n", fname, (int)nmax);
    SHD->npts = nmax;
  }

  fread(sig,sizeof(float),(int)(SHD->npts),fsac);
  fclose (fsac);

/*-------------  calcule de t0  ----------------*/
   {
     int eh, em ,i;
     float fes;
     char koo[9];

     for ( i = 0; i < 8; i++ ) koo[i] = SHD->ko[i];
     koo[8] = '\0';

     SHD->o = SHD->b + SHD->nzhour*3600. + SHD->nzmin*60 +
     SHD->nzsec + SHD->nzmsec*.001;

     sscanf(koo,"%d%*[^0123456789]%d%*[^.0123456789]%g",&eh,&em,&fes);

     SHD->o  -= (eh*3600. + em*60. + fes);
   /*-------------------------------------------*/}

   return SHD;
}


void write_sac (char *fname, float *sig, SAC_HD *SHD)
/*--------------------------------------------------------------------------*/
{
  FILE *fsac;
  int i;
  if((fsac = fopen(fname, "wb"))==NULL) {
    fprintf(stderr,"write_sac: Could not open %s to write\n", fname);
  }
  else {

    if ( !SHD ) {
//      SHD = &SAC_HEADER;
    }

    SHD->iftype = (int)ITIME;
    SHD->leven = (int)TRUE;
    SHD->lovrok = (int)TRUE;
    SHD->internal4 = 6L;
    SHD->depmin = sig[0];
    SHD->depmax = sig[0];

    for ( i = 0; i < SHD->npts ; i++ ) {
      if ( SHD->depmin > sig[i] ) {
        SHD->depmin = sig[i];
      }
      if ( SHD->depmax < sig[i] ) {
        SHD->depmax = sig[i];
      }
    }

    fwrite(SHD,sizeof(SAC_HD),1,fsac);
    fwrite(sig,sizeof(float),(int)(SHD->npts),fsac);

    fclose (fsac);
  }
}
////////////////////////////////////////////////////////////////////////////////////////////

void sac_add ( SAC_HD *shd, float *fsig, char *filename, int flag ) {
  int ns,i,day_num;
  float tsig[1000000];
  SAC_HD shd1;

  day_num = (int)(shd->user0);
  if (read_sac(filename, tsig, &shd1, 100001) == NULL) {
    fprintf(stderr,"cannot open file %s\n",filename);
    }
  ns =(int)shd1.npts;
 
  if (flag == 1) {
    for (i=1;i<=3000;i++) {
	fsig[i+3000] = fsig[i+3000] + tsig[(ns-1)/2 + i];
	fsig[3000-i] = fsig[3000-i] + tsig[(ns-1)/2 - i];
        }
    fsig[3000] = fsig[3000] + tsig[(ns-1)/2];
    *shd = shd1;
    shd->user0 += day_num;
    }
  else if (flag == 0 ) {
    fprintf (stderr,"%s reversed!\n",filename);
    for (i=1;i<=3000;i++) {
        fsig[i+3000] = fsig[i+3000] + tsig[(ns-1)/2 - i];
        fsig[3000-i] = fsig[3000-i] + tsig[(ns-1)/2 + i];
        }
    fsig[3000] = fsig[3000] + tsig[(ns-1)/2];    
    *shd = shd1;
    shd->user0 += day_num;
    }
  }

////////////////////////////////////////////////////////////////////////////////////////////

#define NSTA 5000

int main (int argn, char *argv[]) {
 
  int i,j,k,nsta,ndir,npath,ii,jj,jjj,jjjj,flag;
  short day_fs[32];
  char staname[NSTA][10],outname[100],dirlist[100][300],outfname[300],lst_name[300],outdir1[300],tstr[300],buff[100],dstr[32],pstr[NSTA*NSTA/2][100];
  char tname1[300],tname2[300],tname3[300],tname4[300], stnm1[6], stnm2[6];
  float stalat[NSTA],stalon[NSTA],fsig[6001];
  FILE *ff1, *flst;  

  struct stat st;
  SAC_HD shd,shd1;

  if (argn != 5) {
    fprintf (stderr,"input [station.lst with centre_sta in the 1st line] [dir.lst] [out.dir] [lag]\n");
    return 0;
    }

  sprintf(tstr,"mkdir -p %s\0" , argv[3]);
  system(tstr);

  if ((ff1 = fopen(argv[1],"r"))==NULL) {
    fprintf(stderr,"cannot open file %s\n",argv[1]);
    return 0;
    }
  for (i=0;;i++) {
    if (fscanf(ff1,"%s %g %g", &staname[i], &stalon[i], &stalat[i]) != 3) break;
//    if (strcmp(staname[i],argv[4])==0)icent=i;
    }
  nsta = i;
  fclose(ff1);  

  if ((ff1 = fopen(argv[2],"r"))==NULL) {
    fprintf(stderr,"cannot open file %s\n",argv[2]);
    return 0;
    }
  for (i=0;;i++) {
    if (fscanf(ff1,"%s", &dirlist[i]) != 1) break;
    }
  ndir = i;
  fclose(ff1); 
  
  fprintf (stderr,"%d %d\n",nsta,ndir);

  char tmpc[9], sflag[nsta], day_fc[nsta][9*ndir+1];
  short day_num[nsta];

  for(i=0;i<nsta;i++) day_num[i]=0;
  for (i=0;i<ndir;i++) {
     for(k=0;k<nsta;k++) sflag[k] = '0';
     sprintf(tname1,"%s/Cor_dayflag.lst",dirlist[i]);
     if ((ff1 = fopen(tname1,"r"))==NULL) {
        fprintf(stderr,"cannot locate the dayflag file %s\n",tname1);
        return 0;
     }
     for(j=0;;j++) {
        if( (fgets(pstr[j], 100, ff1)) == NULL ) break; }
     fclose(ff1);
     npath=j;
     for(j=0;j<npath;j++){
        sscanf(pstr[j],"%s %s",&buff,&dstr);
        sscanf(buff,"%[^_]_%[^_]_%[^.]",&tstr,&stnm1,&stnm2);
        if(strcmp(staname[0],stnm1) == 0)
          for(ii=0;ii<nsta;ii++) if(strcmp(stnm2,staname[ii])==0) break;
        else if(strcmp(staname[0],stnm2) == 0)
          for(ii=0;ii<nsta;ii++) if(strcmp(stnm1,staname[ii])==0) break;
        else continue;
        if( ii==nsta ) {
           cout<<"Station "<<stnm1<<" or "<<stnm2<<" not found in "<<argv[1]<<endl;
           continue;
        }
        if( sflag[ii] == '1' ) continue;
        day_fs[0]=0;
        sscanf(dstr,"%1d%1d%1d%1d%1d%1d%1d%1d%1d%1d%1d%1d%1d%1d%1d%1d%1d%1d%1d%1d%1d%1d%1d%1d%1d%1d%1d%1d%1d%1d%1d",&day_fs[1],&day_fs[2],&day_fs[3],&day_fs[4],&day_fs[5],&day_fs[6],&day_fs[7],&day_fs[8],&day_fs[9],&day_fs[10],&day_fs[11],&day_fs[12],&day_fs[13],&day_fs[14],&day_fs[15],&day_fs[16],&day_fs[17],&day_fs[18],&day_fs[19],&day_fs[20],&day_fs[21],&day_fs[22],&day_fs[23],&day_fs[24],&day_fs[25],&day_fs[26],&day_fs[27],&day_fs[28],&day_fs[29],&day_fs[30],&day_fs[31]);
        for(k=1;k<32;k++) day_num[ii] += day_fs[k];
        bin2hex( day_fs, tmpc );
        sprintf(&day_fc[ii][i*9], "%s ", tmpc);
        sflag[ii] = '1';
     }
     for(k=0;k<nsta;k++) {
        if( sflag[k] == '0' ) sprintf(&day_fc[k][i*9],"00000000 "); }
  }

  sprintf(tname1,"%s/Cor_dayflag.lst\0",argv[3]);
  if ((flst = fopen(tname1,"w"))==NULL) {
     fprintf(stderr,"Cannot open Cor_dayflag.lst to write\n");
     return 0;
  }
  for(j=1;j<nsta;j++){
     if( day_num[j] == 0 ) continue;
     sprintf(lst_name,"%s/COR_%s_%s.SAC",staname[0],staname[0],staname[j]);
     fprintf(flst,"%s\t%d\t%s\n",lst_name, day_num[j], day_fc[j]);
  }
  fclose(flst);


  jj = 0; k = 0;

//  for (i=0;i<nsta;i++) {
     sprintf(outdir1,"%s/%s\0",argv[3],staname[0]);
     if (access(outdir1,F_OK) != 1) {
	fprintf (stderr,"open directory %s\n",outdir1);
	sprintf(tstr,"mkdir -p %s\0" , outdir1);
	system(tstr);
	}
     jjj = 0;
     jjjj = 0;

     for (j=0;j<nsta;j++) {
       sprintf(outname,"COR_%s_%s.SAC",staname[0],staname[j]);
       sprintf(outfname,"%s/%s/COR_%s_%s.SAC",argv[3],staname[0],staname[0],staname[j]);
       flag = 0;

       jj = jj +1 ;
       jjj = jjj + 1;
       if (fmod(jj,50) == 0) fprintf(stderr,"%d %d %s %d\n",jj,jjj,outname,jjjj);


       if (access(outfname,F_OK) == 0) continue;

       jjjj = 0;

       for (k=0;k<ndir;k++) {
        if(stat(dirlist[k],&st) != 0){
           printf("Can't access path: %s\n",dirlist[k]);
           continue;
          }

         sprintf(tname1,"%s/%s/COR_%s_%s.SAC\0",dirlist[k],staname[0],staname[0],staname[j]);
         sprintf(tname2,"%s/%s/COR_%s_%s.SAC\0",dirlist[k],staname[j],staname[j],staname[0]);
//	 sprintf(tname3,"%s/COR_%s_%s.SAC\0",dirlist[k],staname[ista1],staname[ista2]);
//	 sprintf(tname4,"%s/COR_%s_%s.SAC\0",dirlist[k],staname[ista2],staname[ista1]);

	 if (access(tname1,F_OK) == 0) {
//	    fprintf (stderr,"yes!! file %s\n",tname1);
	    if (flag == 0) {
		shd = sac_null;
		for (ii=0;ii<6001;ii++) fsig[ii] = 0.;
		sac_add (&shd,fsig,tname1,1);
		flag = 1;
		jjjj ++;
		continue;
		}
	    else {
		sac_add(&shd,fsig,tname1,1);
		jjjj ++;
		continue;
		}
	    }
	 if (access(tname2,F_OK) == 0) {
//	    fprintf (stderr,"yes!! file %s\n",tname2);
            if (flag == 0) {
                shd = sac_null;
                for (ii=0;ii<6001;ii++) fsig[ii] = 0.;
                sac_add (&shd,fsig,tname2,0);
                flag = 1;
		jjjj ++;
		continue;
                }
            else {
                sac_add(&shd,fsig,tname2,0);
		jjjj ++;
		continue;
                }
            }

/*
	 if (access(tname3,F_OK) == 0) {
//	    fprintf (stderr,"yes!! file %s\n",tname3);
            if (flag == 0) {
                shd = sac_null;
                for (ii=0;ii<6001;ii++) fsig[ii] = 0.;
                sac_add (&shd,fsig,tname3,1);
                flag = 1;
		jjjj ++;
		continue;
                }
            else {
                sac_add(&shd,fsig,tname3,1);
		jjjj ++;
		continue;
                }
            }

	  if (access(tname4,F_OK) == 0) {
//	    fprintf (stderr,"yes!! file %s\n",tname4);
            if (flag == 0) {
                shd = sac_null;
                for (ii=0;ii<6001;ii++) fsig[ii] = 0.;
                sac_add (&shd,fsig,tname4,0);
                flag = 1;
		jjjj ++;
                }
            else {
                sac_add(&shd,fsig,tname4,0);
		jjjj ++;
                }
            }
*/
         }
       if (flag==1) {
	// shd.evlo = stalon[0];
	// shd.evla = stalat[0];
	 sprintf(shd.kevnm,"%s\0",staname[0]);
	 shd.stlo = stalon[j];
	 shd.stla = stalat[j];
	 sprintf(shd.kstnm,"%s\0",staname[j]);
	 shd.npts = 6001;
	 shd.b = -3000;
         shd.e = shd.b+shd.npts-1;
	 shd.lcalda = 1;
	 shd.dist = 100.;
	 write_sac(outfname,fsig,&shd);
	 }
//       }
     }

  return 1;
  }
