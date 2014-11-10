#define MAIN

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../inc/mysac64.h"

/* Function prorotypes */


void filter4_(double *f1,double *f2,double *f3,double *f4,int *npow,
              double *dt,int *n, float seis_inZ[],float seis_outZ[],
              float seis_outampZ[],
              float seis_outphZ[],
	      float seis_inE[],float seis_outE[],
              float seis_outampE[],
              float seis_outphE[],
	      float seis_inN[],float seis_outN[],
              float seis_outampN[],
              float seis_outphN[],
	      int *ns,double *dom);

void swapn(unsigned char *b, int N, int n);




/*c/////////////////////////////////////////////////////////////////////////*/
/*--------------------------------------------------------------------------*/
	SAC_HD *read_sac (char *fname, float *sig, SAC_HD *SHD, int nmax)
/*----------------------------------------------------------------------------
----------------------------------------------------------------------------*/
{
 FILE *fsac;
/*..........................................................................*/
	fsac = fopen(fname, "rb");
	if ( !fsac )
	{
	  //fclose (fsac);
	 return NULL;
	}

	if ( !SHD ) SHD = &SAC_HEADER;

	 fread(SHD,sizeof(SAC_HD),1,fsac);

	 if ( SHD->npts > nmax )
	 {
	  fprintf(stderr,
	   "ATTENTION !!! dans le fichier %s npts est limite a %d",fname,nmax);

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
	koo[8] = NULL;

	SHD->o = SHD->b + SHD->nzhour*3600. + SHD->nzmin*60 +
	 SHD->nzsec + SHD->nzmsec*.001;

	sscanf(koo,"%d%*[^0123456789]%d%*[^.0123456789]%g",&eh,&em,&fes);

	SHD->o  -= (eh*3600. + em*60. + fes);
   /*-------------------------------------------*/}

	return SHD;
}

/*c/////////////////////////////////////////////////////////////////////////*/
/*--------------------------------------------------------------------------*/
	void write_sac (char *fname, float *sig, SAC_HD *SHD)
/*----------------------------------------------------------------------------
----------------------------------------------------------------------------*/
{
 FILE *fsac;
 int i;
/*..........................................................................*/
	fsac = fopen(fname, "wb");

	if ( !SHD ) SHD = &SAC_HEADER;


        SHD->iftype = (int)ITIME;
        SHD->leven = (int)TRUE;

        SHD->lovrok = (int)TRUE;
        SHD->internal4 = 6L;



  /*+++++++++++++++++++++++++++++++++++++++++*/
     SHD->depmin = sig[0];
     SHD->depmax = sig[0];
 
   for ( i = 0; i < SHD->npts ; i++ )
   {
    if ( SHD->depmin > sig[i] ) SHD->depmin = sig[i];
    if ( SHD->depmax < sig[i] ) SHD->depmax = sig[i];
   }

	 fwrite(SHD,sizeof(SAC_HD),1,fsac);

	 fwrite(sig,sizeof(float),(int)(SHD->npts),fsac);


	 fclose (fsac);
}



/*c/////////////////////////////////////////////////////////////////////////*/
float sigZ[1000000],sigE[1000000],sigN[1000000];
SAC_HD shdZ,shdE,shdN;
float smoothZ[1000000],smoothE[1000000],smoothN[1000000],smooth[1000000];
SAC_HD sZ,sE,sN;

/*c/////////////////////////////////////////////////////////////////////////*/

int main (int argc, char *argv[])
{
static int n, ns,npow;
static double f1, f2, f3, f4, dt,dom;
static float seis_inZ[400000],seis_outZ[400000];
static float seis_outampZ[400000],seis_outphZ[400000];
static float seis_inE[400000],seis_outE[400000];
static float seis_outampE[400000],seis_outphE[400000];
static float seis_inN[400000],seis_outN[400000];
static float seis_outampN[400000],seis_outphN[400000];
double temp;
double t1,t2,t3,t4,tbeg,tend;
char  name[160], name1[160];
 char nameZ[160],nameE[160],nameN[160];
char  nameamp[160],nameph[160];
FILE  *in, *ff;
int   i, j, nn;


  if( argc != 2) {
      printf("Usage: whiten_phamp  parameter_file\n");
      exit(-1);
  }

// open and read parameter file param.dat
  if((in = fopen(argv[1],"r")) == NULL) {
      printf("Can not find file %s.\n",argv[1]);
      exit(1);
  }

  while((nn = fscanf(in,"%lf %lf %lf %lf %lf %d %s %lf %lf",&t1,&t2,&t3,&t4,
            &dt,&npow,name,&tbeg,&tend)) != EOF) { // start main loop
      if(nn == 0 || nn != 9) break;
      printf("Corners periods. Low: %f - %f, High: %f - %f\n",t1, t2, t3, t4);
      printf("Step: %f, Cosine power: %d\n",dt, npow);
      printf("time segment: from %lf to %lf\n",tbeg, tend);

// remove quotes from name
      j = 0;
      for(i = 0; i < strlen(name); i++) {
          if(name[i] == '\'' || name[i] == '\"') continue;
          name[j] = name[i]; j++;
      }
      name[j] = '\0';
      sprintf(nameZ,"%s.LHZ.SAC",name);
      sprintf(nameE,"%s.LHE.SAC",name);
      sprintf(nameN,"%s.LHN.SAC",name);

// do running average before whitening

  ff = fopen("sac_one_cor","w");
  fprintf(ff,"rm aZ.avg aE.avg aN.avg\n");
  fprintf(ff,"sac << END\n");
  fprintf(ff,"r eqf%s eqf%s eqf%s\n",nameZ,nameE,nameN);
  fprintf(ff,"abs\n");
  fprintf(ff,"smooth mean h 128\n");
  fprintf(ff,"w aZ.avg aE.avg aN.avg\n");
  fprintf(ff,"quit\n");
  fprintf(ff,"END\n");
  fclose(ff);
  system("sh sac_one_cor");
  if ( !read_sac("aZ.avg", smoothZ, &sZ, 1000000 ) )
    {
      fprintf(stderr,"file %s did not found\n", "aZ.avg" );
      continue;
    }
  if ( !read_sac("aE.avg", smoothE, &sE, 1000000 ) )
    {
      fprintf(stderr,"file %s did not found\n", "aE.avg" );
      continue;
    }
  if ( !read_sac("aN.avg", smoothN, &sN, 1000000 ) )
    {
      fprintf(stderr,"file %s did not found\n", "aN.avg" );
      continue;
    }
  n  = sZ.npts;
  if(sZ.npts!=sE.npts||sZ.npts!=sN.npts)
    {
      fprintf(stderr,"format error ZEN %s!!\n",name);
      continue;
    }
  for(i=0;i<n;i++)
    {
      temp=smoothZ[i];
      if(smoothE[i]>temp)
	temp=smoothE[i];
      if(smoothN[i]>temp)
	temp=smoothN[i];
      smooth[i]=temp;
    }
  write_sac("a1.avg",smooth, &sZ);
    
  //fprintf(ff,"r aZ.avg \n");
  //fprintf(ff,"mulf aN.avg \n");
  //fprintf(ff,"mulf aE.avg \n");
  //fprintf(ff,"mulf aE.avg \n");
  //fprintf(ff,"w a1.avg\n");
  //fprintf(ff,"cut 500 82500\n");
  //fprintf(ff,"r a1.avg\n");
  //fprintf(ff,"w a1.avg\n");
  // a2.avg a3.avg\n");
  //fprintf(ff,"r a1.avg\n");
  //fprintf(ff,"mulf a2.avg\n");
  //fprintf(ff,"mulf a3.avg\n");
  //fprintf(ff,"smooth h 128\n");
  //fprintf(ff,"w a4.avg\n");
  ff = fopen("sac_one_cor","w");
  fprintf(ff,"rm smoothZ.sac smoothE.sac smoothN.sac\n");
  fprintf(ff,"sac << END\n");
  fprintf(ff,"r %s %s %s\n",nameZ,nameE,nameN);
  fprintf(ff,"divf a1.avg a1.avg a1.avg\n");
  fprintf(ff,"w smoothZ.sac smoothE.sac smoothN.sac\n");
  fprintf(ff,"cut %lf %lf\n",tbeg,tend);
  fprintf(ff,"r smoothZ.sac smoothE.sac smoothN.sac\n");
  fprintf(ff,"w over\n");
  fprintf(ff,"cut off\n");
  fprintf(ff,"quit\n");
  fprintf(ff,"END\n");
  fclose(ff);
  system("sh sac_one_cor");

  // end of running average


  if ( !read_sac("smoothZ.sac", sigZ, &shdZ, 1000000 ) )
    {
      fprintf(stderr,"file %s did not found\n", nameZ );
      continue;
    }

  if ( !read_sac("smoothE.sac", sigE, &shdE, 1000000 ) )
    {
      fprintf(stderr,"file %s did not found\n", nameE );
      continue;
    }
  if ( !read_sac("smoothN.sac", sigN, &shdN, 1000000 ) )
    {
      fprintf(stderr,"file %s did not found\n", nameN );
      continue;
    }
  n  = shdZ.npts;
  dt = shdZ.delta;
  if(shdZ.npts!=shdE.npts||shdZ.npts!=shdN.npts||shdZ.delta!=shdE.delta||shdZ.delta!=shdN.delta)
    {
      fprintf(stderr,"format error ZEN %s!!\n",name);
      continue;
    }
     for( i =0; i< n; i++)
     {  
     seis_inZ[i] = sigZ[i];
     seis_inE[i] = sigE[i];
     seis_inN[i] = sigN[i];
     //     printf(" seis_in1  %d %f\n", i,sig[i]);
     }

      printf(" Dt1= %f, Nsamples1= %d\n",dt, n);


      f1 = 1.0/t1; f2 = 1.0/t2; f3 = 1.0/t3; f4 = 1.0/t4;
      filter4_(&f1,&f2,&f3,&f4,&npow,&dt,&n,seis_inZ,seis_outZ,seis_outampZ,seis_outphZ,seis_inE,seis_outE,seis_outampE,seis_outphE,seis_inN,seis_outN,seis_outampN,seis_outphN,&ns,&dom);
      
      sprintf(name1,"wt%s",nameZ);
      write_sac(name1,seis_outZ, &shdZ);
      strcpy(nameamp,name1);
      strcpy(nameph,name1);
      strcat(nameamp,".am");
      strcat(nameph, ".ph");
      shdZ.npts = ns/2 + 1;
      shdZ.delta = dom;
      shdZ.b = 0;
      shdZ.iftype = IXY;
      write_sac(nameamp,seis_outampZ, &shdZ );
      write_sac(nameph, seis_outphZ,  &shdZ );
      
      sprintf(name1,"wt%s",nameE);
      write_sac(name1,seis_outE, &shdE);
      strcpy(nameamp,name1);
      strcpy(nameph,name1);
      strcat(nameamp,".am");
      strcat(nameph, ".ph");
      shdE.npts = ns/2 + 1;
      shdE.delta = dom;
      shdE.b = 0;
      shdE.iftype = IXY;
      write_sac(nameamp,seis_outampE, &shdE );
      write_sac(nameph, seis_outphE,  &shdE );

      sprintf(name1,"wt%s",nameN);
      write_sac(name1,seis_outN, &shdN);
      strcpy(nameamp,name1);
      strcpy(nameph,name1);
      strcat(nameamp,".am");
      strcat(nameph, ".ph");
      shdN.npts = ns/2 + 1;
      shdN.delta = dom;
      shdN.b = 0;
      shdN.iftype = IXY;
      write_sac(nameamp,seis_outampN, &shdN );
      write_sac(nameph, seis_outphN,  &shdN );

  }
  
   return 0;
  }
