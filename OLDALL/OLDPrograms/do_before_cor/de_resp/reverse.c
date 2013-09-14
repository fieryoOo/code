/*
 *  reverse.c
 *  
 *
 *  Created by Jiayi Xie on 12/1/10.
 *  Copyright 2010 CU-Boulder. All rights reserved.
 *
 */
#define MAIN
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include "/Users/jiayixie/progs/jy/HEAD/head_noise/64_koftan.h"
#include "/Users/jiayixie/progs/jy/HEAD/head_noise/64_sac_db.h"
#include "/Users/jiayixie/progs/jy/HEAD/head_noise/64_mysac.h"

/*--------------------------------------------------------------------------*/
SAC_HD *read_sac (char *fname, float *sig, SAC_HD *SHD, int nmax)
/*----------------------------------------------------------------------------
 ----------------------------------------------------------------------------*/
{
	FILE *fsac;
	/*..........................................................................*/
	if((fsac = fopen(fname, "rb")) == NULL) {
		printf("could not open sac file %s to write\n",fname);
		return NULL;
	}
	
	if ( !fsac )
	{
		/*fprintf(stderr,"file %s not find\n", fname);*/
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
        koo[8] = 0;
		
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
	if((fsac = fopen(fname, "wb"))==NULL) {
		printf("could not open sac file to write\n");
		exit(1);
	}
	
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

#define SLEN 200000
int reverse(char *sacnm )
{
	float sig0[SLEN],sig1[SLEN];
	int i,leng;
	char exe[200];
	SAC_HD SAC_HEADER;
	
	if(read_sac(sacnm,sig0,&SAC_HEADER,SLEN)== NULL)
	{
		fprintf(stderr,"file %s not found\n",sacnm);
	}
	leng = SAC_HEADER.npts;
	sprintf(exe,"cp %s old_%s",sacnm,sacnm);
	system(exe);
	
	for(i=0;i<leng;i++)
	 {
		 sig1[i]=sig0[i]*(-1);
	 }
	
	write_sac(sacnm,&sig1[0],&SAC_HEADER);
	
	return(1);
}
#define N 60;
int main(int argc, char *argv[])

{
	char sacnm[200],stnm[200];
	FILE *filein;
	int i;
	if(argc!=2)
	{
		fprintf(stderr,"please enter: the stnm filelist contains the stationame whose response file has pi shift\n");
		fprintf(stderr," the reverse file formate is : stnm; the seis_sac name format is ft_stnm.BHZ.SAC\n");
		exit(1);
	}
	
	if((filein=fopen(argv[1],"r"))==NULL)
		{
			fprintf(stderr,"cannot open file for reverse %s\n",argv[1]);
			exit(1);
		}

	for(i=0;;i++)
	{
	if(fscanf(filein,"%s",&(stnm[0]))==EOF)
		break;
		
	/******change here if the name format is different*********/
	sprintf(sacnm,"ft_%s.BHZ.SAC",stnm);
	reverse(sacnm);
	fprintf(stderr,"reversed station:%s SAC:%s\n",stnm,sacnm);
	}
	
	return(1);
}















