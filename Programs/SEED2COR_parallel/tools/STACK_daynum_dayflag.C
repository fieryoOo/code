#define MAIN

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

#include <iostream>
#include <sys/stat.h>
#include "/projects/yeti4009/code/Programs/head/mysac.h"
//#include "/home/tianye/code/Programs/head/mysac.h"

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


/*c/////////////////////////////////////////////////////////////////////////*/
/*--------------------------------------------------------------------------*/
SAC_HD *read_sac (char *fname, float **sig, SAC_HD *SHD) {
	FILE *fsac;
	if((fsac = fopen(fname, "r"))==NULL) return NULL;
	if ( !SHD ) SHD = &SAC_HEADER;
	fread(SHD,sizeof(SAC_HD),1,fsac);
	*sig = (float *) malloc (SHD->npts * sizeof(float));
	fread(*sig,sizeof(float),SHD->npts,fsac);
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
	fsac = fopen(fname, "wb");

	if ( !SHD ) SHD = &SAC_HEADER;


	SHD->iftype = (int)ITIME;
	SHD->leven = (int)TRUE;

	SHD->lovrok = (int)TRUE;
	SHD->internal4 = 6L;



	/*+++++++++++++++++++++++++++++++++++++++++*/
	SHD->depmin = sig[0];
	SHD->depmax = sig[0];

	for ( i = 0; i < SHD->npts ; i++ ) {
		if ( SHD->depmin > sig[i] ) SHD->depmin = sig[i];
		if ( SHD->depmax < sig[i] ) SHD->depmax = sig[i];
	}
	fwrite(SHD,sizeof(SAC_HD),1,fsac);

	fwrite(sig,sizeof(float),(int)(SHD->npts),fsac);


	fclose (fsac);
}



////////////////////////////////////////////////////////////////////////////////////////////

void sac_add ( SAC_HD *shd, float **fsig, char *filename, int flag ) {
	int ns,i, lag, lag1;
	float *tsig;
	SAC_HD shd1;

	cout<<"sac_add: "<<filename<<endl;
	if(*fsig==NULL) {
		if(flag==1) {
			read_sac(filename, fsig, shd);
			if( shd->user0 == 0 ) {
				free(*fsig);
				*fsig=NULL;
				return;
			}
		}
		else {
			fprintf (stderr,"Warning: Reversed signal, skipped!\n");
			/*
			fprintf (stderr,"Warning: Signal reversed!\n");
			float *fsigtmp, ftmp;
			char ctmp[10];
			read_sac(filename, &fsigtmp, shd);
			if( shd->user0 == 0 ) {
				free(fsigtmp);
				return;
			}
			*fsig = (float *) malloc (shd->npts * sizeof(float));
			sprintf(ctmp, "%s", shd->kstnm); sprintf(shd->kstnm, "%s", shd->kevnm); sprintf(shd->kevnm, "%s", ctmp);
			ftmp = shd->stla; shd->stla = shd->evla; shd->evla = ftmp;
			ftmp = shd->stlo; shd->stlo = shd->evlo; shd->evlo = ftmp;
			ftmp = shd->az; shd->az = shd->baz; shd->baz = ftmp;
			lag = (shd->npts-1)/2;
			for (i=1;i<=lag;i++) {
				(*fsig)[lag+i] = fsigtmp[lag-i];
				(*fsig)[lag-i] = fsigtmp[lag+i];
			}
			(*fsig)[lag] = fsigtmp[lag];
			free(fsigtmp);
			*/
		}
		if( shd->user0 == 0 ) {
			free(*fsig);
			*fsig = NULL;
		}
		return;
	}

	if (read_sac(filename, &tsig, &shd1) == NULL) {
		fprintf(stderr,"Error: cannot open file %s\n",filename);
		return;
	}
	if( shd1.user0 == 0 ) {
		free(tsig);
		return;
	}
	ns = shd1.npts;

	if(shd->delta != shd1.delta) {
		fprintf(stderr,"Error: incompatible dt in file %s\n",filename);
		return;
	}

	lag = (shd->npts-1)/2;
	lag1 = (ns-1)/2;
	if(lag!=lag1) fprintf(stderr,"Warning: incompatible lag time in file %s\n",filename);
	if (flag == 1) {
		for (i=1;i<=lag;i++) {
			(*fsig)[lag+i] = (*fsig)[lag+i] + tsig[lag1+i];
			(*fsig)[lag-i] = (*fsig)[lag-i] + tsig[lag1-i];
		}
		(*fsig)[lag] = (*fsig)[lag] + tsig[lag1];
		shd->user0 += shd1.user0;
	}
	else if (flag == 0 ) {
		fprintf (stderr,"Warning: Reversed signal, skipped!\n");
		/*
		fprintf (stderr,"Warning: signal reversed!\n");
		for (i=1;i<=lag;i++) {
			(*fsig)[lag+i] = (*fsig)[lag+i] + tsig[lag1-i];
			(*fsig)[lag-i] = (*fsig)[lag-i] + tsig[lag1+i];
		}
		(*fsig)[lag] = (*fsig)[lag] + tsig[lag1];    
		shd->user0 += shd1.user0;
		*/
	}
	free(tsig);
}

////////////////////////////////////////////////////////////////////////////////////////////

#define NSTA 2000

int main (int argn, char *argv[]) {

	int i,j,k,nsta,ndir,npath,ii,jj,jjj,jjjj,flag;
	short day_fs[32];
	char staname[NSTA][10],outname[100],dirlist[100][300],outfname[300],lst_name[300],outdir1[300],tstr[300],buff[300],dstr[32];
	char tname1[300],tname2[300],tname3[300],tname4[300], stnm1[6], stnm2[6];
	float stalat[NSTA],stalon[NSTA],*fsig;
	FILE *ff1, *flst;  

	struct stat st;
	SAC_HD shd,shd1;

	if (argn != 5) {
		fprintf (stderr,"input [sta.lst (sta lon lat)] [dir.lst (path/dir)] [type(COR or DCV)] [out.dir]\n");
		return 0;
	}

	char type[4];
	sprintf(type,"%s",argv[3]);
	if(strcmp(type,"COR")!=0 && strcmp(type,"DCV")!=0) {
		cout<<"Unrecognized data type: "<<type<<endl;
		exit(0);
	}

	sprintf(tstr,"mkdir -p %s\0" , argv[4]);
	system(tstr);

	// read in stations
	if ((ff1 = fopen(argv[1],"r"))==NULL) {
		fprintf(stderr,"cannot open file %s\n",argv[1]);
		return 0;
	}

	for (i=0;i<NSTA&&fgets(buff, 300, ff1)!=NULL;i++)
		if (sscanf(buff,"%s %g %g", &staname[i], &stalon[i], &stalat[i]) != 3) {
			cout<<"Warning: format error in stations list: "<<buff<<". Skipped!"<<endl;
			continue;
		}
	nsta = i;
	fclose(ff1);
	if( i >= NSTA ) {
		fprintf(stderr, "num of stations exceeds the limit! Increase NSTA!\n");
		exit(0);
	}

	// read in directories
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

	char tmpc[9], day_fc[nsta][nsta][9*ndir+1];
	int sflag[nsta][nsta], daynum[nsta][nsta];

	char *pstr[NSTA*NSTA]; // define 2 times of the possible max # of paths to be safe
	for(i=0; i<NSTA*NSTA; i++) pstr[i] = new char[300];



	/*
		char ***day_fc = (char ***) malloc(nsta * sizeof(char **));
		for(i=0;i<nsta;i++) {day_fc[i] = (char **) malloc(nsta * sizeof(char *));}
		for(i=0;i<nsta;i++) for(j=0;j<nsta;j++) {
		day_fc[i][j] = (char *) malloc((9 * ndir + 1) * sizeof(char)); }
		char **sflag = (char **) malloc(nsta * sizeof(char *));
		for(i=0;i<nsta;i++) {sflag[i] = (char *) malloc(nsta * sizeof(char));}
		if( sflag == NULL || day_fc == NULL) return 0;
		for(i=0;i<nsta;i++) {if( sflag[i] == NULL || day_fc[i] == NULL )return 0;}
		for(i=0;i<nsta;i++) for(j=0;j<nsta;j++) {if( day_fc[i][j] == NULL )return 0;}
		*/
	for(j=0;j<nsta;j++)for(k=j;k<nsta;k++) daynum[j][k] = 0;

	for (i=0;i<ndir;i++) {
		for(j=0;j<nsta;j++)for(k=j;k<nsta;k++) sflag[j][k] = -1;
		sprintf(tname1,"%s/Cor_dayflag.lst",dirlist[i]);
		if ((ff1 = fopen(tname1,"r"))==NULL) {
			fprintf(stderr,"cannot locate the dayflag file %s\n",tname1);
			return 0;
		}
		for(j=0;;j++) {
			if( (fgets(pstr[j], 300, ff1)) == NULL ) break; 
		}
		fclose(ff1);
		npath=j;
		for(j=0;j<npath;j++){
			sscanf(pstr[j],"%s%s",&buff,&dstr);
			sscanf(buff,"%[^_]_%[^_]_%[^.]",&tstr,&stnm1,&stnm2);
			for(ii=0;ii<nsta;ii++) if(strcmp(stnm1,staname[ii])==0) break;
			for(jj=0;jj<nsta;jj++) if(strcmp(stnm2,staname[jj])==0) break;
			if( ii==nsta || jj==nsta ) {
				cout<<"Station "<<stnm1<<" or "<<stnm2<<" not found in "<<argv[1]<<endl;
				continue;
			}
			if(ii>jj) {
				k=ii; ii=jj; jj=k;
			}
			if( sflag[ii][jj] != -1 ) continue;
			day_fs[0]=0;
			sscanf(dstr,"%1d%1d%1d%1d%1d%1d%1d%1d%1d%1d%1d%1d%1d%1d%1d%1d%1d%1d%1d%1d%1d%1d%1d%1d%1d%1d%1d%1d%1d%1d%1d",&day_fs[1],&day_fs[2],&day_fs[3],&day_fs[4],&day_fs[5],&day_fs[6],&day_fs[7],&day_fs[8],&day_fs[9],&day_fs[10],&day_fs[11],&day_fs[12],&day_fs[13],&day_fs[14],&day_fs[15],&day_fs[16],&day_fs[17],&day_fs[18],&day_fs[19],&day_fs[20],&day_fs[21],&day_fs[22],&day_fs[23],&day_fs[24],&day_fs[25],&day_fs[26],&day_fs[27],&day_fs[28],&day_fs[29],&day_fs[30],&day_fs[31]);
			bin2hex( day_fs, tmpc );
			sprintf(&day_fc[ii][jj][i*9], "%s ", tmpc);
			for(k=1,sflag[ii][jj]=0;k<32;k++) sflag[ii][jj] += day_fs[k];
			daynum[ii][jj] += sflag[ii][jj];
			/*
				if(strlen(day_fc[ii][jj])!=(i+1)*9){
				for(k=0;k<(i+1)*9;k++)if(day_fc[ii][jj][k]=='\0')day_fc[ii][jj][k]='!';
				cout<<strlen(day_fc[ii][jj])<<" ?= "<<(i+1)*9<<endl;
				ntmp=day_fc[ii][jj];
				cout<<"after: "<<tmpc<<" "<<ntmp<<endl;
				}
				*/
		}
		for(j=0;j<nsta;j++) for(k=j+1;k<nsta;k++)
			if( sflag[j][k] == -1 ) sprintf(&day_fc[j][k][i*9],"00000000 "); 
	}
	for(i=0; i<nsta*nsta; i++) delete [] pstr[i];

	sprintf(tname1,"%s/Cor_dayflag.lst\0",argv[4]);
	if ((flst = fopen(tname1,"w"))==NULL) {
		fprintf(stderr,"Cannot open Cor_dayflag.lst to write\n");
		return 0;
	}
	for(i=0;i<nsta;i++) for(j=i+1;j<nsta;j++){
		sprintf(lst_name,"%s/%s_%s_%s.SAC",staname[i],type,staname[i],staname[j]);
		fprintf(flst,"%-30s\t%s\t%4d\n",lst_name, day_fc[i][j], daynum[i][j]);
	}
	fclose(flst);

	//  sprintf(tname1,"%s/file_daynum.lst",argv[4]);
	//  if ((flst = fopen(tname1,"w"))==NULL) {
	//     fprintf(stderr,"Cannot open list file file_daynum.lst to write\n");
	//     return 0;
	//    }

	jj=0; k=0;
	for (i=0;i<nsta;i++) {
		sprintf(outdir1,"%s/%s\0",argv[4],staname[i]);
		if (access(outdir1,F_OK) != 1) {
			fprintf (stderr,"open directory %s\n",outdir1);
			sprintf(tstr,"mkdir -p %s\0" , outdir1);
			system(tstr);
		}
		jjj = 0;
		jjjj = 0;

		for (j=i+1;j<nsta;j++) {
			sprintf(outname,"%s_%s_%s.SAC",type,staname[i],staname[j]);
			sprintf(lst_name,"%s/%s_%s_%s.SAC",staname[i],type,staname[i],staname[j]);
			sprintf(outfname,"%s/%s/%s_%s_%s.SAC",argv[4],staname[i],type,staname[i],staname[j]);
			//flag = 0;
			fsig = NULL;

			jj = jj +1 ;
			jjj = jjj + 1;
			if (fmod(jj,50) == 0) fprintf(stderr,"%d %d %s %d\n",jj,jjj,outname,jjjj);


			//if (access(outfname,F_OK) == 0) continue;

			jjjj = 0;
			for (k=0;k<ndir;k++) {
				if(stat(dirlist[k],&st) != 0){
					printf("Can't access path: %s\n",dirlist[k]);
					continue;
				}
				sprintf(tname1,"%s/%s/%s_%s_%s.SAC\0",dirlist[k],staname[i],type,staname[i],staname[j]);
				sprintf(tname2,"%s/%s/%s_%s_%s.SAC\0",dirlist[k],staname[j],type,staname[j],staname[i]);
				//	 sprintf(tname3,"%s/%s_%s_%s.SAC\0",dirlist[k],type,staname[i],staname[j]);
				//	 sprintf(tname4,"%s/%s_%s_%s.SAC\0",dirlist[k],type,staname[j],staname[i]);

				if (access(tname1,F_OK) == 0) {
					sac_add(&shd,&fsig,tname1,1);
					jjjj++;
					/*	    fprintf (stderr,"yes!! file %s\n",tname1);
							 if (flag == 0) {
							 shd = sac_null;
							 shd.user0=0;
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
							 */
				}
				if (access(tname2,F_OK) == 0) {
					sac_add (&shd,&fsig,tname2,0);
					jjjj++;
					/*	    fprintf (stderr,"yes!! file %s\n",tname2);
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
							 */
				}

			}
			if (fsig!=NULL) {
				//shd.evlo = stalon[i];
				//shd.evla = stalat[i];
				//sprintf(shd.kevnm,"%s\0",staname[i]);
				//shd.stlo = stalon[j];
				//shd.stla = stalat[j];
				//sprintf(shd.kstnm,"%s\0",staname[j]);
				//shd.npts = 6001;
				//shd.b = -3000;
				shd.lcalda = 1;
				//shd.dist = 100.;
				write_sac(outfname,fsig,&shd);
				//         fprintf(flst,"%s   %d\n",lst_name,(int)shd.user0);
				free(fsig);
			}

		} //sta2
	} //sta1
	//  fclose(flst);

	return 1;
}
