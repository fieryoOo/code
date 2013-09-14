#define MAIN
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include "/home/weisen/progs/NOISE_CODA/HEADFILE/mysac.h"
//#include "/Users/jiayixie/progs/NOISE_CODA/HEAD_NOISE/64_mysac.h"
#include "/home/tianye/code/Programs/head/64_mysac.h"
#include "string.h"
/* a program to take the list files, for up to 30 files, 
   and stack all the files that are present in all 
   directories.  The list files must have the absolute 
   path as the first line.  The other lines should be 
   a file name followed by a tab then the number of days 
   of raw data that lead to this correlation.  		*/


  /* information about all files in all directories are 
     read into this af[] structure.   */
  // change the delta of sac file when it is not 1.0
/* _cv vision dosn't use the SAC program but use the pure C progrrame.
 *
 */

  struct onefile { /* file name and number */
    char name[100]; /* complete file name */
    int  ndays; /* number of days for the given correlation */
    int  onoff; /* tells if file has been stacked (1) or not (0) */
    float depmin; /* min signal value, from sac header depmin value */
    float depmax; /* like depmin, max signal value */
  };

  struct allfiles { /*for whole filelist with dir and len*/
    char dir[100]; /* absolute path, ending in a "/" */
    struct onefile fname[400000]; /* structure for each file */
    int len; /* total number of files in a given directory */
  } af[100]; /* up to 30 month.lst files can be used, this can be changed */

/**/

        void write_sac (char *fname, float *sig, SAC_HD *SHD)
/*----------------------------------------------------------------------------
 * ----------------------------------------------------------------------------*/
{
 FILE *fsac;
 int i;
/*..........................................................................*/
        fsac = fopen(fname, "wb");

        if ( !SHD ) SHD = &SAC_HEADER;


//       SHD->iftype = (long)ITIME;
       SHD->iftype = (int)ITIME;
 //       SHD->leven = (long)TRUE;
        SHD->leven = (int)TRUE;

//        SHD->lovrok = (long)TRUE;
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

 /**/
/*c/////////////////////////////////////////////////////////////////////////*/
/*--------------------------------------------------------------------------*/
        SAC_HD *read_sac (char *fname, float *sig, SAC_HD *SHD, int nmax)
/*----------------------------------------------------------------------------
----------------------------------------------------------------------------*/
{
 FILE *fsac;
/*..........................................................................*/
        if((fsac = fopen(fname, "rb")) == NULL) {
          printf("could not open sac file to read%s \n", fname);
          exit(1);
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
        koo[8] = 0 ;

        SHD->o = SHD->b + SHD->nzhour*3600. + SHD->nzmin*60 +
         SHD->nzsec + SHD->nzmsec*.001;

        sscanf(koo,"%d%*[^0123456789]%d%*[^.0123456789]%g",&eh,&em,&fes);

        SHD->o  -= (eh*3600. + em*60. + fes);
   /*-------------------------------------------*/}

        return SHD;
}

#define SLEN 400000
float sig0[SLEN];
float cvsig[SLEN];
float tempsig[SLEN];

int main(int argc, char *argv[])
{

  FILE *fp1, *fp0;
  char ch;
  int i = 0, j = 0, k = 0, quit = 0, nfiles = 0;
  int inotbase = 0, jbase = 0, jnotbs = 0;
  int div_by = 0, mindays = 0, ibase = 0;
  char dir[200], filename[200], nombre[200];
  char dump[200], outdir[200], fname[200];

  float minv = 2000000;
  SAC_HD cvSAC_HEADER, tempSAC_HEADER;
  /* check for comand line arg */
  if(argc!=4) {
    printf("Usage: enter filelist outfile_directory mindays\n");
    exit(1);
  }
 
  strcpy(outdir, argv[2]); 

  /* open file for reading */
  if((fp0 = fopen(argv[1], "r"))==NULL) {
    printf("cannot open infile.\n");
    exit(1);
  }
  fgets(filename, 60, fp0); /*name of month.lst */  

  /* fill the structure with all the file information for up to 
	30 files */  //read from the inputfile to get the the month.lst(af[k])
  do{             
    i = strlen(filename);
    filename[i-1] ='\0';
    printf("file name %d is %s\n", k, filename);

    /* open file for reading */
    if((fp1 = fopen(filename, "r"))==NULL) {    //read the month.lst
      printf("cannot open infile number %d.\n", k);
      exit(1);
    }

    /* first line of file must be the absolute path*/
    fgets(dir, 100, fp1);
    i = strlen(dir);
    dir[i-1] = '\0';
    strcpy(af[k].dir, dir);
    //fprintf(stderr,"the SAC files are at : %s \n",dir);

    if( fgets(nombre,100,fp1)==NULL) //80  //read from month.lst to get name of the COR
        goto next;
    do {								// initialize the af[k].fname[j].ndays/onoff/depmin/depmax
      i = strlen(nombre);
      nombre[i-1] = '\0';
//	printf("/// dir %s //// nombre %s\n",af[k].dir,nombre);
      sprintf(fname, "%s%s", af[k].dir, nombre);
      if (fmod(j,50.0)==0.0)                    //????j
      fprintf(stderr,"fname : %s%s\n",af[k].dir,nombre);
      /*---------------- reading sac file  -------------------*/
      if ( read_sac (fname, sig0, &SAC_HEADER, SLEN) == NULL ){
          fprintf(stderr,"file %s not found\n", fname);
          return 0;
        }

      /* judge if the sac file is a null file with only headfile no signal */     
     if( SAC_HEADER.depmin < minv ) {           //???? if no signal, the depmin will be very big? It seems if no sig the depmin will be -12345.0
       strcpy(af[k].fname[j].name, nombre);  //name of the COR data
      /* get the number of days */
       //      fgets(dump, 30, fp1);
       //af[k].fname[j].ndays = atoi(dump);
      if (SAC_HEADER.user0<=31 && SAC_HEADER.user0 >= 0)  //????user0 and the following 20
       af[k].fname[j].ndays = SAC_HEADER.user0;
      else
       af[k].fname[j].ndays = 20;
      af[k].fname[j].onoff = 0;
      af[k].fname[j].depmin = ((-1)*SAC_HEADER.depmin);  //why -1
      af[k].fname[j].depmax = SAC_HEADER.depmax;
      j++;
     } 
     else {
       //fgets(dump, 30, fp1);
     }   
     fgets(nombre,100,fp1);//80
    } while(!feof(fp1));
    next:

    /* j is number of files for a given dir */ //the month.lst
    af[k].len = j; 
    fprintf(stderr, "length of file %d is %d\n", k, af[k].len);
    j = 0;
    k++;
  
    fgets(filename, 100, fp0); //60
  } while(!feof(fp0));   
  fprintf(stderr, "!!!all info read fine\n");
	fprintf(stderr,"next sentence!!!\n");
//	abort();
	
  nfiles = k; /* number of files as determined from loop*/  //file month.lst

  fclose(fp1);
  fclose(fp0);
  /*next section checks for file existence in all directories
    and writes script to add the files together.  It will only be 
    executed if the file is in all directories.    		*/  
 
  /* inotbase is file num.  inotbase=0 is reference file*/
  jbase = 0;/* j is line of first (ibase) file	*/
  jnotbs = 0;/* jnotbs is line of second file 		*/
  char cvsacname[300];
  char tempsacname[300];
  char cvobsac[300];
  int cvi=0;
  /* increment over nfiles so if a file is in the second, third
	etc but not the first it will not be missed */
  for(ibase = 0; ibase < nfiles; ibase++){
    if (af[ibase].len==0)   //num of COR in the month.lst
      continue;
   
do{  /* outer do while */
	inotbase = ibase + 1;
	
      if(af[ibase].fname[jbase].onoff==0){  //have not been stacked
        /* write sac script header */
        /*if((fp1 = fopen("dostack.csh", "w"))==NULL) {
          fprintf(stderr, "cannot open dostack.csh \n");
          exit(1);
        }
        fprintf(stderr, "filename is %s\n", af[ibase].fname);*/
        // read sac file to cvsig[200000], SAC_HD cvsac_hd;
		  sprintf(cvsacname,"%s%s",af[ibase].dir, af[ibase].fname[jbase].name);  //pass to the CORsac file
	fprintf(stderr,"cv !! %s\n",cvsacname);
//		  abort();
		  
	if (  read_sac (cvsacname, cvsig, &cvSAC_HEADER, SLEN) == NULL ){
          fprintf(stderr,"file %s not found\n", af[ibase].fname);
          return 0;
        }
		  
        if (cvSAC_HEADER.delta!=1.0) cvSAC_HEADER.delta=1.0;
        af[ibase].fname[jbase].onoff=1;
        //fprintf(fp1, "#!/bin/csh\n");
        //fprintf(fp1, "sac2000 << END\n");
 
        //fprintf(fp1, "r %s%s\n", af[ibase].dir, af[ibase].fname[jbase].name);

	//        if(af[ibase].fname[jbase].depmin > af[ibase].fname[jbase].depmax){
	//   fprintf(fp1, "div %f\n", af[ibase].fname[jbase].depmin);
	// }
	// else fprintf(fp1, "div %f\n", af[ibase].fname[jbase].depmax);
	//fprintf(fp1, "mul %d\n", af[ibase].fname[jbase].ndays);

        //fprintf(fp1, "w a\n\n");
        /* div_by will sum up total number of days from each month	*/
        div_by = af[ibase].fname[jbase].ndays; 
  
        do{ /* inner do while */
		
          /* if strings are different, increment jnotbs, quit if jnotbs
     	  gets to reach the end of the file			*/
          if(strcmp(af[ibase].fname[jbase].name, af[inotbase].fname[jnotbs].name)) { 
			if(jnotbs == af[inotbase].len) {
//              fprintf(stderr,"%d %d\n",jnotbs, af[inotbase].len);
	      fprintf(stderr, "match wrong\n");
	      if(inotbase != nfiles) inotbase++;
	      jnotbs = 0;
  	    }
            else jnotbs++;
          }

          /* if file names are the same, write this part of the script*/
          else {
    	    //fprintf(fp1, "r %s%s\n", af[inotbase].dir, af[ibase].fname[jbase].name);
    	    sprintf(tempsacname,"%s%s", af[inotbase].dir, af[ibase].fname[jbase].name);
	    fprintf(stderr,"%s\n",tempsacname);
    	    if (  read_sac (tempsacname, tempsig, &tempSAC_HEADER, SLEN) == NULL ){
              fprintf(stderr,"file %s not found\n", af[ibase].fname); //????why not fanme[].name
              return 0;
            }
            if (tempSAC_HEADER.delta!=1.0) tempSAC_HEADER.delta=1.0;
	    af[inotbase].fname[jnotbs].onoff = 1;

	    //          if(af[inotbase].fname[jnotbs].depmin > af[inotbase].fname[jnotbs].depmax){
	    //	 fprintf(fp1, "div %f\n", af[inotbase].fname[jnotbs].depmin);
	    // }
            //else fprintf(fp1, "div %f\n", af[inotbase].fname[jnotbs].depmax);
	    //fprintf(fp1, "mul %d\n", af[inotbase].fname[jnotbs].ndays);

            //fprintf(fp1, "w aa\n");
	    //fprintf(fp1, "r a\n");
	    //fprintf(fp1, "addf aa\n");
	    //fprintf(fp1, "w a\n\n");
	    fprintf(stderr,"npts: %d\n",tempSAC_HEADER.npts);
	    if (tempSAC_HEADER.npts == cvSAC_HEADER.npts) {
		fprintf(stderr,"now do one!!!\n");
              for (cvi=0; cvi<tempSAC_HEADER.npts; cvi++) {
                    cvsig[cvi]=cvsig[cvi]+tempsig[cvi];
	      } 
              div_by += af[inotbase].fname[jnotbs].ndays;   ///????
	    }
	    else fprintf(stderr,"NPTS the sac file %s is not consistant with sac file %s\n", cvsacname,tempsacname);
            inotbase++;  //go to the next month.lst to find the match CORsac
    	    jnotbs = 0;
          }
        } while(inotbase != nfiles); /* inner do while */
        //fprintf(fp1, "r a\n");
        //fprintf(fp1, "ch user0 %d\n", div_by);
        //fprintf(fp1, "w %s%s\n", outdir, af[ibase].fname[jbase].name); 
        //fprintf(fp1, "END\n\n");
        //fclose(fp1);

        /* if inotbase = nfiles, file present in all directories */
        /* this could be made into something less than nfiles
    	  for example, if nfiles = 12, this could be executed
  	  for inotbase > 10 or something. 			*/
        if(inotbase == nfiles){  
          //fprintf(stderr, "file %s present in all dirs \n", af[ibase].fname[jbase].name);
        }
        strcpy(dump, argv[3]);
        mindays=atoi(dump);
        if(div_by >= mindays){
          //system("csh dostack.csh"); 
	  //sprintf(dump, "mv dostack.csh scripts/%d%s.csh", ibase, af[ibase].fname[jbase].name);
	  //system(dump);
	  //fprintf(stderr,"outdir : %s\n",outdir);
	  sprintf(cvobsac,"%s%s",outdir, af[ibase].fname[jbase].name);
          fprintf(stderr,"%s\n",af[ibase].fname[jbase].name);
          fprintf(stderr,"%s\n",cvobsac);
//          abort();
	  cvSAC_HEADER.delta=1.0;		//????why the delta is so important?
          write_sac (cvobsac,&(cvsig[0]),&cvSAC_HEADER);
        }
    } /* end if loop */

    /* reset the indices */
    quit = 0;
    jnotbs = 0;
    jbase++;  //move to the next CORsac in af[ibase]
    } while(jbase != af[ibase].len); /* outer do while */
    jbase = 0;
  } /* end for loop */

  return 0;
}

