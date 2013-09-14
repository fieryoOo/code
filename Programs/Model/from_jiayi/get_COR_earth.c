#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* The gen_cor_pred_files.c program generates the 
   COR_per1.per2_STA1_STA2.SAC_PRED files needed for the gmt plotting
   scripts for dispersion curves. The program reads the PREDICTION_R 
   file and generates files for all station pairs withi the *_R file */

/* FUNCTION PROTOTYPES */
void trailing_slash(char *string_in, char *string_out);


/*------------------------------------------------------------------------*/
void trailing_slash(char *string_in, char *string_out)
/*------------------------------------------------------------------------*/
/* the function checks the last character of the string, string_in, and adds a forward slash to the end if it is not already there. Used for appending directory names that will be used in an absolute path call */
{
  int len_str_in;

  len_str_in=strlen(string_in);
  if (string_in[len_str_in-1] == '/') {
    strcpy(string_out, string_in);
  }
  else {
    strcpy(string_out, string_in);
    string_out[len_str_in] = '/';
    string_out[len_str_in+1] = '\0';
  }
}


/*--------------------------------------------------------------------------*/
int main(int argc, char *argv[])
/*--------------------------------------------------------------------------*/
{
  FILE *fp1, *fp2;  /* fp1-PREDICTION_? file, fp2-output COR_*.SAC_PRED files */
  int cnt1, len_path;
  int no_sta1, no_sta2, no_gpvel;
  float lon1, lat1, lon2, lat2;
  float per, gpvel;
  char ch1;
  char buff[200], sta1[18], sta2[18];
  char filename[100], path_out[100], fullfilename[200];

/* CHECK FOR COMMAND LINE ARG */
  if(argc!=3) {
    printf("USAGE: gen_cor_pred_files [PREDICTION_R file] [destination dir] >& err_log_cor_pred\n");
    exit(1);
  }

/* requires that the path have a trailing "/" */
  trailing_slash(argv[2],path_out);

/* OPENING FILES */
  if((fp1 = fopen(argv[1], "r"))==NULL) {
    fprintf(stderr,"cannot open prediction file: %s\n", argv[1]);
    exit(1);
  }

/* READ ALL STATION PAIRS AND GROUP VELOCITIES TO _PRED files */
  while ( !feof(fp1) ){
    fgets(buff, 180, fp1);
    sscanf(buff,"%d %d %d %s %s %f %f %f %f", &no_sta1, &no_sta2, &no_gpvel, sta1, sta2, &lat1, &lon1, &lat2, &lon2);

    sprintf(filename,"%s.%s.dat\0", sta1,sta2);
    strcpy(fullfilename,path_out);
    strcat(fullfilename,filename);
    if((fp2 = fopen(fullfilename, "w"))==NULL) {
      fprintf(stderr,"cannot open _PRED file: %s\n", fullfilename);
      exit(1);
    }

 //   fprintf(fp2,"  %d %d %d %s %s %f %f %f %f\n", no_sta1, no_sta2, no_gpvel, sta1, sta2, lat1, lon1, lat2, lon2);

    for(cnt1=0; cnt1<no_gpvel; cnt1++) {
      fgets(buff, 180, fp1);
      sscanf(buff,"%f %f", &per, &gpvel);
      fprintf(fp2,"    %5.1f	%6f\n", per, gpvel);
    }
  fclose(fp2);
  ch1=getc(fp1);
  }

  fclose(fp1);
  return 0;
}


