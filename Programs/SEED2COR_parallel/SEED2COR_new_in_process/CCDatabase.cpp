#include "CCDatabase.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cstring>
#include <cmath>

static int isTermi(int deno) {
   while(deno%2==0) deno /= 2;
   while(deno%5==0) deno /= 5;
   return deno==1;
}
/* Constructor: read in parameters from the inputfile */
CCDatabase::CCDatabase( char *fname ) {
   FILE *fparam, *filetmp;
   char buff[300], ctmp[200];
   int ERR=0, itmp;
   float ftmp;
   if((fparam=fopen(fname,"r"))==NULL) {
      std::cerr<<"Error(GetParam): Cannot open parameter file "<<fname<<std::endl;
      exit(0);
   }
   std::cout<<"---------------------------Checking input parameters---------------------------"<<std::endl;
   //CCParams.rdsexe
   if( (fgets(buff, 300, fparam)) == NULL ) {std::cerr<<"   Error: No enough parameters read in from the input file!"<<std::endl; exit(0);}
   sscanf(buff, "%s", CCParams.rdsexe);
   std::cout<<"rdseed excutable:\t"<<CCParams.rdsexe<<std::endl;
   if( ! strstr(CCParams.rdsexe, "rdseed") ) std::cout<<"   Warning: Are you sure this is an rdseed excutable?"<<std::endl;
   if( access(CCParams.rdsexe, R_OK)!=0 ) {
     std::cerr<<"   Error: cannot access rdseed through "<<CCParams.rdsexe<<std::endl;
     ERR = 1;
   }
   //CCParams.evrexe
   if( (fgets(buff, 300, fparam)) == NULL ) {std::cerr<<"   Error: No enough parameters read in from the input file!"<<std::endl; exit(0);}
   sscanf(buff, "%s", CCParams.evrexe);
   std::cout<<"evalresp excutable:\t"<<CCParams.evrexe<<std::endl;
   if( ! strstr(CCParams.evrexe, "evalresp") ) std::cout<<"   Warning: Are you sure this is an evalresp excutable?"<<std::endl;
   if( access(CCParams.evrexe, R_OK)!=0 ) {
     std::cerr<<"   Error: cannot access evalresp through "<<CCParams.evrexe<<std::endl;
     ERR = 1;
   }
   //CCParams.stalst
   if( (fgets(buff, 300, fparam)) == NULL ) {std::cerr<<"   Error: No enough parameters read in from the input file!"<<std::endl; exit(0);}
   sscanf(buff, "%s", CCParams.stalst);
   std::cout<<"station list:\t\t"<<CCParams.stalst<<std::endl;
   if((filetmp=fopen(CCParams.stalst,"r"))==NULL) {
      std::cerr<<"   Error: cannot access file "<<CCParams.stalst<<std::endl;
      ERR = 1;
   }
   else { 
      fgets(buff, 300, filetmp);
      if(sscanf(buff,"%s %f %f", ctmp, &ftmp, &ftmp)!=3) {
         std::cerr<<"   Error: incorrect format in file "<<CCParams.stalst<<"! Shoud be (sta lon lat)"<<std::endl;
         ERR = 1;
      }
      fclose(filetmp);
   }
   //CCParams.seedlst
   if( (fgets(buff, 300, fparam)) == NULL ) {std::cerr<<"   Error: No enough parameters read in from the input file!"<<std::endl; exit(0);}
   sscanf(buff, "%s", CCParams.seedlst);
   std::cout<<"seed list:\t\t"<<CCParams.seedlst<<std::endl;
   if((filetmp=fopen(CCParams.seedlst,"r"))==NULL) {
      std::cerr<<"   Error: cannot access file "<<CCParams.seedlst<<std::endl;
      ERR = 1;
   }
   else {
      fgets(buff, 300, filetmp);
      if(sscanf(buff,"%s %d %d %d", ctmp, &itmp, &itmp, &itmp)!=4) {
         std::cerr<<"   Error: incorrect format in file "<<CCParams.seedlst<<"! (path/seedfile yyyy mm dd) is expected"<<std::endl;
         ERR = 1;
      }
      fclose(filetmp);
   }
   //CCParams.ch
   if( (fgets(buff, 300, fparam)) == NULL ) {std::cerr<<"   Error: No enough parameters read in from the input file!"<<std::endl; exit(0);}
   sscanf(buff, "%s", CCParams.ch);
   std::cout<<"channel:\t\t"<<CCParams.ch<<std::endl;
   //CCParams.sps
   if( (fgets(buff, 300, fparam)) == NULL ) {std::cerr<<"   Error: No enough parameters read in from the input file!"<<std::endl; exit(0);}
   sscanf(buff, "%f", &ftmp);
   CCParams.sps = (int)ftmp;
   std::cout<<"target sampling rate:\t"<<ftmp<<std::endl;
   if( CCParams.sps!=ftmp || CCParams.sps <= 0 ) {
      std::cerr<<"   Error: a positive integer is expected!"<<std::endl;
      ERR = 1;
   }
   else if( !isTermi(CCParams.sps) ) {
      std::cerr<<"   Error: 1/"<<CCParams.sps<<" isn't a terminating decimal, will cause rounding error!"<<std::endl;
      ERR = 1;
   }
   //CCParams.gapfrac
   if( (fgets(buff, 300, fparam)) == NULL ) {std::cerr<<"   Error: No enough parameters read in from the input file!"<<std::endl; exit(0);}
   sscanf(buff, "%f", &CCParams.gapfrac);
   std::cout<<"max gap fraction:\t"<<CCParams.gapfrac*100<<"%"<<std::endl;
   if( CCParams.gapfrac<0 || CCParams.gapfrac>1 ) {
      std::cerr<<"   Error: a number between 0. - 1. is expected!"<<std::endl;
      ERR = 1;
   }
   //CCParams.t1
   if( (fgets(buff, 300, fparam)) == NULL ) {std::cerr<<"   Error: No enough parameters read in from the input file!"<<std::endl; exit(0);}
   sscanf(buff, "%f", &CCParams.t1);
   std::cout<<"cutting begining:\t"<<CCParams.t1<<"sec"<<std::endl;
   if( fabs(CCParams.t1) > 86400 ) std::cout<<"   Warning: "<<CCParams.t1<<"sec exceeds one day."<<std::endl;
   //CCParams.tlen
   if( (fgets(buff, 300, fparam)) == NULL ) {std::cerr<<"   Error: No enough parameters read in from the input file!"<<std::endl; exit(0);}
   sscanf(buff, "%f", &CCParams.tlen);
   std::cout<<"time-rec length:\t"<<CCParams.tlen<<"sec"<<std::endl;
   ftmp = CCParams.t1+CCParams.tlen;
   if( ftmp>86400 || ftmp<0 ) std::cout<<"   Warning: ending time '"<<ftmp<<"sec' out of range."<<std::endl;
   //CCParams.perl CCParams.perh
   if( (fgets(buff, 300, fparam)) == NULL ) {std::cerr<<"   Error: No enough parameters read in from the input file!"<<std::endl; exit(0);}
   sscanf(buff, "%f", &CCParams.perl);
   if( (fgets(buff, 300, fparam)) == NULL ) {std::cerr<<"   Error: No enough parameters read in from the input file!"<<std::endl; exit(0);}
   sscanf(buff, "%f", &CCParams.perh);
   std::cout<<"signal per band:\t"<<CCParams.perl<<" - "<<CCParams.perh<<"sec"<<std::endl;
   if( CCParams.perl<0 || CCParams.perh<0 ) {
      std::cerr<<"   Error: period band can not go below 0."<<std::endl;
      ERR = 1;
   }
   if(CCParams.perl<2./CCParams.sps) {
      std::cerr<<"   Error: "<<CCParams.perl<<"sec is out of the lower limit at sps "<<CCParams.sps<<std::endl;
      ERR = 1;
   }
   if(CCParams.perl<5./CCParams.sps) std::cout<<"   Warning: signal at "<<CCParams.perl<<"sec might be lost at a sampling rate of "<<CCParams.sps<<std::endl;
   //CCParams.tnorm_flag
   if( (fgets(buff, 300, fparam)) == NULL ) {std::cerr<<"   Error: No enough parameters read in from the input file!"<<std::endl; exit(0);}
   sscanf(buff, "%d", &CCParams.tnorm_flag);
   std::cout<<"t-norm method:\t\t";
   if(CCParams.tnorm_flag==0) std::cout<<"none"<<std::endl;
   else if(CCParams.tnorm_flag==1) std::cout<<"One-bit"<<std::endl;
   else if(CCParams.tnorm_flag==2) std::cout<<"Running average"<<std::endl;
   else if(CCParams.tnorm_flag==3) std::cout<<"Earthquake cutting"<<std::endl;
   else {
      std::cerr<<std::endl<<"   Error: Unknow method. integer between 0 - 3 is expected"<<std::endl;
      ERR = 1;
   }
   //CCParams.Eperl CCParams.Eperh
   if( (fgets(buff, 300, fparam)) == NULL ) {std::cerr<<"   Error: No enough parameters read in from the input file!"<<std::endl; exit(0);}
   sscanf(buff, "%f", &CCParams.Eperl);
   if( (fgets(buff, 300, fparam)) == NULL ) {std::cerr<<"   Error: No enough parameters read in from the input file!"<<std::endl; exit(0);}
   sscanf(buff, "%f", &CCParams.Eperh);
   if(CCParams.tnorm_flag!=1) {
      if( CCParams.Eperl == -1. ) std::cout<<"Eqk filter:\t\toff"<<std::endl;
      else {
         std::cout<<"Eqk filter:\t\t"<<CCParams.Eperl<<" - "<<CCParams.Eperh<<"sec"<<std::endl;
         if( CCParams.Eperl<0 || CCParams.Eperh<0 ) {
            std::cerr<<"   Error: period band can not go below 0."<<std::endl;
            ERR = 1;
         }
      }
   }
   //CCParams.timehlen
   if( (fgets(buff, 300, fparam)) == NULL ) {std::cerr<<"   Error: No enough parameters read in from the input file!"<<std::endl; exit(0);}
   sscanf(buff, "%f", &CCParams.timehlen);
   if(CCParams.tnorm_flag!=1) {
      std::cout<<"t-len for run-avg:\t"<<CCParams.timehlen<<std::endl;
      if(CCParams.timehlen<0) {
         std::cerr<<"   Error: positive number is expected"<<std::endl;
         ERR = 1;
      }
   }
   //CCParams.frechlen
   if( (fgets(buff, 300, fparam)) == NULL ) {std::cerr<<"   Error: No enough parameters read in from the input file!"<<std::endl; exit(0);}
   sscanf(buff, "%f", &CCParams.frechlen);
   std::cout<<"t-len for whitening:\t"<<CCParams.frechlen<<"  ";
   if(CCParams.frechlen==-1) std::cout<<"input smoothing file will be used";
   else if(CCParams.frechlen<0) {
      std::cerr<<std::endl<<"   Error: non-negative number is expected";
      ERR = 1;
   }
   std::cout<<std::endl;
   //CCParams.fwname
   if( (fgets(buff, 300, fparam)) == NULL ) {std::cerr<<"   Error: No enough parameters read in from the input file!"<<std::endl; exit(0);}
   sscanf(buff, "%s", CCParams.fwname);
   if(CCParams.frechlen==-1) {
      std::cout<<"spec reshaping file:\t"<<CCParams.fwname<<std::endl;
      if( access(CCParams.fwname, R_OK)!=0 ) {
         std::cerr<<"   Error: cannot access file "<<CCParams.fwname<<std::endl;
         ERR = 1;
      }  
   }
   //CCParams.ftlen
   if( (fgets(buff, 300, fparam)) == NULL ) {std::cerr<<"   Error: No enough parameters read in from the input file!"<<std::endl; exit(0);}
   sscanf(buff, "%d", &CCParams.ftlen);
   if(CCParams.ftlen<=0) std::cout<<"cor-time-len correction\toff"<<std::endl;
   else {
      std::cout<<"cor-time-len correction\ton"<<std::endl;
      if(CCParams.tnorm_flag==3) CCParams.ftlen = 2;
   }
   //CCParams.fprcs
   if( (fgets(buff, 300, fparam)) == NULL ) {std::cerr<<"   Error: No enough parameters read in from the input file!"<<std::endl; exit(0);}
   sscanf(buff, "%d", &CCParams.fprcs);
   if(CCParams.fprcs<=0) std::cout<<"prcsr-signal checking\toff"<<std::endl;
   else std::cout<<"prcsr-signal checking\ton"<<std::endl;
   //CCParams.memomax
   if( (fgets(buff, 300, fparam)) == NULL ) {std::cerr<<"   Error: No enough parameters read in from the input file!"<<std::endl; exit(0);}
   sscanf(buff, "%f", &CCParams.memomax);
   std::cout<<"memory fraction:\t"<<CCParams.memomax*100<<"%"<<std::endl;
   if( CCParams.memomax<0 || CCParams.memomax>1 ) {
      std::cerr<<"   Error: a number between 0. - 1. is expected!"<<std::endl;
      ERR = 1;
   }
   //CCParams.lagtime
   if( (fgets(buff, 300, fparam)) == NULL ) {std::cerr<<"   Error: No enough parameters read in from the input file!"<<std::endl; exit(0);}
   sscanf(buff, "%d", &CCParams.lagtime);
   std::cout<<"cor lag time:\t\t"<<CCParams.lagtime<<"sec"<<std::endl;
   if(CCParams.lagtime<0) {
      std::cerr<<"   Error: negative CCParams.lagtime!"<<std::endl;
      ERR = 1;
   }
   else if(CCParams.lagtime>524288) std::cout<<"   Warning: lag time exceeds the maximum"<<std::endl;
   //CCParams.mintlen
   if( (fgets(buff, 300, fparam)) == NULL ) {std::cerr<<"   Error: No enough parameters read in from the input file!"<<std::endl; exit(0);}
   sscanf(buff, "%d", &CCParams.mintlen);
   std::cout<<"min time len:\t\t"<<CCParams.mintlen<<"sec"<<std::endl;
   if( CCParams.mintlen>86400 ) std::cout<<"   Warning: allowed minimum time length larger than a day."<<std::endl;
   //CCParams.fdelosac
   if( (fgets(buff, 300, fparam)) == NULL ) {std::cerr<<"   Error: No enough parameters read in from the input file!"<<std::endl; exit(0);}
   sscanf(buff, "%d", &CCParams.fdelosac);
   std::cout<<"Delete orig sacs?\t";
   if(CCParams.fdelosac==1) std::cout<<"Yes"<<std::endl;
   else std::cout<<"No"<<std::endl;
   //CCParams.fdelamph
   if( (fgets(buff, 300, fparam)) == NULL ) {std::cerr<<"   Error: No enough parameters read in from the input file!"<<std::endl; exit(0);}
   sscanf(buff, "%d", &CCParams.fdelamph);
   std::cout<<"Delete am&ph files?\t";
   if(CCParams.fdelamph==1) std::cout<<"Yes"<<std::endl;
   else std::cout<<"No"<<std::endl;
   //CCParams.fskipesac
   if( (fgets(buff, 300, fparam)) == NULL ) {std::cerr<<"   Error: No enough parameters read in from the input file!"<<std::endl; exit(0);}
   sscanf(buff, "%d", &CCParams.fskipesac);
   std::cout<<"ExtractSac()\t\t";
   if(CCParams.fskipesac==2) std::cout<<"skip"<<std::endl;
   else if(CCParams.fskipesac==1) std::cout<<"skip if target file exists"<<std::endl;
   else std::cout<<"overwrite if target file exists"<<std::endl;
   //CCParams.fskipresp
   if( (fgets(buff, 300, fparam)) == NULL ) {std::cerr<<"   Error: No enough parameters read in from the input file!"<<std::endl; exit(0);}
   sscanf(buff, "%d", &CCParams.fskipresp);
   std::cout<<"RmRESP()\t\t";
   if(CCParams.fskipresp==2) std::cout<<"skip"<<std::endl;
   else if(CCParams.fskipresp==1) std::cout<<"skip if target file exists"<<std::endl;
   else std::cout<<"overwrite if target file exists"<<std::endl;
   //CCParams.fskipamph
   if( (fgets(buff, 300, fparam)) == NULL ) {std::cerr<<"   Error: No enough parameters read in from the input file!"<<std::endl; exit(0);} 
   sscanf(buff, "%d", &CCParams.fskipamph);
   std::cout<<"TempSpecNorm()\t\t";
   if(CCParams.fskipamph==2) std::cout<<"skip"<<std::endl;
   else if(CCParams.fskipamph==1) std::cout<<"skip if target file exists"<<std::endl;
   else std::cout<<"overwrite if target file exists"<<std::endl;
   //CCParams.fskipcrco
   if( (fgets(buff, 300, fparam)) == NULL ) {std::cerr<<"   Error: No enough parameters read in from the input file!"<<std::endl; exit(0);}
   sscanf(buff, "%d", &CCParams.fskipcrco);
   std::cout<<"CrossCorr()\t\t";
   if(CCParams.fskipcrco==2) std::cout<<"skip"<<std::endl;
   else if(CCParams.fskipcrco==1) std::cout<<"skip if target file exists"<<std::endl;
   else std::cout<<"overwrite if target file exists"<<std::endl;
   //CCParams.CorOutflag
   if( (fgets(buff, 300, fparam)) == NULL ) {std::cerr<<"   Error: No enough parameters read in from the input file!"<<std::endl; exit(0);}
   sscanf(buff, "%d", &CCParams.CorOutflag);
   std::cout<<"CC Output: \t\t";
   switch( CCParams.CorOutflag ) {
      case 0:
	 std::cout<<"monthly"<<std::endl;
	 break;
      case 1:
	 std::cout<<"daily"<<std::endl;
	 break;
      case 2:
	 std::cout<<"monthly + daily"<<std::endl;
	 break;
      default:
	 std::cerr<<std::endl<<"   Error: unknow outflag("<<CCParams.CorOutflag<<")!"<<std::endl;
	 ERR = 1;
   }
   //check for EOF
   if( (fgets(buff, 300, fparam)) != NULL ) std::cout<<"   Warning: End of file not reached!"<<std::endl;

   std::cout<<"-------------------------------Checking completed-------------------------------"<<std::endl;
   fclose(fparam);
   if(ERR==1) exit(0);
   char cin;
   std::cout<<"Continue?  ";
   fgets(buff, 300, stdin);
   sscanf(buff, "%c", &cin);
   if( cin!='Y' && cin!='y' ) exit(0);

}

