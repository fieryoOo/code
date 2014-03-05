#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>

extern"C" {
   void rad_pattern_r_(char *eig_fname, int *eig_namelen, char *phv_fname, int *phv_namelen, 
	               float *strike, float *dip, float *rake, float *depth, float *per, int *nper, 
		       float *azi, float grT[][181], float phT[][181], float amp[][181]);
}

int main() {
   float per[3] = {10, 20, 30};
   int nper = sizeof(per)/sizeof(per[0]);
   float azi[181];
   float grT[nper][181], phT[nper][181], amp[nper][181];
   char eig_fname[] = "245.25_41.R", phv_fname[] = "245.25_41.R.phv";
   int eig_namelen=strlen(eig_fname), phv_namelen=strlen(phv_fname);
   float strike = 192, dip = 34, rake = -101, depth = 10;
   rad_pattern_r_( eig_fname, &eig_namelen, phv_fname, &phv_namelen,
		   &strike, &dip, &rake, &depth, per, &nper, azi, grT, phT, amp);
   for(int iper=0; iper<nper; iper++) for(int i=0; i<181; i++) std::cerr<<azi[i]<<" "<<phT[iper][i]<<std::endl;

   return 0;
}
