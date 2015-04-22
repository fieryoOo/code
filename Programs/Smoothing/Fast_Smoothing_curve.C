#define MAIN
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
using namespace std;

#define NMOD 5000

int Smoothing (char *infile, double hwidth) {
	FILE *ff;
	int i, j, ndata, ism, nmod=0;
	char buff[300];
	double *datax = NULL, *datay = NULL, xtmp, ytmp, midx;
	double dx=hwidth,  alpha=0.5/(hwidth*hwidth);

	if((ff=fopen(infile,"r"))==NULL) {
		cout<<"Can't open file "<<infile<<endl;
		return 0;
	}

	for(ndata=0;;ndata++) {
		if( nmod*NMOD <= ndata ) {
			datax = (double *) realloc (datax, (++nmod)*NMOD * sizeof(double));
			datay = (double *) realloc (datay, nmod*NMOD * sizeof(double));
		}
		if((fgets(buff, 300, ff))==NULL) break;
		sscanf(buff,"%lf %lf", &xtmp, &ytmp);
		for(i=ndata-1;i>=0&&xtmp<datax[i];i--) { 
			datax[i+1] = datax[i];
			datay[i+1] = datay[i];
		}
		datax[i+1] = xtmp; datay[i+1] = ytmp;
	}
	fclose(ff);

	//double hwidth2 = hwidth * 2;
	i = (int)floor((datax[0]+hwidth)/dx); // min midx from datax[0]
	sprintf(buff,"%s_sm\0",infile);
	ff=fopen(buff,"w");
	// loop through each avg point: midx=i*dx
	int jstart = 0, jend = 0;
	for(j=0; ;i++) {
		midx = i*dx;
		if(midx>datax[ndata-1]-hwidth) break;		// max midx from datax[last]
		for(j=jstart;datax[j]<midx-hwidth;j++);		// starting j for current midx
		jstart = j;					// next jstart cannot be smaller than current starting j
		for(j=jend; j<ndata; j++) { if( datax[j]>midx+hwidth) break; }
		jend = j;						// next jend cannot be smaller than current starting j
		//std::cerr<<jstart<<" "<<jend<<" "<<datax[jstart]<<" "<<datax[jend]<<" "<<midx<<"\n";
		if( jend - jstart == 0 ) continue;
		double mean = 0., V1 = 0.,  weight[jend-jstart];
		for(j=jstart; j<jend ;j++) {			// loop through data in the window for mean
			//weight[j-jstart]=exp(-alpha*pow((datax[j]-midx),2));
			weight[j-jstart] = 1.;
			mean += datay[j] * weight[j-jstart];
			V1 += weight[j-jstart];
		}
		mean /= V1;
		double var = 0., V2 = 0.;
		for(j=jstart; j<jend ;j++) {			// loop through data 2nd time for standard deviation
			double dtmp = datay[j] - mean;
			var += weight[j-jstart] * dtmp * dtmp;
			V2 += weight[j-jstart] * weight[j-jstart];
		}
		var *= V1/(V1*V1-V2);
		fprintf(ff, "%lf %lf %lf %d\n", midx, mean, sqrt(var), jend-jstart );
	}
	fclose(ff);
	free(datax); free(datay);

	return 1;
}

int main (int argc, char *argv[])
{ 
	if (argc != 3) {
		cout<<argv[0]<<" [input_file] [half_width]"<<endl;
		return 0;
	}

	Smoothing(argv[1], atof(argv[2]));

	return 1;
}
