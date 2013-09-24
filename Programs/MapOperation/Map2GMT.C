/* Read in from a file one value at a time assuming each is a data point on a map with lon/lat increased constantly
 starting at (lon0, lat0). Lon goes first and will step nlon times before it falls back to lon0 and carry/borrow one
step to/from lat.
*/

#include <stdio.h>
#include <iostream>
using namespace std;

int main(int argc, char *argv[]) {
   if( argc != 8 ) {
      cerr<<"Usage: "<<argv[0]<<" [input_map_file] [lon0] [lat0] [lon_step] [lat_step] [nlon] [nlat]"<<endl;
      exit(-1);
   }
   // read in data
   int i, ndat;
   int nlon = atoi(argv[6]), nlat = atoi(argv[7]), slon;
   float *data = new float[nlon*nlat+1000];
   FILE *fin, *fout;
   if( (fin=fopen(argv[1],"r")) == NULL ) {
      cerr<<"Error(main): Cannot access file "<<argv[1]<<endl;
      exit(0);
   }
   for(ndat=0;fscanf(fin, "%f", &data[ndat])==1;ndat++) {}
   fclose(fin);
   // check for number of data points
   if( ndat != nlon*nlat ) {
      char cin, buff[100];
      cerr<<"Warning(main): input data number("<<ndat<<") doesn't match nlon*nlat("<<nlon*nlat<<")! Continue? ";
      fgets(buff, 300, stdin);
      sscanf(buff, "%c", &cin);
      if( cin!='Y' && cin!='y' ) exit(0);
   }
   // assign location and output
   char outname[150];
   sprintf(outname, "%s.gmt", argv[1]);
   if( (fout=fopen(outname,"w")) == NULL ) {
      cerr<<"Error(main): Cannot open file "<<outname<<endl;
      exit(0);
   }
   float lon0 = atof(argv[2]), lat0 = atof(argv[3]), lonstep = atof(argv[4]), latstep = atof(argv[5]);
   float lon=lon0, lat=lat0;
   for( slon=0,i=0; i<ndat; i++,slon++,lon+=lonstep ) {
      if(slon == nlon) { slon = 0; lon = lon0; lat += latstep; }
      fprintf(fout, "%f %f %f\n", lon, lat, data[i]);
   }
   fclose(fout);
   delete [] data;

   return 0;
}
