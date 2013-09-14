#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
using namespace std;

#define BLCS 5000

struct MapData {
   float lon, lat;
   double key;
   struct MapData *next;
};

int nint(float a) {
   return (int)floor(a+0.5);
}

main(int argc, char *argv[]) 
{
   if (argc!=5) {
      cout<<"Usage: "<<argv[0]<<" [Point.loc (lon lat)] [Input_Map (lon lat value)] [avg_dist] [Output_file]"<<endl;
      return 0;
   }

   struct MapData *mdata=NULL, Point;
   FILE *fin, *fout;
   int i, nblock, ndat;
   char buff[300];
   float lonmin, lonmax, latmin, latmax;
//read in map
   if( (fin = (fopen(argv[2],"r"))) == NULL ) {
      cout<<"Cannot access file "<<argv[2]<<endl;
      exit(0);
   }
   for(nblock=0,i=0;fgets(buff, 300, fin)!=NULL;) {
      if(nblock*BLCS<=i)
         mdata = (struct MapData *) realloc (mdata, (++nblock)*BLCS * sizeof(struct MapData));
      if(sscanf(buff, "%f %f %lf", &(mdata[i].lon), &(mdata[i].lat), &(mdata[i].key))!=3) {
	 cout<<"Warning: format error: "<<buff<<endl;
	 continue;
      }
      if(i==0) { 
	 lonmin = mdata[i].lon; lonmax = mdata[i].lon;
	 latmin = mdata[i].lat; latmax = mdata[i].lat;
      }
      else {
	 if( lonmin > mdata[i].lon ) lonmin = mdata[i].lon;
	 else if( lonmax < mdata[i].lon ) lonmax = mdata[i].lon;
	 if( latmin > mdata[i].lat ) latmin = mdata[i].lat;
	 else if( latmax < mdata[i].lat ) latmax = mdata[i].lat;
      }
      i++;
   }
   ndat = i;
   fclose(fin);

//produce Hash table
   latmin = nint(latmin);
   int nlat = nint(latmax)-(int)latmin+1;
   lonmin = nint(lonmin);
   int nlon = nint(lonmax)-(int)lonmin+1;
   int ilon, ilat, j;
   struct MapData *Hash[nlon][nlat];
   for(i=0;i<nlon;i++) for(j=0;j<nlat;j++) Hash[i][j] = NULL;
   for(i=0;i<ndat;i++) {
      ilon = nint(mdata[i].lon);
      ilat = nint(mdata[i].lat);
      mdata[i].next = Hash[ilon][ilat];
      Hash[ilon][ilat] = &(mdata[i]);
   }

//read in and process points
   if( (fin = (fopen(argv[1],"r"))) == NULL ) {
      cout<<"Cannot access file "<<argv[2]<<endl;
      exit(0);
   }
   fout = fopen(argv[4],"w");
   for(i=0;fgets(buff, 300, fin)!=NULL;) {
      if(sscanf(buff, "%f %f %lf", &(Point.lon), &(Point.lat))!=2) {
         cout<<"Warning: format error: "<<buff<<endl;
         continue;
      }
      {
	 ilon = nint(Point.lon);
	 ilat = nint(Point.lat);
         Point.key = Point.lon Point.lat mdata[i].lon mdata[i].lat mdata[i].key
      }
      fprintf(fout, "%f %f %lf\n", Point.lon, Point.lat, Point.key);
      i++;
   }
   fclose(fin); fclose(fout);
   cout<<i<<" points processed."<<endl;
   
}
