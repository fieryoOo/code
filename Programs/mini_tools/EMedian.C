#include <stdio.h>
#include <math.h>
#include <iostream>
#include <algorithm>
using namespace std;

float Median4( float *arr, int npt ) {
   if(npt == 1) return arr[0];
   if(npt == 2) return (arr[0]+arr[1]) / 2.;
   float a, b, c, d;
   if( arr[0]<arr[1] ) { a=arr[0]; b=arr[1]; }
   else { a=arr[1]; b=arr[0]; }
   c = arr[2];
   if(npt == 3) {
      if( c<=a ) return a;
      if( c>=b ) return b;
      return c;
   }
   if( c<arr[3] ) { d=arr[3]; }
   else { d=c; c=arr[3]; }
   if( a < c ) return min(b,c);
   else return min(a,d);
}
float Median5( float *arr ) {
   float a, b, c, d;

   if( arr[0]<arr[1] ) { a=arr[0]; b=arr[1]; }
   else { a=arr[1]; b=arr[0]; } //sort a b
   if( arr[2]<arr[3] ) { c=arr[2]; d=arr[3]; }
   else { c=arr[3]; d=arr[2]; } //sort c d
   // eliminate the smallest of 4 and get arr[4] in
   if( a < c ) {
      if( arr[4] < b ) { a = arr[4]; }
      else { a = b; b = arr[4]; }
   }
   else {
      if( arr[4] < d ) { c = arr[4]; }
      else { c = d; d = arr[4]; }
   }
   // eliminate another smallest of 4 and median is the 3rd smallest
   if( a < c ) return min(b,c);
   else return min(a,d);
}
float MedianOfMedian( float *arr, int npt ) {
   if( npt < 5 ) return Median4(arr, npt);
   if( npt == 5) return Median5(arr);
   int i, nstep = npt/5;
   float arr2[nstep+1];
   for(i=0; i<nstep; i++) arr2[i] = Median5( &(arr[i*5]) );
   if(npt==nstep*5) return MedianOfMedian( arr2, nstep );
   arr2[nstep] = Median5( &(arr[npt-5]) );
   return MedianOfMedian( arr2, nstep+1 );
}

int main() {
   float arr[20];
   arr[0]=3; arr[1]=5; arr[2]=1; arr[3]=2; arr[4]=4;
   arr[5]=10; arr[6]=7; arr[7]=8; arr[8]=9; arr[9]=6;
   arr[10]=15; arr[11]=12; arr[12]=11; arr[13]=13; arr[14]=14;
   arr[15]=18; arr[16]=19; arr[17]=20; arr[18]=17; arr[19]=15;
   cout<<MedianOfMedian( &arr[0], 16 )<<endl;
}
