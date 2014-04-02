#include <cstdio>
#include <cstdlib>
#include <iostream>

int main() {
   const int size = 16*1024*1024;
   char *c = (char *) malloc (size);
   for( int i=0; i<0xffff; i++ ) {
      for( int j=0; j<size; j++ ) {
	 c[j] = j;
      }
   }
   return 0;
}
