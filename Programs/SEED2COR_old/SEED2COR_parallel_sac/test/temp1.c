#include <stdio.h>
#include <stdlib.h>
//#include <string.h>
//#include <iostream>
#include <errno.h>
//using namespace std;

int main() {
   FILE *fin = NULL;
   if( (fin=fopen("asda","r")) == NULL ) {
      perror("fopen");
      exit(0);
   }

   return 1;
}
