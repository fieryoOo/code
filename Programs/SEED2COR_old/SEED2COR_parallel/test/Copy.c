#include <stdio.h>
#define BSZ 8192
void Copy(char *oldname, char *newname) {
   FILE *fin, *fou;
   if( (fin = fopen(oldname, "r")) == NULL ) {
      perror("### Warning from Copy (fopen)");
      return;
   }
   if( (fou = fopen(newname, "w")) == NULL ) {
      perror("### Warning from Copy (fopen)");
      fclose(fin);
      return;
   }
   int result;
   char buff[BSZ];
   while ( result = fread(buff, 1, BSZ, fin) ) {
       fwrite(buff, 1, result, fou);
   }
   fclose(fin); fclose(fou);
}

int main() {
   Copy("ft_2012.FEB.11.G03A.BHZ.SAC", "temp3");
   return 1;
}
