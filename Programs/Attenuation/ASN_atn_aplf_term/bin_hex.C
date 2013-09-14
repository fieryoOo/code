#include <math.h>

void bin2hex( short* bin, char* hex )
{
   short i, j, hv;
   for(i=0;i<8;i++){
      hv = 0;
      for(j=0;j<4;j++){
         hv += bin[i*4+3-j]*(int)pow(2,j);
      }
      if( hv >= 0 && hv <=9 ) hex[i] = hv + '0';
      else hex[i] = hv - 10 + 'A';
   }
   hex[i]='\0';
}

void hex2bin( char* hex, short* bin )
{
   short i, j, hv, len = 8;
   for(i=0;i<len;i++){
      if (hex[i] >= '0' && hex[i] <= '9') hv = hex[i] - '0';
      else hv = hex[i] - 'A' + 10;
      for(j=i*4+3;j>i*4;j--){
         bin[j]= hv % 2;
         hv = hv/2;
      }
      bin[i*4]=hv;
   }
}

