#include<stdio.h>
#include <stdlib.h>
#include<iostream>
#include<iomanip>
#include<math.h>
#include <string>
#include <unistd.h>
using namespace std;

#define max(a,b) (((a) > (b)) ? (a) : (b))
int main()
{
   double a=max(13.,18);
   cout<<a<<endl;
   return 1;
}
