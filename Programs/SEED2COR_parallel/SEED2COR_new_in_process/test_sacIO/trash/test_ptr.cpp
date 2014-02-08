#include <cstdio>
#include <iostream>
#include <algorithm>
#include <memory>
#include <cstring>

int main() {
   std::unique_ptr<float[]> sig0(new float[10]), sig1(new float[5]);
   std::fill(&(sig0[0]), &(sig0[0])+10, 1.e30);
   std::fill(&(sig1[0]), &(sig1[0])+5, -1.e1);
   std::copy(&(sig1[0]), &(sig1[0])+5, &(sig0[0]));
   
   for(int i=0;i<10;i++) std::cout<<sig0[i]<<std::endl;
   return 0;
}
