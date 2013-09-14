#include <stdio.h>
#include <string.h>
#include <iostream>
#include <unistd.h>
#include <stdlib.h>
using namespace std;

struct nd {
   double data;
   struct nd ** pnext;
};
typedef struct nd node ;

int main () {
   node *p;   
}
