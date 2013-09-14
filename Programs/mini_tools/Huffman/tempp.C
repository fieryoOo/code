#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <unistd.h>
#include <stdlib.h>
#include <limits.h>
#include "boost/random.hpp"
#include "boost/generator_iterator.hpp"
#include <sys/time.h>
using namespace std;

void gotoxy(int x, int y)
{
    printf("%c[%d;%df",0x1B,y,x);
}

int main()
{
       gotoxy(200,5);
       cout<<"wawa";
       gotoxy(5,15);
       cout<<"hha";
       gotoxy(100,2);
       cout<<"xixi";

gotoxy(1,15);
cout<<endl;
       return 0;
}

