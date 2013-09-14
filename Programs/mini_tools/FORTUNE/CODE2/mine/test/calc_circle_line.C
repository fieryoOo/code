#include <math.h>
#include <iostream.h>

int main(int argc, char ** argv)
{
   if( argc != 7 ) {
      cout<<"Usage: "<<argv[0]<<" [al] [bl] [cl] [x1] [y1] [y3]"<<endl;
      exit(-1);
   }

   double al = atof(argv[1]), bl = atof(argv[2]), cl = atof(argv[3]);
   double x1 = atof(argv[4]), y1 = atof(argv[5]), y3 = atof(argv[6]);
   double a = bl*bl;
   double b = 2*x1*al*bl+2*cl*bl-2*y1+2*y3;
   double c = al*al*x1*x1+2*x1*al*cl+cl*cl-y3*y3+y1*y1;
   double k = b*b-4*a*c;
cout<<a<<" "<<b<<" "<<c<<" "<<k<<endl;
   if( k<0 ) return 0;
   double bc = (-b+sqrt(k))/2./a;
   double ac = (-bl*bc-cl)/al;
   cout<<ac<<" "<<bc<<endl;
   bc = (-b-sqrt(k))/2./a;
   ac = (-bl*bc-cl)/al;
   cout<<ac<<" "<<bc<<endl;

   return 1;
}
