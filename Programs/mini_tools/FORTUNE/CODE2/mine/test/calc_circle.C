#include <math.h>
#include <iostream.h>

int main(int argc, char ** argv)
{
   if( argc != 6 ) {
      cout<<"Usage: "<<argv[0]<<" [x1] [y1] [x2] [y2] [y3]"<<endl;
      exit(-1);
   }

   double x1 = atof(argv[1]), y1 = atof(argv[2]), x2 = atof(argv[3]), y2 = atof(argv[4]), y3 = atof(argv[5]);
   double a = y2-y1;
   double b = 2*(x1*(y3-y2)+x2*(y1-y3));
   double c = -(y3-y1)*(y3-y2)*(y2-y1)+x1*x1*(y2-y3)-x2*x2*(y1-y3);
   double k = b*b-4*a*c;
cout<<a<<" "<<b<<" "<<c<<" "<<k<<endl;
   if( k<0 ) return 0;
   double ac = (-b+sqrt(k))/2./a;
   cout<<ac<<endl;
   ac = (-b-sqrt(k))/2./a;
   cout<<ac<<endl;

   return 1;
}
