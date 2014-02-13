#include <cstdio>
#include <iostream>
#include <memory>

class Foo {
   float *array;

public:
   Foo() {
      std::cout<<" Allocating array... "<<std::endl;
      array = new float[100];
   }
   ~Foo() {
      std::cout<<" Deallocating array... "<<std::endl;
      delete [] array;
   }
};

int main() {
   std::cout<<" define p1 pointing to new Foo[1]: "<<std::endl;
   std::unique_ptr<Foo[]> p1(new Foo[1]);
   std::cout<<" define newfoo pointing to new Foo[5]: "<<std::endl;
   Foo *newfoo = new Foo[5];
   std::cout<<" let p2 own newfoo: "<<std::endl;
   std::unique_ptr<Foo[]> p2(newfoo);
   std::cout<<" move p1 to p2: "<<std::endl;
   p2 = std::move(p1);
   std::cout<<" end of main... "<<std::endl;
}
