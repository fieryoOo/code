#include <iostream>
#include "List.h"


int main() {
   List<int> A;

   A.Append(-183);
   A.Append('a');
   A.Append('b');
   std::cout<<"1st: "<<A.GetFirst()<<std::endl;
   std::cout<<"last: "<<A.GetLast()<<std::endl;
   return 0;
}
