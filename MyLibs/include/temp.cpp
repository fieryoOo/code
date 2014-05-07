#include "DisAzi.h"

int main() {
   std::cerr<<std::setprecision(15)<<Path<float>(0.,1.,9.,10.).Dist()<<std::endl;
   std::cerr<<std::setprecision(15)<<Path<double>(0.,1.,9.,10.).Dist()<<std::endl;
   std::cerr<<std::setprecision(15)<<Path<double>(Point<double>(0.,1.),Point<double>(9.,10.)).Dist()<<std::endl;
   std::cerr<<std::setprecision(15)<<Path<double>(0.,1.,9.,10.).Azi1()<<std::endl;
   std::cerr<<std::setprecision(15)<<Path<double>(0.,1.,9.,10.).Azi2()<<std::endl;

   return 0;
}
