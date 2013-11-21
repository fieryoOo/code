#include "2DArray.h"
#include <iostream>
#include <cstdlib>

int main() {

   /* initialize */
   Array2D<float> A(5,23,-1.);
   std::cerr<<A.NumRows()<<" "<<A.NumCols()<<" "<<A.Size()<<std::endl;

   /* access elements */
   A(1, 3) = 15.8;
   A(2, 3) = -1.3;
   Array2D<float> B(A);
   Array2D<float> C(A); C.CopySub( Array2D<float>(2,2,-3.), 2, 2 );
   std::cerr<<"A: "<<A(0,0)<<" "<<A(1,3)<<" "<<A(2,3)<<std::endl;
   std::cerr<<"B: "<<B(0,0)<<" "<<B(1,3)<<" "<<B(2,3)<<std::endl;
   std::cerr<<"contents of C:"<<std::endl;
   for(int irow=0; irow<C.NumRows(); irow++) {
      for(int icol=0; icol<C.NumCols(); icol++) std::cerr<<C(irow, icol)<<" ";
      std::cerr<<std::endl;
   }

   /* get a row or a column */
   std::vector<float> rowcur = A.GetRow(2);
   std::vector<float> colcur = A.GetCol(3, 1, 4);
   for(int i=0; i<rowcur.size(); i++) { std::cerr<<rowcur.at(i)<<", "; } std::cerr<<std::endl;
   for(int i=0; i<colcur.size(); i++) std::cerr<<colcur.at(i)<<std::endl;
}
