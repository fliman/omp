#include <iostream>
#include <vector>
#include <stdlib.h>
#include <cstring>
#include <cmath>



int main(){
   std::cout<<"OMP test\n";
   size_t nrow = 256;
   size_t ncol = 256 * 2;
   size_t dim = nrow * ncol;	  
   double *D = new double[dim];
   memset(D, 0, sizeof(double) * dim);	
   for(int i = 0; i < nrow; i++){
	D[i*ncol] = 1 / sqrt(nrow);
        for(int j = 1; j < ncol / 2; j++ ){
            double v = cos(i * M_PI * j /nrow);
            D[i*ncol + j] = v;
        }
	D[i*ncol + nrow + i] = 1.0;
   }
    	

   

   delete [] D;

}
