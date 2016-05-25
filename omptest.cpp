#include <iostream>
#include <vector>
#include <stdlib.h>
#include <cstring>
#include <cmath>
#include <iomanip>



template<typename type>
class matrix{
public:
  matrix(size_t nrow, size_t ncol){
    data = new type[nrow * ncol];
    _ncol = ncol;
    _nrow = nrow;
  }
  ~matrix(){
    delete [] data;
  }
  size_t _nrow;
  size_t _ncol;
  type *data;
  type * operator[](size_t row)
  { return &data[row*_ncol];}
};

int main(){
  std::cout<<"OMP test\n";
  size_t nrow = 5;
  size_t ncol = nrow * 2;
  size_t dim = nrow * ncol;	  
  matrix<double> mat(nrow, ncol);
  memset(mat.data, 0, sizeof(double) * dim);	
  for(int i = 0; i < nrow; i++){
    mat[i][0] = 1.0/sqrt(nrow);
    for(int j = 1; j < ncol / 2; j++ ){
      double v = cos(i * M_PI * j /nrow);
      mat[i][j] = v;
    }
    mat[i][i+nrow] = 1.0; 
  }
    	
//zero mean and normalize it

  for(int i = 1; i < ncol; i++){
    double sum = 0.0;
    double mean = 0.0;
    double norm = 0.0;
    for(int j = 0; j < nrow; j++){
      sum += mat[j][i];
    }
    mean = sum / nrow;
    std::cout<<mean<<" mean\n";
    for(int j = 0; j < nrow; j++){
      mat[j][i] -= mean;
      norm += mat[j][i] * mat[j][i];
    }
    norm = sqrt(norm);
    for(int j = 0; j < nrow; j++){
      mat[j][i] /= norm;
    }
  }


  for(int i = 0; i < nrow; i++){
    for(int j = 0; j < ncol; j++){
      std::cout.width(8);
      std::cout<<std::right<<std::setw(10)<<mat[i][j]<<" ";
    }
    std::cout<<"\n";
  } 

   



}
