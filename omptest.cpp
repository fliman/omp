#include <iostream>
#include <vector>
#include <stdlib.h>
#include <cstring>
#include <cmath>
#include <iomanip>
#include <vector>
#include "omp.hpp"
#include <vector>



int main(){
  std::cout<<"OMP test\n";
  size_t nrow = 256;
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

  std::vector<double> g(2*nrow, 0.0);  

  int lowfreq = 20;

  int p[] = {13, 17,  1, 18,  7,  4, 10,  8,  2,  9,  6,  3,  0, 19, 11, 14,  5,
       16, 12, 15}; 

  std::vector<double> v1(2,0.0);
  v1[0] = (0.0606794*2  -0.5) *  -1.0;
  v1[1] = (0.41122722*2 -0.5 ) * 1.0;     

  g[p[0]] = v1[0];
  g[p[1]] = v1[1];

  int p2[] = {230,  92, 152,  19, 144,  49,  55, 156,  29,  86,  27,  98, 110,
       136, 174,  89, 234, 120, 210, 178, 181,  80, 205, 154,  63,  25,
        31, 111,   3, 219, 141, 106, 189,  96, 149, 117, 177,  59,  82,
       172,  45,  65, 183, 218, 203, 184, 103,  51, 199,  85, 142, 112,
        14,  12, 130, 128, 215, 134, 186,  60, 180,  39, 150, 201,  15,
         4, 224, 209,  36, 191, 135,  79, 168,  73, 188,  93,  47,  32,
         0, 200, 148,   2, 226, 104,  78, 202,  83,  43, 225, 222, 131,
       206, 157, 187, 166, 122, 207,  97,  91, 221,  68, 212,  67, 153,
        62,  81, 123, 151, 109,  94, 232,   6, 208, 124, 155,   9,  77,
        16, 115,  76, 227, 220, 163, 108, 121, 216, 231, 185,  72, 102,
       132, 162,   7,  57,   8, 167, 126,  99, 164, 195,  42,  34, 179,
        28, 133,  52,  90,  66,  54, 125, 145,  53,  61, 143,  69, 194,
       137, 192, 161,  37, 100, 233, 129,  64,  10, 171, 158, 165,  46,
       169,  41, 105, 170, 160, 217,  22,  48,  74,  26,  13,  40, 147,
       229, 193, 140, 204,  24, 107, 138,  84,  56, 139,  38,  75, 197,
       190,  21,  20, 196,  23,  33, 119, 214,  50, 175, 211, 114,  17,
       113, 235,  58, 116, 213,  44, 146,  11, 182,  18,  35, 159,  30,
         5, 176, 223,  71, 228,  95,   1, 127, 118, 198,  87,  70,  88,
       173, 101};

  std::cout<<"sizeof p2 "<<sizeof(p2)/sizeof(float)<<"\n";
  
  g[20+p2[0]] = (0.1232521 * 2 - 0.5) * 1.0;
  g[20+p2[1]] = (0.9081769 * 2 - 0.5) * 1.0;    

  int p3[] = {190,  51, 127, 187, 167, 209, 251, 124, 198,  76,  87,  84,  49,
       242, 121,  62, 118,  17, 100,  36,  83, 168,  15,  96, 141, 102,
        81,  24, 230,  27,  41, 222,   9,  45, 145,  72,  20, 241, 238,
        19,  80, 115,  21,  31, 192,  92, 218,  46, 153,  35, 219,  54,
       207,  78, 169, 159,  75,  65, 214,  70, 158,  61, 227, 245, 211,
        22, 221,  53,  38, 188, 137, 212,  42, 132,  64, 183, 226, 133,
        50,  79, 107, 254,  95, 163,  55, 213,  56, 250, 113,  16,  33,
        26, 164, 130, 205,  58, 112,   0, 233, 123,  98,  18, 180,  52,
       152, 154, 128,  60, 147,   2, 181,  66,  99, 249, 191,  11, 165,
       216,  25,   3, 224,  73, 196, 143,  34, 156, 110,  13,  93, 220,
       252,   4,  71, 232, 166,  39, 199,  67, 125,  97,  59, 122, 126,
       235,  28, 236, 114, 197,  48, 151, 134,  63, 103, 109, 200, 225,
        47,  23, 240, 176, 223, 120, 210,  74, 111,  85,  43, 247, 144,
        69,   8, 194,   5,  44, 129,  68,  40,  29, 117, 136, 228,  12,
        57, 239, 208,  90,  89,  82, 150, 138, 248,  77, 162, 255, 160,
       139, 177,  32, 106, 173, 171, 244, 243, 195, 172, 140, 148, 105,
       108, 146, 179, 237,   6, 246, 215, 170,  88, 155, 135, 182, 231,
       119, 186, 185,   1, 175, 189, 253,  14,  37, 178, 201, 202, 229,
         7, 204, 157, 149,  30, 104, 161,  91, 174,  86, 206, 217, 101,
       203, 234,  10, 116, 142, 193,  94, 184, 131};


  g[nrow+p3[0]] = (0.46045077  - 0.25) * 1.0;
  g[nrow+p3[1]] = (0.10070385  - 0.25) * -1.0; 
  g[nrow+p3[2]] = (0.93948024  - 0.25) * 1.0;

  std::vector<double> x = mat_vec<double>(mat, g);   

  double rand[] = {2.70337337e-01,  -3.12445315e-01,   6.70319619e-01,
         -1.54748962e+00,  -2.18957346e-01,   3.25104887e-01,
         -3.15215569e-01,  -7.23852991e-01,  -9.04395656e-01,
          1.44351126e-01,  -4.79160470e-01,  -7.91927387e-01,
          3.88352866e-01,   1.60149289e-01,   2.31674459e-01,
         -3.11753947e-01,  -1.52018347e+00,   1.17767768e+00,
          2.09952271e-01,  -7.54079675e-01,   2.22036538e-01,
          3.30139686e-01,  -7.05126185e-01,  -4.61087809e-02,
         -3.92389677e-01,  -4.00938718e-02,   1.01833944e-01,
          9.91582613e-02,  -6.82665209e-01,  -3.24065747e+00,
          9.79539165e-01,   6.03392331e-01,   1.73925721e+00,
         -9.62657915e-01,   1.12128512e+00,  -3.19556504e-01,
          7.18890924e-01,   1.91340975e+00,  -5.15273432e-01,
         -1.34125538e-02,  -1.09569409e+00,  -2.24790338e-01,
          1.95899014e-02,  -3.41419819e-01,   6.20796189e-01,
          5.09284838e-01,   7.79287310e-01,   1.79474671e+00,
          6.02970933e-01,   1.27364791e+00,   8.25984080e-01,
         -4.44075483e-01,   6.91537829e-01,   1.12963541e+00,
         -1.29173219e+00,   3.45500165e-01,   2.60733758e-01,
          2.44268430e-01,   1.35723192e-01,  -7.88049311e-01,
         -3.35490460e-01,  -3.49883492e-01,   6.32723934e-01,
          7.27235806e-01,   5.16886250e-01,  -2.11031904e-01,
          7.34324465e-01,  -8.84530649e-01,  -2.01087224e-01,
          1.61429173e-01,   3.81457515e-01,  -1.44314641e+00,
          1.18746409e+00,   7.94929608e-01,   2.81894978e-02,
         -1.37840235e+00,  -8.17886931e-01,   2.36331822e-01,
         -8.41297550e-01,  -5.01367688e-02,   2.30530921e+00,
         -9.75636475e-01,  -8.88099870e-01,   7.25399316e-01,
         -6.45621192e-01,   8.06054901e-01,   1.03797618e+00,
          1.48439522e-01,  -1.15997847e+00,  -8.40280194e-01,
         -3.27890593e-01,   1.81727950e+00,  -1.35792573e+00,
          1.61825811e+00,  -2.89402915e-01,  -5.44882473e-01,
          1.19756544e+00,  -1.24223211e+00,   4.19321351e-01,
          1.24477041e+00,  -1.08117893e+00,   5.04169315e-01,
          7.60830419e-01,  -3.60323846e-01,  -3.85075136e-01,
         -1.76219393e+00,   5.01216665e-01,   1.57753581e+00,
          1.66605635e-02,   9.52737627e-01,  -1.49049183e-01,
         -1.12974108e+00,  -5.75269712e-01,  -2.24575819e+00,
          1.79586668e-01,  -1.03244273e+00,   2.48610325e-01,
          2.37788910e+00,   3.14954996e-03,  -7.36481429e-01,
          1.16222379e+00,  -2.39384251e-01,   4.99021744e-01,
         -8.38539225e-01,   4.20740980e-01,   5.06603061e-01,
         -8.04034004e-01,  -1.47417665e+00,  -1.59799653e-01,
          2.13576671e-01,  -3.60021653e-01,  -1.41390214e+00,
         -9.39527572e-02,   3.37387479e-01,   3.14652209e-01,
          1.45830305e+00,  -1.88445104e+00,   1.75268349e+00,
          2.21967574e-01,  -4.75181505e-01,   2.30356562e-01,
         -7.27660727e-01,   8.86723957e-01,   4.79378736e-01,
          7.24450599e-02,  -2.49199804e-01,  -3.75728531e-01,
          1.47544181e+00,   1.20724318e+00,   1.20480652e-01,
         -4.00498361e-01,   6.13863922e-01,  -3.83658606e-01,
          1.59618635e+00,  -7.73497249e-01,   1.51800161e+00,
         -6.53377043e-01,  -8.55723059e-01,  -2.55134059e-01,
          1.92049603e+00,  -3.59715410e-01,  -1.45550811e+00,
          8.00544067e-02,   7.66414528e-01,   1.95249589e-01,
          3.31176450e-01,   5.36982117e-01,   4.61460451e-01,
         -7.19892665e-01,  -6.64100712e-01,   1.84087833e+00,
          9.05631741e-01,   7.29791744e-01,   3.54122008e-02,
          1.56684603e+00,   3.38447827e-01,  -7.95108919e-02,
         -7.62517526e-01,   9.55105723e-01,   1.31791850e-01,
         -1.43753678e-01,  -7.82486650e-01,   2.59323943e-01,
          8.44864811e-01,  -4.50740999e-01,  -5.56659068e-01,
          1.15571156e+00,   2.05859934e+00,   3.79555087e-01,
          3.18917493e-01,  -1.99134589e+00,   4.63197206e-01,
         -1.69953276e+00,   8.48927319e-01,  -1.93617735e+00,
          1.97758414e+00,  -3.51367614e-01,   4.48509606e-02,
          1.26150320e+00,   3.63513328e-01,   1.55683021e+00,
         -6.47013998e-01,   9.84042748e-01,   9.67572289e-01,
          2.31221068e-01,  -1.25433178e-01,   1.26029762e-01,
         -1.70487354e-02,  -9.89908330e-01,  -1.57800842e-01,
          2.88200518e-01,   5.86333603e-01,  -1.21109163e+00,
         -4.91727028e-01,  -4.50114462e-01,   4.01724130e-01,
         -5.68167103e-01,   2.23077149e+00,  -1.70417322e-01,
         -1.96893511e-01,   1.88369155e+00,   9.26793565e-01,
          8.15471641e-01,   1.73730852e-01,   1.11797986e+00,
          1.64696526e+00,   6.88746400e-01,  -1.57388177e+00,
         -1.88131746e+00,   7.71730264e-02,  -2.32997071e-01,
         -4.67220202e-01,  -1.24402816e+00,   1.62876321e+00,
         -8.43903188e-02,  -7.02804912e-01,   1.13514601e+00,
          1.12871887e+00,  -2.56352897e-01,  -3.16443834e-01,
          1.89102158e-01,   7.37319201e-01,  -1.69482370e+00,
          2.24235737e+00,   1.38437018e-01,  -1.44400447e-01,
          2.28158242e-01,  -7.63898285e-01,  -2.83090755e-01,
          3.47560525e-01,   7.97712320e-01,  -1.73680781e+00,
         -1.49569098e-01,  -1.20322410e+00,  -5.74700011e-01,
          1.78813995e-01};

    double scl = l2norm(&x[0], nrow) / 4.0 / l2norm(&rand[0], nrow);
    
    scale(rand, scl, nrow);      

    std::vector<double> vrand (rand, rand + sizeof(rand) / sizeof(double));
    std::vector<double> vy = vec_plus(x, vrand);
    matrix<double> DtD(ncol, ncol);
    std::vector<double> XtX(1, 0.0);

    XtX[0] = l2norm(&vy[0], vy.size());
    XtX[0] *= XtX[0];

    matrix<double> gamma(nrow, 1);
    matrix<double> DtX(ncol, 1);
    std::vector dtx = mat_tran_vec(mat, vy);
    memcpy(DtX.data, &dtx[0], dtx.size()*sizeof(double));
    mat_tran_mat(mat, mat, DtD);
    omp_ec(DtX, XtX, DtD, 0.01, 256, gamma);


}
