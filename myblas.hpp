/**************************************************************************
 *
 * File name: myblas.c
 *
 * Ron Rubinstein
 * Computer Science Department
 * Technion, Haifa 32000 Israel
 * ronrubin@cs
 *
 * Version: 1.1
 * Last updated: 13.8.2009
 *
 *************************************************************************/


/* find maximum of absolute values */
template<typename type>
size_t maxabs(type c[], size_t m)
{
  size_t maxid=0, k;
  type absval, maxval = SQR(*c);   /* use square which is quicker than absolute value */

  for (k=1; k<m; ++k) {
    absval = SQR(c[k]);
    if (absval > maxval) {
      maxval = absval;
      maxid = k;
    }
  }
  return maxid;
}


/* compute y := alpha*x + y */
template<typename type>
void vec_sum(type alpha, type x[], type y[], size_t n)
{
  size_t i;

  for (i=0; i<n; ++i) {
    y[i] += alpha*x[i];
  }
}


/* compute y := alpha*A*x */
template<typename type>
void mat_vec(type alpha, type A[], type x[], type y[], size_t n, size_t m)
{
  size_t i, j, i_n;
  type *Ax;

  Ax = new type[n];

  for (i=0; i<m; ++i) {
    i_n = i*n;
    for (j=0; j<n; ++j) {
      Ax[j] += A[i_n+j] * x[i];
    }
  }

  for (j=0; j<n; ++j) {
    y[j] = alpha*Ax[j];
  }

  delete [] Ax;
}


/* compute y := alpha*A'*x */
template<typename type>
void matT_vec(type alpha, type A[], type x[], type y[], size_t n, size_t m)
{
  size_t i, j, n_i;
  type sum0, sum1, sum2, sum3;

  for (j=0; j<m; ++j) {
    y[j] = 0;
  }

  /* use loop unrolling to accelerate computation */

  for (i=0; i<m; ++i) {
    n_i = n*i;
    sum0 = sum1 = sum2 = sum3 = 0;
    for (j=0; j+4<n; j+=4) {
      sum0 += A[n_i+j]*x[j];
      sum1 += A[n_i+j+1]*x[j+1];
      sum2 += A[n_i+j+2]*x[j+2];
      sum3 += A[n_i+j+3]*x[j+3];
    }
    y[i] += alpha * ((sum0 + sum1) + (sum2 + sum3));
    while (j<n) {
      y[i] += alpha*A[n_i+j]*x[j];
      j++;
    }
  }
}


/* matrix-matrix multiplication */
template<typename type>
void mat_mat(type alpha, type A[], type B[], type X[], size_t n, size_t m, size_t k)
{
  size_t i1, i2, i3, iX, iA, i2_n;
  type b;
  
  for (i1=0; i1<n*k; i1++) {
    X[i1] = 0;
  }

  for (i2=0; i2<m; ++i2) {
    i2_n = i2*n;
    iX = 0;
    for (i3=0; i3<k; ++i3) {
      iA = i2_n;
      b = B[i2+i3*m];
      for (i1=0; i1<n; ++i1) {
        X[iX++] += A[iA++]*b;
      }
    }
  }
  
  for (i1=0; i1<n*k; i1++) {
    X[i1] *= alpha;
  }
}


/* matrix-transpose-matrix multiplication */
template<typename type>
void matT_mat(type alpha, type A[], type B[], type X[], size_t n, size_t m, size_t k)
{
  size_t i1, i2, i3, iX, iA, i2_n;
  type *x, sum0, sum1, sum2, sum3;

  for (i2=0; i2<m; ++i2) {
    for (i3=0; i3<k; ++i3) {
      sum0 = sum1 = sum2 = sum3 = 0;
      for (i1=0; i1+4<n; i1+=4) {
        sum0 += A[i1+0+i2*n]*B[i1+0+i3*n];
        sum1 += A[i1+1+i2*n]*B[i1+1+i3*n];
        sum2 += A[i1+2+i2*n]*B[i1+2+i3*n];
        sum3 += A[i1+3+i2*n]*B[i1+3+i3*n];
      }
      X[i2+i3*m] = (sum0+sum1) + (sum2+sum3);
      while(i1<n) {
        X[i2+i3*m] += A[i1+i2*n]*B[i1+i3*n];
        i1++;
      }
    }
  }
  
  for (i1=0; i1<m*k; i1++) {
    X[i1] *= alpha;
  }
}


/* tensor-matrix product */
template<typename type>
void tens_mat(type alpha, type A[], type B[], type X[], size_t n, size_t m, size_t k, size_t l)
{
  size_t i1, i2, i3, i4, i2_n, nml;
  type b;
  
  nml = n*m*l;
  for (i1=0; i1<nml; ++i1) {
    X[i1] = 0;
  }

  for (i2=0; i2<m; ++i2) {
    i2_n = i2*n;
    for (i3=0; i3<k; ++i3) {
      for (i4=0; i4<l; ++i4) {
        b = B[i4+i3*l];
        for (i1=0; i1<n; ++i1) {
          X[i1 + i2_n + i4*n*m] += A[i1 + i2_n + i3*n*m] * b;
        }
      }
    }
  }
  
  for (i1=0; i1<nml; ++i1) {
    X[i1] *= alpha;
  }
}


/* tensor-matrix-transpose product */
template<typename type>
void tens_matT(type alpha, type A[], type B[], type X[], size_t n, size_t m, size_t k, size_t l)
{
  size_t i1, i2, i3, i4, i2_n, nml;
  type b;
  
  nml = n*m*l;
  for (i1=0; i1<nml; ++i1) {
    X[i1] = 0;
  }

  for (i2=0; i2<m; ++i2) {
    i2_n = i2*n;
    for (i4=0; i4<l; ++i4) {
      for (i3=0; i3<k; ++i3) {
        b = B[i3+i4*k];
        for (i1=0; i1<n; ++i1) {
          X[i1 + i2_n + i4*n*m] += A[i1 + i2_n + i3*n*m] * b;
        }
      }
    }
  }
  
  for (i1=0; i1<nml; ++i1) {
    X[i1] *= alpha;
  }
}


/* dot product */
template<typename type>
type dotprod(type a[], type b[], size_t n)
{
  type sum = 0;
  size_t i;
  for (i=0; i<n; ++i)
    sum += a[i]*b[i];
  return sum;
}


/* find maximum of vector */
template<typename type>
size_t maxpos(type c[], size_t m)
{
  size_t maxid=0, k;
  type val, maxval = *c;

  for (k=1; k<m; ++k) {
    val = c[k];
    if (val > maxval) {
      maxval = val;
      maxid = k;
    }
  }
  return maxid;
}


/* solve L*x = b */
template<typename type>
void backsubst_L(type L[], type b[], type x[], size_t n, size_t k)
{
  size_t i, j;
  type rhs;

  for (i=0; i<k; ++i) {
    rhs = b[i];
    for (j=0; j<i; ++j) {
      rhs -= L[j*n+i]*x[j];
    }
    x[i] = rhs/L[i*n+i];
  }
}


/* solve L'*x = b */
template<typename type>
void backsubst_Lt(type L[], type b[], type x[], size_t n, size_t k)
{
  size_t i, j;
  type rhs;

  for (i=k; i>=1; --i) {
    rhs = b[i-1];
    for (j=i; j<k; ++j) {
      rhs -= L[(i-1)*n+j]*x[j];
    }
    x[i-1] = rhs/L[(i-1)*n+i-1];
  }
}


/* solve U*x = b */
template<typename type>
void backsubst_U(type U[], type b[], type x[], size_t n, size_t k)
{
  size_t i, j;
  type rhs;

  for (i=k; i>=1; --i) {
    rhs = b[i-1];
    for (j=i; j<k; ++j) {
      rhs -= U[j*n+i-1]*x[j];
    }
    x[i-1] = rhs/U[(i-1)*n+i-1];
  }
}


/* solve U'*x = b */
template<typename type>
void backsubst_Ut(type U[], type b[], type x[], size_t n, size_t k)
{
  size_t i, j;
  type rhs;

  for (i=0; i<k; ++i) {
    rhs = b[i];
    for (j=0; j<i; ++j) {
      rhs -= U[i*n+j]*x[j];
    }
    x[i] = rhs/U[i*n+i];
  }
}


/* back substitution solver */
template<typename type>
void backsubst(char ul, type A[], type b[], type x[], size_t n, size_t k)
{
  if (tolower(ul) == 'u') {
    backsubst_U(A, b, x, n, k);
  }
  else if (tolower(ul) == 'l') {
    backsubst_L(A, b, x, n, k);
  }
  else {
    std::cout<<"Invalid triangular matrix type: must be ''U'' or ''L''\n";
  }
}


/* solve equation set using cholesky decomposition */
template<typename type>
void cholsolve(char ul, type A[], type b[], type x[], size_t n, size_t k)
{
  type *tmp;

  tmp = new type[k];
  if (tolower(ul) == 'l') {
    backsubst_L(A, b, tmp, n, k);
    backsubst_Lt(A, tmp, x, n, k);
  }
  else if (tolower(ul) == 'u') {
    backsubst_Ut(A, b, tmp, n, k);
    backsubst_U(A, tmp, x, n, k);
  }
  else {
    std::cout<<"Invalid triangular matrix type: must be either ''U'' or ''L''\n";
  }

  delete [] tmp;
}


/* perform a permutation assignment y := x(ind(1:k)) */
template<typename type>
void vec_assign(type y[], type x[], size_t ind[], size_t k)
{
  size_t i;

  for (i=0; i<k; ++i)
    y[i] = x[ind[i]];
}


/* matrix transpose */
template<typename type>
void transpose(type X[], type Y[], size_t n, size_t m)
{
  size_t i, j, i_m, j_n;
  
  if (n<m) {
    for (j=0; j<m; ++j) {
      j_n = j*n;
      for (i=0; i<n; ++i) {
        Y[j+i*m] = X[i+j_n];
      }
    }
  }
  else {
    for (i=0; i<n; ++i) {
      i_m = i*m;
      for (j=0; j<m; ++j) {
        Y[j+i_m] = X[i+j*n];
      }
    }
  }
}






/* matrix multiplication using Winograd's algorithm */

/*
void mat_mat2(double alpha, double A[], double B[], double X[], mwSize n, mwSize m, mwSize k)
{
  
  size_t i1, i2, i3, iX, iA, i2_n;
  double b, *AA, *BB;
  
  AA = mxCalloc(n,sizeof(double));
  BB = mxCalloc(k,sizeof(double));
  
  for (i1=0; i1<n*k; i1++) {
    X[i1] = 0;
  }
  
  for (i1=0; i1<n; ++i1) {
    for (i2=0; i2<m/2; ++i2) {
      AA[i1] += A[i1+2*i2*n]*A[i1+(2*i2+1)*n];
    }
  }

  for (i2=0; i2<k; ++i2) {
    for (i1=0; i1<m/2; ++i1) {
      BB[i2] += B[2*i1+i2*m]*B[2*i1+1+i2*m];
    }
  }

  for (i2=0; i2<k; ++i2) {
    for (i3=0; i3<m/2; ++i3) {
      for (i1=0; i1<n; ++i1) {
        X[i1+i2*n] += (A[i1+(2*i3)*n]+B[2*i3+1+i2*m])*(A[i1+(2*i3+1)*n]+B[2*i3+i2*m]);
      }
    }
  }
  
  if (m%2) {
    for (i2=0; i2<k; ++i2) {
      for (i1=0; i1<n; ++i1) {
        X[i1+i2*n] += A[i1+(m-1)*n]*B[m-1+i2*m];
      }
    }
  }
  
  for (i2=0; i2<k; ++i2) {
    for (i1=0; i1<n; ++i1) {
      X[i1+i2*n] -= (AA[i1] + BB[i2]);
      X[i1+i2*n] *= alpha;
    }
  }
  
  mxFree(AA);
  mxFree(BB);
}
*/




