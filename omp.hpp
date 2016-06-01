#ifndef __OMP_HPP__
#define __OMP_HPP__
#include <vector>
#include "omputil.hpp"
#include "myblas.hpp"
#include <math.h>
#include <string.h>

template<typename type>
void
ompcore(
	type *D,
	type *x,
	type *DtX, 
	type *XtX, 
	type *G, 
	size_t n, 
	size_t m, 
	size_t L,
	size_t maxatoms, 
    type eps,
    size_t T,
    int erroromp,
    type *gamma
){
  	

  size_t i, j, signum, pos, *ind, *gammaIr, *gammaJc, gamma_count;
  size_t allocated_coefs, allocated_cols;
  size_t *selected_atoms;
  type *alpha, *r, *Lchol, *c, *Gsub, *Dsub, sum, *gammaPr, *tempvec1, *tempvec2; 
  type eps2, resnorm, delta, deltaprev, secs_remain;
  int mins_remain, hrs_remain;
  
 
  /*** status flags ***/
  
  bool DtX_specified = (DtX!=0);   /* indicates whether D'*x was provided */
  bool XtX_specified = (XtX!=0);   /* indicates whether sum(x.*x) was provided */
  
  bool standardomp = (G==0);       /* batch-omp or standard omp are selected depending on availability of G */
  bool batchomp = !standardomp;
  
  
  //gamma must allocated outside
  
  
  /*** helper arrays ***/
  
  alpha = new type[m];        /* contains D'*residual */
  ind = new size_t[n];        /* indices of selected atoms */
  selected_atoms = new size_t[m];    /* binary array with 1's for selected atoms */
  c = new type[n];            /* orthogonal projection result */
  
  /* current number of columns in Dsub / Gsub / Lchol */
  allocated_cols = erroromp ? (size_t)(ceil(sqrt((double)n)/2.0) + 1.01) : T;
  
  /* Cholesky decomposition of D_I'*D_I */
  Lchol = new type[n*allocated_cols];
  /* temporary vectors for various computations */
  tempvec1 = new type[m];
  tempvec2 = new type[m];
  
  if (batchomp) {
    /* matrix containing G(:,ind) - the columns of G corresponding to the selected atoms, in order of selection */
    Gsub = new type[m * allocated_cols];
  }else {
    /* matrix containing D(:,ind) - the selected atoms from D, in order of selection */
    Dsub = new type[n * allocated_cols];
    /* stores the residual */
    r = new type[n];        
  }
  
  if (!DtX_specified) {
    /* contains D'*x for the current signal */
    DtX = new type[m];  
  }
  
  
  
  /*** initializations for error omp ***/
  
  if (erroromp) {
    eps2 = eps*eps;        /* compute eps^2 */
    if (T == 0 || T>n) {      /* unspecified max atom num - set max atoms to n */
      T = n;
    }
  }
  
  /**********************   perform omp for each signal   **********************/
  
  
  
  for (signum=0; signum<L; ++signum) {
       
    /* initialize residual norm and deltaprev for error-omp */
    
    if (erroromp) {
      if (XtX_specified) {
        resnorm = XtX[signum];
      }
      else {
        resnorm = dotprod(x+n*signum, x+n*signum, n);
      }
      deltaprev = 0;     /* delta tracks the value of gamma'*G*gamma */
    }
    else {
      /* ignore residual norm stopping criterion */
      eps2 = 0;
      resnorm = 1;
    }
    
    
    if (resnorm>eps2 && T>0) {
      
      /* compute DtX */
      
      if (!DtX_specified) {
        matT_vec((double)1.0, D, x+n*signum, DtX, n, m);
      }
      
      
      /* initialize alpha := DtX */
      
      memcpy(alpha, DtX + m*signum*DtX_specified, m*sizeof(type));
      
      
      /* mark all atoms as unselected */
      
      for (i=0; i<m; ++i) {
        selected_atoms[i] = 0;
      }
      
    }
    

    /* main loop */
    
    i=0;
    while (resnorm>eps2 && i<T) {

      /* index of next atom */
      
      pos = maxabs(alpha, m);
           
      /* stop criterion: selected same atom twice, or inner product too small */
      
      if (selected_atoms[pos] || alpha[pos]*alpha[pos]<1e-14) {
        break;
      }
      
      
      /* mark selected atom */
      
      ind[i] = pos;
      selected_atoms[pos] = 1;
      
      
      /* matrix reallocation */
      
      if (erroromp && i>=allocated_cols) {
        
        allocated_cols = (size_t)(ceil(allocated_cols*MAT_INC_FACTOR) + 1.01);
        
        //need to improve
        delete [] Lchol;
        Lchol = new type[n*allocated_cols];
        
        batchomp ? (Gsub = new type[m*allocated_cols]) :
                   (Dsub = new type[n*allocated_cols]) ;
      }
      
      
      /* append column to Gsub or Dsub */
      
      if (batchomp) {
        memcpy(Gsub+i*m, G+pos*m, m*sizeof(type));
      }
      else {
        memcpy(Dsub+i*n, D+pos*n, n*sizeof(type));
      }
      
      
      /*** Cholesky update ***/
      
      if (i==0) {
        *Lchol = 1;
      }
      else {
        
        /* incremental Cholesky decomposition: compute next row of Lchol */
        
        if (standardomp) {
          matT_vec(1.0, Dsub, D+n*pos, tempvec1, n, i);      /* compute tempvec1 := Dsub'*d where d is new atom */
        }
        else {
          vec_assign(tempvec1, Gsub+i*m, ind, i);          /* extract tempvec1 := Gsub(ind,i) */
        }
        backsubst('L', Lchol, tempvec1, tempvec2, n, i);   /* compute tempvec2 = Lchol \ tempvec1 */
        for (j=0; j<i; ++j) {                              /* write tempvec2 to end of Lchol */
          Lchol[j*n+i] = tempvec2[j];
        }
        
        /* compute Lchol(i,i) */
        sum = 0;
        for (j=0; j<i; ++j) {         /* compute sum of squares of last row without Lchol(i,i) */
          sum += SQR(Lchol[j*n+i]);
        }
        if ( (1-sum) <= 1e-14 ) {     /* Lchol(i,i) is zero => selected atoms are dependent */
          break;
        }
        Lchol[i*n+i] = sqrt(1-sum);
      }
      

      i++;
      
      
      /* perform orthogonal projection and compute sparse coefficients */
      
      vec_assign(tempvec1, DtX + m*signum*DtX_specified, ind, i);   /* extract tempvec1 = DtX(ind) */
      cholsolve('L', Lchol, tempvec1, c, n, i);                     /* solve LL'c = tempvec1 for c */
      

      /* update alpha = D'*residual */
      
      if (standardomp) {
        mat_vec(-1.0, Dsub, c, r, n, i);             /* compute r := -Dsub*c */
        vec_sum(1.0, x+n*signum, r, n);              /* compute r := x+r */
        
        
        //memcpy(r, x+n*signum, n*sizeof(double));   /* assign r := x */
        //mat_vec1(-1, Dsub, c, 1, r, n, i);         /* compute r := r-Dsub*c */
        
        matT_vec(1.0, D, r, alpha, n, m);            /* compute alpha := D'*r */
       
        /* update residual norm */
        if (erroromp) {
          resnorm = dotprod(r, r, n);
          
        }
      }
      else {
        mat_vec(1.0, Gsub, c, tempvec1, m, i);                              /* compute tempvec1 := Gsub*c */
        memcpy(alpha, DtX + m*signum*DtX_specified, m*sizeof(double));    /* set alpha = D'*x */
        vec_sum(-1.0, tempvec1, alpha, m);                                  /* compute alpha := alpha - tempvec1 */
        
        
        /* update residual norm */
        if (erroromp) {
          vec_assign(tempvec2, tempvec1, ind, i);      /* assign tempvec2 := tempvec1(ind) */
          delta = dotprod(c,tempvec2,i);               /* compute c'*tempvec2 */
          resnorm = resnorm - delta + deltaprev;       /* residual norm update */
          deltaprev = delta;
          
        }
      }
    }
    
    
    /*** generate output vector gamma ***/

    //if (gamma_mode == FULL_GAMMA) {    /* write the coefs in c to their correct positions in gamma */
      for (j=0; j<i; ++j) {
	gamma[m*signum + ind[j]] = c[j];
      }
    //}
    // else {
       /* sort the coefs by index before writing them to gamma */
    //   quicksort(ind,c,i);
      
      
       /* gamma is full - reallocate */
    //   if (gamma_count+i >= allocated_coefs) {
        
    //     while(gamma_count+i >= allocated_coefs) {
    //       allocated_coefs = (size_t)(ceil(GAMMA_INC_FACTOR*allocated_coefs) + 1.01);
    //     }
        
    //     mxSetNzmax(Gamma, allocated_coefs);
    //     mxSetPr(Gamma, mxRealloc(gammaPr, allocated_coefs*sizeof(double)));
    //     mxSetIr(Gamma, mxRealloc(gammaIr, allocated_coefs*sizeof(mwIndex)));
        
    //     gammaPr = mxGetPr(Gamma);
    //     gammaIr = mxGetIr(Gamma);
    //   }
      
      /* append coefs to gamma and update the indices */
    //   for (j=0; j<i; ++j) {
    //     gammaPr[gamma_count] = c[j];
    //     gammaIr[gamma_count] = ind[j];
    //     gamma_count++;
    //   }
    //   gammaJc[signum+1] = gammaJc[signum] + i;
    // }
    
    
    
    /*** display status messages ***/
    
    // if (msg_delta>0 && (clock()-lastprint_time)/(double)CLOCKS_PER_SEC >= msg_delta)
    // {
    //   lastprint_time = clock();
      
    //   /* estimated remainig time */
    //   secs2hms( ((L-signum-1)/(double)(signum+1)) * ((lastprint_time-starttime)/(double)CLOCKS_PER_SEC) ,
    //     &hrs_remain, &mins_remain, &secs_remain);
      
    //   mexPrintf("omp: signal %d / %d, estimated remaining time: %02d:%02d:%05.2f\n",        
    //     signum+1, L, hrs_remain, mins_remain, secs_remain);
    //   mexEvalString("drawnow;");
    // }
    
  }
  
  /* end omp */
  
  
  
  /* free memory */
  
  if (!DtX_specified) {
    delete [] DtX;
  }
  if (standardomp) {
    delete [] r;
    delete [] Dsub;
  }
  else {
    delete [] Gsub;
  }  
  delete [] tempvec2;
  delete [] tempvec1;
  delete [] Lchol;
  delete [] c;
  delete [] selected_atoms;
  delete [] ind;
  delete [] alpha;
}

/*
Error-constrained Orthogonal Matching Pursuit
Fast version use more memory
*/
template<typename type>
void
omp_ec(
       //	matrix<type> &DtX,
        type *DtX,  
	std::vector<type> &XtX,
	//	matrix<type> &G,
	type *G,
	type epsilon,
	size_t maxatoms,
	size_t m,
	size_t n,
	size_t L,
	//	matrix<type> &gamma_ou
	type *gamma_ou
){
  /*
	size_t m = DtX._nrow;
	size_t n;
	size_t L = DtX._ncol;
	if (maxatoms>0) {
		n = maxatoms;
	}else {
		n = m;
	}
  */
	ompcore(
		(double*)NULL, 
		(double*)NULL, 
		DtX, 
		(type *)&XtX[0], 
		G, 
		n, 
		m, 
		L, 
		maxatoms, 
		epsilon,
		0, 
		1,
		gamma_ou);


}

#endif
