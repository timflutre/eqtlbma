/** \file IRLS.cc
*
* `IRLS' is a C++ implementation of the IRLS algorithm for GLM
* Copyright (C) 2013 Xioaquan Wen
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <cstring>

#include <gsl/gsl_multifit.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include "IRLS.h"
#include "LogLink.h"

using namespace std;


IRLS::IRLS(const char * link_type){
  if(strcmp(link_type,"log-link")==0){
    link = new LogLink();
  }
  link->quasi = false;
  fit_coef = 0;
  VB = 0;
}

// Xv should not contain the intercept
void IRLS::load_data(vector<double> & yv, vector<vector<double> > &Xv){
  
  free_data = true;
  
  n = yv.size();
  p = Xv.size()+1;
  
  Y = gsl_vector_calloc(n);
  X = gsl_matrix_calloc(n,p);
  
  for(int i=0;i<n;i++){
    gsl_vector_set(Y,i,yv[i]);
    for(int j=0;j<p;j++){
      double val = 1;
      if(j>=1)
	val = Xv[j-1][i];
      
      gsl_matrix_set(X,i,j,val);
    }
  }
  
  return;
}

// Xv should contain the intercept
void IRLS::set_data(gsl_vector * yv, gsl_matrix * Xv){

  free_data = false;
  
  n = yv->size;
  p = Xv->size2;
  
  Y = yv;
  X = Xv;
}
 

void IRLS::fit_model(){
  
  gsl_vector *mv = link->init_mv(Y);
  
  double old_chisq = -1;
  while(1){
    
    gsl_vector *z = link->compute_Z(Y,mv);
    gsl_vector *w = link->compute_weights(mv);
  
    // weighted least square fitting
    gsl_vector *bv = gsl_vector_alloc(p);
    gsl_matrix *cov = gsl_matrix_alloc(p, p);
    gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc (n, p);
    double chisq;
    gsl_multifit_wlinear_svd (X, w, z, GSL_DBL_EPSILON, &rank, bv, cov, &chisq, work);
    
    gsl_vector_free(z);
    
    gsl_matrix_free (cov);
    gsl_multifit_linear_free (work);
    gsl_vector_free(mv);


    //convergence diagnosis here
    if(fabs(chisq - old_chisq)<1e-6){
 
      if(fit_coef !=0){
	gsl_vector_free(fit_coef);
      }
      
      psi = link->compute_dispersion(Y, X, bv, rank, link->quasi);
      compute_variance(w);
      fit_coef = bv;
      gsl_vector_free (w);
      break;
    }
    
    old_chisq = chisq;
    mv = link->compute_mv(bv,X);
    gsl_vector_free(bv);
    gsl_vector_free (w);
  }
  

}


void IRLS::compute_variance(gsl_vector *w){
  
  if(VB!=0)
    gsl_matrix_free(VB);
  
  VB = gsl_matrix_calloc(p,p);
  gsl_matrix *W = gsl_matrix_calloc(n,n);
  for(int i=0;i<n;i++){
    gsl_matrix_set(W,i,i,gsl_vector_get(w,i));
  }

  gsl_matrix *t1 = gsl_matrix_calloc(p,n);
  gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,X,W,0,t1);
  gsl_matrix *t2 = gsl_matrix_calloc(p,p);
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,t1,X,0,t2);
  
  // invert t2
  int ss;
  gsl_permutation * pp = gsl_permutation_alloc(p);
  gsl_linalg_LU_decomp (t2, pp, &ss);
  gsl_linalg_LU_invert (t2, pp, VB);
 
  gsl_permutation_free(pp);
 
  gsl_matrix_free(W);
  gsl_matrix_free(t1);
  gsl_matrix_free(t2);
 
}



vector<double> IRLS::get_fit_coef(){
  
  vector<double> coev;
  for(int i=0;i<p;i++){
    coev.push_back(gsl_vector_get(fit_coef,i));
  }
  return coev;
}


vector<double> IRLS::get_stderr(){
  
  vector<double> sev;
  for(int i=0;i<p;i++){
    sev.push_back(sqrt(psi * gsl_matrix_get(VB,i,i)));
  }
  return sev;
}
  


IRLS::~IRLS(){
  
  delete link;
  
  if(free_data){
    gsl_vector_free(Y);
    gsl_matrix_free(X);
  }
  
  if(fit_coef !=0){
    gsl_vector_free(fit_coef);
  }
  
  if(VB != 0)
    gsl_matrix_free(VB);
}
