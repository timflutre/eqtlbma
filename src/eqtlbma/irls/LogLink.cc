/** \file LogLink.cc
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

#include <cmath>

#include <gsl/gsl_blas.h>

#include "LogLink.h"

using namespace std;

gsl_vector * LogLink::init_mv(gsl_vector *Y){
  size_t n = Y->size;
  
  gsl_vector *mv = gsl_vector_calloc(n);
  
  for(int i=0;i<n;i++){
    double val = gsl_vector_get(Y,i);
    double mval = val;
    if(val == 0){
      mval = 0.01;
    }
    
    gsl_vector_set(mv,i,mval);
  }
  
  return mv;
}



gsl_vector *LogLink::compute_Z(gsl_vector *Y, gsl_vector *mv){
  
  size_t n = Y->size;
  gsl_vector *Z = gsl_vector_calloc(n);
  for(int i=0;i<n;i++){
    double mu = gsl_vector_get(mv,i);
    double y = gsl_vector_get(Y,i);
    double val = log(mu)+(1.0/mu)*(y-mu);
    gsl_vector_set(Z,i,val);
  }

  return Z;
}

gsl_vector *LogLink::compute_weights(gsl_vector *mv){
  
  size_t n = mv->size;
  gsl_vector *w = gsl_vector_calloc(n);
  for(int i=0;i<n;i++){
    double val = gsl_vector_get(mv,i);
    gsl_vector_set(w,i,val);
  }
  return w;

}



gsl_vector *LogLink::compute_mv(gsl_vector *bv, gsl_matrix *Xv){
  size_t n = Xv->size1, p = Xv->size2;
  
  gsl_matrix *beta = gsl_matrix_calloc(p,1);
  
  for(int i=0;i<p;i++){
    gsl_matrix_set(beta,i,0,gsl_vector_get(bv,i));
  }

  gsl_matrix *fit = gsl_matrix_calloc(n,1);
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,Xv,beta,0,fit);
  
  gsl_vector *mv = gsl_vector_calloc(n);

  for(int i=0;i<n;i++){
    double val = gsl_matrix_get(fit,i,0);
    gsl_vector_set(mv,i,exp(val));
  }
  
  gsl_matrix_free(fit);
  gsl_matrix_free(beta);

  return mv;
}




double LogLink::compute_dispersion(gsl_vector *Y, gsl_matrix *Xv, gsl_vector *bv, double rank, bool quasi_lik)
{
  double psi;
  if(! quasi_lik){
    psi = 1.0;
  }
  else{
    gsl_vector *mv = compute_mv(bv, Xv);
    double wtss = 0.0;
    for(size_t i = 0; i < Y->size; ++i)
      wtss += pow(gsl_vector_get(Y,i) - gsl_vector_get(mv,i), 2) / 
	gsl_vector_get(mv,i);
    psi = (1/(Y->size-rank)) * wtss;
    gsl_vector_free(mv);
  }
  return psi;
}
