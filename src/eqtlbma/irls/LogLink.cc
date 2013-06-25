/** \file LogLink.cc
*
* `IRLS' is a C++ implementation of the IRLS algorithm for GLM
* Copyright (C) 2013 Xioaquan Wen, Timothee Flutre
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

void LogLink::init_mv(gsl_vector * y, gsl_vector * mv)
{
  size_t n = y->size;
  for(size_t i = 0; i < n; ++i){
    if(gsl_vector_get(y, i) == 0)
      gsl_vector_set(mv, i, 0.01);
    else
      gsl_vector_set(mv, i, gsl_vector_get(y, i));
  }
}

void LogLink::compute_z(gsl_vector * y, gsl_vector * mv,
			gsl_vector * offset, gsl_vector * z)
{
  size_t n = y->size;
  double mv_i, y_i, val;
  for(size_t i = 0; i < n; ++i){
    mv_i = gsl_vector_get(mv, i);
    y_i = gsl_vector_get(y, i);
    val = log(mv_i) + (1.0/mv_i)*(y_i-mv_i) - gsl_vector_get(offset, i);
    gsl_vector_set(z, i, val);
  }
}

void LogLink::compute_weights(gsl_vector * mv, gsl_vector * w)
{
  size_t n = mv->size;
  for(size_t i = 0; i < n; ++i)
    gsl_vector_set(w, i, gsl_vector_get(mv, i));
}

void LogLink::compute_mv(gsl_vector * bv, gsl_matrix * Xv,
			 gsl_vector * offset, gsl_vector * mv)
{
  size_t n = Xv->size1, p = Xv->size2;
  
  gsl_matrix * B = gsl_matrix_calloc(p, 1);
  for(size_t i = 0; i < p; ++i)
    gsl_matrix_set(B, i, 0, gsl_vector_get(bv, i));
  
  gsl_matrix * fit = gsl_matrix_calloc(n, 1);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, Xv, B, 0, fit);
  
  for(size_t i = 0; i < n; ++i)
    gsl_vector_set(mv, i, exp(gsl_matrix_get(fit, i, 0) 
			      + gsl_vector_get(offset, i)));
  
  gsl_matrix_free(B);
  gsl_matrix_free(fit);
}

double LogLink::compute_dispersion(gsl_vector * y, gsl_matrix * Xv,
				   gsl_vector * bv, gsl_vector * offset,
				   gsl_vector * mv, double rank,
				   bool quasi_lik)
{
  double psi;
  if(! quasi_lik){
    psi = 1.0;
  }
  else{
    compute_mv(bv, Xv, offset, mv);
    double wtss = 0.0;
    for(size_t i = 0; i < y->size; ++i)
      wtss += pow(gsl_vector_get(y,i) - gsl_vector_get(mv,i), 2) / 
	gsl_vector_get(mv,i);
    psi = (1/(y->size-rank)) * wtss;
  }
  return psi;
}
