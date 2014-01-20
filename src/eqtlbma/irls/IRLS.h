/** \file IRLS.h
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

#ifndef _IRLS_H_
#define _IRLS_H_

#include <vector>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include "LinkFunc.h"

class IRLS
{
private:
  gsl_vector * y;      // response vector
  gsl_matrix * X;      // design matrix
  gsl_vector * offset; // offset vector
  bool free_data;      // depends on load_data() or set_data()
  
  size_t n;            // sample size;
  size_t p;            // number of parameters (including intercept)
  size_t rank;         // of X (useful for p-values)
  
  gsl_vector * bv;     // vector of estimated effect sizes
  gsl_matrix * VB;     // covariance matrix of estimated effect sizes
  double psi;          // dispersion
  
public:
  
  LinkFunc * link;
  
  IRLS(const char * link_type);
  ~IRLS();
  void load_data(const std::vector<double> & yv,
		 const std::vector<std::vector<double> > &Xv,
		 const std::vector<double> & offv);
  void set_data(gsl_vector * yv,
		gsl_matrix * Xv,
		gsl_vector * offv);
  void fit_model();
  std::vector<double> get_coef();
  std::vector<double> get_stderr();
  size_t get_rank_X() { return rank; };
  double get_dispersion() { return psi; };
  
private:
  void compute_variance(gsl_vector *w);
  
}; // IRLS

#endif // _IRLS_H_
