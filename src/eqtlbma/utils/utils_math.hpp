/** \file utils_math.hpp
 *
 *  `utils_math' is a set of functions useful in statistical genetics.
 *  Copyright (C) 2013 Timothee Flutre
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef UTILS_UTILS_MATH_HPP
#define UTILS_UTILS_MATH_HPP

#include <cstdlib>

#include <string>
#include <limits>
#include <iostream>
#include <vector>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

namespace utils {

  const double NaN = std::numeric_limits<double>::quiet_NaN();

  bool isNonZero(size_t i);

  bool isNonNpos(size_t i);

  bool isNan(double i);

  size_t getSeed(void);

  double round(double x);
  
  size_t sum_bool(const std::vector<bool> & vec);
  
  void qqnorm(double * ptData, const size_t n);

  double log10_weighted_sum(const double * vec, const size_t size);

  double log10_weighted_sum(const double * vec, const double * weights,
			    const size_t size);

  void FitSingleGeneWithSingleSnp(const gsl_matrix * X,
				  const gsl_vector * y,
				  double & pve,
				  double & sigmahat,
				  double & betahat_geno,
				  double & sebetahat_geno,
				  double & betapval_geno);

  double mygsl_vector_sum(const gsl_vector * vec);

  void mygsl_vector_pow(gsl_vector * vec, const double exponent);

  void mygsl_matrix_pow(gsl_matrix * mat, const double exponent);

  gsl_matrix * mygsl_matrix_diagalloc(const gsl_vector * vec, const double x);

  gsl_matrix * mygsl_matrix_diagalloc(const gsl_matrix * mat, const double x);

  void mygsl_linalg_pseudoinverse(const gsl_matrix * X, gsl_matrix * X_ps);

  gsl_vector * mygsl_vector_alloc(const gsl_vector * vec);

  gsl_matrix * mygsl_matrix_alloc(const gsl_matrix * src);

  void mygsl_linalg_invert(const gsl_matrix * A, gsl_matrix * A_inv);

  void CalcMleErrorCovariance(const gsl_matrix * Y, const gsl_matrix * X,
			      gsl_matrix * XtX, gsl_matrix * Sigma_hat);

  void print_matrix(const gsl_matrix * A, const size_t M, const size_t N);

  void mygsl_linalg_outer(const gsl_vector * vec1, const gsl_vector * vec2,
			  gsl_matrix * mat);

} // namespace utils

#endif // UTILS_UTILS_MATH_HPP
