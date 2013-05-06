/** \file utils_math.cpp
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

#include <cmath>
#include <sys/time.h>

#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include "utils/utils_math.hpp"

using namespace std;

namespace utils {

  bool isNonZero(size_t i) { return (i != 0); };

  bool isNonNpos(size_t i) { return (i != string::npos); };

// http://www.johndcook.com/IEEE_exceptions_in_cpp.html
  bool isNan(double i) { return (! (i == i)); };

/** \brief Return a seed based on microseconds since epoch.
 *  \note http://www.guyrutenberg.com/2007/09/03/seeding-srand/
 */
  size_t getSeed(void)
  {
    timeval t1;
    gettimeofday (&t1, NULL);
    return ((size_t) t1.tv_usec * t1.tv_sec);
  }

  size_t sum_bool(const vector<bool> & vec)
  {
    size_t res = 0;
    for(vector<bool>::const_iterator it = vec.begin(); it != vec.end(); ++it)
      if(*it)
	++res;
    return res;
  }

/** \brief Round the given value.
 *  \note http://stackoverflow.com/a/485549/597069
 */
  double round(double x)
  {
    return (x > 0.0) ? floor(x + 0.5) : ceil(x - 0.5);
  }

/** \brief Quantile-normalize an input vector to a standard normal.
 *  \note Missing values should be removed beforehand.
 *  \note code inspired from "qqnorm" in GNU R.
 */
  void qqnorm(double * ptData, const size_t n)
  {
    size_t * order = (size_t*) calloc(n, sizeof(size_t));
    if (order == NULL) {
      fprintf(stderr, "ERROR: can't allocate memory for order in qqnorm\n");;
      exit(1);
    }
    gsl_sort_index(order, ptData, 1, n);
  
    double q, a = (n <= 10 ? 0.375 : 0.5);
    for (size_t i=0; i<n; ++i) {
      q = (i+1 - a) / (n + 1 - 2 * a);
      ptData[order[i]] = gsl_cdf_ugaussian_Pinv(q);
    }
  
    free(order);
  }

/** \brief Return log_{10}(\sum_{i=1}^n 1/n 10^vec_i)
 */
  double log10_weighted_sum(const double * vec, const size_t size)
  {
    size_t i = 0;
    double res, max = vec[0];
    double * weights = (double*) calloc(size, sizeof(double));
    if (weights == NULL) {
      fprintf(stderr, "ERROR: can't allocate memory for weights\n");
      exit(1);
    }
    for (i = 0; i < size; ++i) {
      if (vec[i] > max)
	max = vec[i];
      weights[i] = (double) (1 / ((double) size));
    }
    double sum = 0.0;
    for (i = 0; i < size; ++i)
      sum += weights[i] * pow(10, vec[i] - max);
    free(weights);
    res = max + log10(sum);
    if (abs(res) <= GSL_DBL_EPSILON)
      res = 0.0;
    return res;
  }

/** \brief Return log_{10}(\sum_i w_i 10^vec_i)
 */
  double log10_weighted_sum(const double * vec, const double * weights,
			    const size_t size)
  {
    size_t i = 0;
    double res, max = vec[0];
    for (i = 0; i < size; ++i)
      if (vec[i] > max)
	max = vec[i];
    double sum = 0.0;
    for (i = 0; i < size; ++i)
      sum += weights[i] * pow(10, vec[i] - max);
    res = max + log10(sum);
    if (abs(res) <= GSL_DBL_EPSILON)
      res = 0.0;
    return res;
  }

/** \brief Estimate by ML the effect size of the genotype, the std deviation 
 *  of the errors and the std error of the estimated effect size in the 
 *  multiple linear regression Y = XB + E with E~MVN(0,sigma^2I)
 *  \note genotype supposed to be 2nd column of X
 */
  void FitSingleGeneWithSingleSnp(const gsl_matrix * X,
				  const gsl_vector * y,
				  double & pve,
				  double & sigmahat,
				  double & betahat_geno,
				  double & sebetahat_geno,
				  double & betapval_geno)
  {
    size_t N = X->size1, P = X->size2, rank;
    double rss;
    gsl_vector * Bhat = gsl_vector_alloc(P);
    gsl_matrix * covBhat = gsl_matrix_alloc(P, P);
    gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc(N, P);
    gsl_multifit_linear_svd(X, y, GSL_DBL_EPSILON, &rank, Bhat, covBhat,
			    &rss, work);
    pve = 1 - rss / gsl_stats_tss(y->data, y->stride, y->size);
    sigmahat = sqrt(rss / (double)(N-rank));
    betahat_geno = gsl_vector_get(Bhat, 1);
    sebetahat_geno = sqrt(gsl_matrix_get (covBhat, 1, 1));
    betapval_geno = 2 * gsl_cdf_tdist_Q(fabs(betahat_geno / sebetahat_geno),
					N-rank);
    gsl_vector_free(Bhat);
    gsl_matrix_free(covBhat);
    gsl_multifit_linear_free(work);
  }

  double mygsl_vector_sum(const gsl_vector * vec)
  {
    double res = 0.0;
    for(size_t i = 0; i < vec->size; ++i)
      res += gsl_vector_get(vec, i);
    return res;
  }

  void mygsl_vector_pow(gsl_vector * vec, const double exponent)
  {
    for(size_t i = 0; i < vec->size; ++i)
      gsl_vector_set(vec, i, pow(gsl_vector_get(vec, i), exponent));
  }

  void mygsl_matrix_pow(gsl_matrix * mat, const double exponent)
  {
    for(size_t i = 0; i < mat->size1; ++i)
      for(size_t j = 0; j < mat->size2; ++j)
	gsl_matrix_set(mat, i, j, pow(gsl_matrix_get(mat, i, j), exponent));
  }

// from http://lists.gnu.org/archive/html/help-gsl/2005-09/msg00007.html
  gsl_matrix * mygsl_matrix_diagalloc(const gsl_vector * vec, const double x)
  {
    gsl_matrix * mat = gsl_matrix_alloc(vec->size, vec->size);
    gsl_matrix_set_all(mat, x);
    gsl_vector_view diag = gsl_matrix_diagonal(mat);
    gsl_vector_memcpy(&diag.vector, vec);
    return mat;
  }

  gsl_matrix * mygsl_matrix_diagalloc(const gsl_matrix * mat, const double x)
  {
    gsl_vector_const_view diag = gsl_matrix_const_diagonal(mat);
    return mygsl_matrix_diagalloc(&diag.vector, x);
  }

// if A is M x N, then A_ps is NxM
  void mygsl_linalg_pseudoinverse(const gsl_matrix * A, gsl_matrix * A_ps)
  {
    size_t M = A->size1, N = A->size2;
  
    // A = U D V' where U is MxN and V is NxN
    gsl_matrix * U = gsl_matrix_alloc(M, N), * V = gsl_matrix_alloc(N, N);
    gsl_vector * D_diag = gsl_vector_alloc(N),
      * work = gsl_vector_alloc(N);
    gsl_matrix_memcpy(U, A);
    gsl_linalg_SV_decomp(U, V, D_diag, work);
  
    // V D^(-1)
    gsl_matrix * VDinv = gsl_matrix_alloc(N, N);
    mygsl_vector_pow(D_diag, -1.0);
    gsl_matrix * D_inv = mygsl_matrix_diagalloc(D_diag, 0.0);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, V, D_inv, 0.0, VDinv);
  
    // A_ps = VDinv U'
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, VDinv, U, 0.0, A_ps);
  
    gsl_matrix_free(U);
    gsl_matrix_free(V);
    gsl_vector_free(D_diag);
    gsl_vector_free(work);
    gsl_matrix_free(VDinv);
    gsl_matrix_free(D_inv);
  }

  gsl_vector * mygsl_vector_alloc(const gsl_vector * src)
  {
    gsl_vector * dst = gsl_vector_alloc(src->size);
    gsl_vector_memcpy(dst, src);
    return dst;
  }

  gsl_matrix * mygsl_matrix_alloc(const gsl_matrix * src)
  {
    gsl_matrix * dst = gsl_matrix_alloc(src->size1, src->size2);
    gsl_matrix_memcpy(dst, src);
    return dst;
  }

  void mygsl_linalg_invert(const gsl_matrix * A, gsl_matrix * A_inv)
  {
    gsl_matrix * tmp = gsl_matrix_alloc(A->size1, A->size2);
    gsl_matrix_memcpy(tmp, A);
    gsl_permutation * perm = gsl_permutation_alloc(A->size1);
    int signum;
    gsl_linalg_LU_decomp(tmp, perm, &signum);
    gsl_linalg_LU_invert(tmp, perm, A_inv);
    gsl_matrix_free(tmp);
    gsl_permutation_free(perm);
  }

/** \brief Estimate by ML the covariance matrix Sigma of the errors in
 *  the multivariate linear regression Y = XB + E with E~MN(0,I,Sigma)
 */
  void CalcMleErrorCovariance(const gsl_matrix * Y, const gsl_matrix * X,
			      gsl_matrix * XtX, gsl_matrix * Sigma_hat)
  {
    size_t N = X->size1, P = X->size2;
  
    bool calc_XtX = false;
    if(XtX == NULL){
      calc_XtX = true;
      XtX = gsl_matrix_alloc(P, P);
      gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, X, X, 0.0, XtX);
    }
    gsl_matrix * XtX_inv = gsl_matrix_alloc(P, P);
    mygsl_linalg_pseudoinverse(XtX, XtX_inv);
  
    gsl_matrix * tmp1 = gsl_matrix_alloc(N, P); // X (X'X)^-1
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, X, XtX_inv, 0.0, tmp1);
  
    gsl_matrix * tmp2 = gsl_matrix_alloc(N, N); // X (X'X)^-1 X'
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, tmp1, X, 0.0, tmp2);
  
    gsl_matrix * T = gsl_matrix_alloc(N, N); // I - X (X'X)^-1 X'
    for(size_t i = 0; i < T->size1; ++i){
      for(size_t j = 0; j < T->size2; ++j){
	double T_ij = - gsl_matrix_get(tmp2, i, j);
	if(i == j)
	  T_ij += 1;
	gsl_matrix_set(T, i, j, T_ij);
      }
    }
  
    gsl_matrix * tmp3 = gsl_matrix_alloc(N, 2); // (I - X (X'X)^-1 X') Y
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T, Y, 0.0, tmp3);
  
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1/(double)N, Y, tmp3, 0.0, Sigma_hat);
  
    if(calc_XtX)
      gsl_matrix_free(XtX);
    gsl_matrix_free(XtX_inv);
    gsl_matrix_free(tmp1);
    gsl_matrix_free(tmp2);
    gsl_matrix_free(T);
    gsl_matrix_free(tmp3);
  }

  void print_matrix(const gsl_matrix * A, const size_t M, const size_t N)
  {
    for(size_t i = 0; i < min(M,A->size1); ++i){
      for(size_t j = 0; j < min(N,A->size2); ++j)
	fprintf(stderr, "%e  ", gsl_matrix_get(A,i,j));
      fprintf(stderr, "\n");
    }
  }

  double mygsl_linalg_det(const gsl_matrix * A)
  {
    double det = NaN;
    gsl_matrix * tmp = gsl_matrix_alloc(A->size1, A->size2);
    gsl_matrix_memcpy(tmp, A);
    gsl_permutation * perm = gsl_permutation_alloc(A->size1);
    int signum;
    gsl_linalg_LU_decomp(tmp, perm, &signum);
    gsl_linalg_LU_lndet(tmp);
    gsl_matrix_free(tmp);
    gsl_permutation_free(perm);
    return det;
  }

/** \brief Fill matrix with the outer product of vec1 and vec2
 *  \note mat = vec1 vec2^T
 */
  void mygsl_linalg_outer(const gsl_vector * vec1, const gsl_vector * vec2,
			  gsl_matrix * mat)
  {
    for(size_t i = 0; i < mat->size1; ++i)
      for(size_t j = 0; j < mat->size2; ++j)
	gsl_matrix_set(mat, i, j, gsl_vector_get(vec1, i) *
		       gsl_vector_get(vec2, j));
  }

} // namespace utils
