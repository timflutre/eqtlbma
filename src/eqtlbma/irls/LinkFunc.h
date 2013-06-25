/** \file LinkFunc.h
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

#ifndef _LINKFUNC_H_
#define _LINKFUNC_H_

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

class LinkFunc
{
public:
  bool quasi;
  virtual void init_mv(gsl_vector * Y, gsl_vector * mv)=0;
  virtual void compute_z(gsl_vector * y, gsl_vector * mv,
			 gsl_vector * offset, gsl_vector * z)=0;
  virtual void compute_weights(gsl_vector * mv, gsl_vector * w)=0;
  virtual void compute_mv(gsl_vector * bv, gsl_matrix * Xv,
			  gsl_vector * offset, gsl_vector * mv)=0;
  virtual double compute_dispersion(gsl_vector * y, gsl_matrix * Xv,
				    gsl_vector * bv, gsl_vector * offset,
				    gsl_vector * mv, double rank,
				    bool quasi_lik)=0;
  virtual ~LinkFunc() {};
  
}; // LinkFunc

#endif // _LINKFUNC_H_
