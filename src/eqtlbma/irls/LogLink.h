/** \file LogLink.h
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

#ifndef _LOGLINK_H_
#define _LOGLINK_H_

#include "LinkFunc.h"



class LogLink : public LinkFunc{
  
  gsl_vector *init_mv(gsl_vector *Y);
  gsl_vector *compute_Z(gsl_vector *Y, gsl_vector *mv);
  gsl_vector *compute_weights(gsl_vector *mv);
  
  gsl_vector *compute_mv(gsl_vector *bv, gsl_matrix *Xv);
  
  bool quasi;
  double compute_dispersion(gsl_vector *Y, gsl_matrix *Xv, gsl_vector *bv, double rank, bool quasi_lik);
};

#endif
