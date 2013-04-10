/** \file covariates.hpp
 *
 *  `Covariates' is a class 
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

#ifndef QUANTGEN_COVARIATES_HPP
#define QUANTGEN_COVARIATES_HPP

#include <vector>
#include <map>
#include <string>
#include <iostream>

#include "quantgen/samples.hpp"
#include "quantgen/utils_io.hpp"
#include "quantgen/utils_math.hpp"

namespace quantgen {
  
  class Covariates {
  private:
    map<string,map<string,vector<double> > > subgroup2covariates_;
    
  public:
    size_t GetNbSubgroups(void) const { return subgroup2covariates_.size(); };
    bool HasCovariates(const string & subgroup) const;
    size_t GetNbCovariates(const string & subgroup) const;
    map<string,vector<double> >::const_iterator begin(
      const string & subgroup) const;
    map<string,vector<double> >::const_iterator end(
      const string & subgroup) const;
    void AddSubgroup(const string & subgroup,
		     const map<string,vector<string> > & covariate2values);
    double GetCovariate(const string & subgroup, const string & covariate,
			const size_t & idx) const;
  };
  
} // namespace quantgen

#endif // QUANTGEN_COVARIATES_HPP
