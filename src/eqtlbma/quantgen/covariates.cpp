/** \file covariates.cpp
 *
 *  `Covariates' is a class 
 *  Copyright (C) 2013 Timothee Flutre
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Snpral Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Snpral Public License for more details.
 *
 *  You should have received a copy of the GNU Snpral Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "quantgen/covariates.hpp"

using namespace std;

using namespace utils;

namespace quantgen {
  
  bool Covariates::HasCovariates(const string & subgroup) const
  {
    return subgroup2covariates_.find(subgroup) != subgroup2covariates_.end();
  }
  
  size_t Covariates::GetNbCovariates(
    const string & subgroup) const
  {
    size_t res = 0;
    map<string,map<string,vector<double> > >::const_iterator it
      = subgroup2covariates_.find(subgroup);
    if (it != subgroup2covariates_.end())
      res = it->second.size();
    return res;
  }
  
  map<string,vector<double> >::const_iterator
  Covariates::begin(const string & subgroup) const
  {
    return subgroup2covariates_.find(subgroup)->second.begin();
  }
  
  map<string,vector<double> >::const_iterator
  Covariates::end(const string & subgroup) const
  {
    return subgroup2covariates_.find(subgroup)->second.end();
  }
  
  void Covariates::AddSubgroup(
    const string & subgroup,
    const map<string,vector<string> > & covariate2values)
  {
    subgroup2covariates_.insert(make_pair(subgroup,
					  map<string,vector<double> >()));
    for (map<string,vector<string> >::const_iterator itC
	   = covariate2values.begin(); itC != covariate2values.end(); ++itC) {
      vector<double> values(itC->second.size(), NaN);
      for (vector<string>::const_iterator itV = itC->second.begin();
	   itV != itC->second.end(); ++itV) {
	if (itV->compare("NA") == 0 || itV->compare("na") == 0
	    || itV->compare("NaN") == 0 || itV->compare("nan") == 0) {
	  cerr << "ERROR: no missing value allowed, see covariate " << itC->first
	       << " in subgroup " << subgroup << endl;
	  exit(EXIT_FAILURE);
	}
	values[itV - itC->second.begin()] = atof(itV->c_str());
      }
      subgroup2covariates_[subgroup].insert(make_pair(itC->first, values));
    }
  }
  
  double Covariates::GetCovariate(const string & subgroup,
				  const string & covariate,
				  const size_t & idx) const
  {
    return subgroup2covariates_.find(subgroup)->second.find(covariate)->second[idx];
  }
  
} // namespace quantgen
