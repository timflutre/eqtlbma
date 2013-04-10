/** \file samples.cpp
 *
 *  `Samples' is a class 
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

using namespace std;

#include "quantgen/samples.hpp"

namespace quantgen {

  Samples::Samples(void)
  {
  }

  bool Samples::IsPresent(const string & sample) const
  {
    return find(all_.begin(), all_.end(), sample) != all_.end();
  }


  bool Samples::IsAbsent(const string & sample) const
  {
    return ! IsPresent(sample);
  }

  void Samples::AddSamplesIfNew(const vector<string> & samples)
  {
    for (vector<string>::const_iterator it = samples.begin();
	 it != samples.end(); ++it) {
      if (find(all_.begin(), all_.end(), *it) == all_.end())
	all_.push_back(*it);
    }
  }

/** \brief subgroup2genotypes_["s2"][0] = 5 means that the 1st sample in all_
 *  corresponds to the 6th genotype sample in subgroup "s2"
 *  \note 'npos' means that the sample is absent from the given subgroup
 */
  const vector<size_t> Samples::MapAllSamplesToTheGivenSubgroup(
    const vector<string> * pt_samples) const
  {
    vector<size_t> indices_of_all_in_subgroup(all_.size(), string::npos);
    vector<string>::const_iterator it_s; // points to a sample in a subgroup
  
    // for each sample in 'all_', record its index in 'pt_samples',
    // which corresponds to the proper column from the input file
    for (vector<string>::const_iterator it_a = all_.begin();
	 it_a != all_.end(); ++it_a) {
      it_s = find(pt_samples->begin(), pt_samples->end(), *it_a);
      if (it_s != pt_samples->end())
	indices_of_all_in_subgroup[it_a - all_.begin()] =
	  it_s - pt_samples->begin();
    }
  
    return indices_of_all_in_subgroup;
  }

  void Samples::AddSamplesFromGenotypes(
    const map<string,vector<string> > & subgroup2samples)
  {
    const vector<string> * pt_samples;
    for (map<string,vector<string> >::const_iterator it = 
	   subgroup2samples.begin(); it != subgroup2samples.end(); ++it) {
      pt_samples = &(it->second);
      subgroup2genotypes_.insert(
	make_pair(it->first, MapAllSamplesToTheGivenSubgroup(pt_samples)));
    }
  }

  void Samples::AddSamplesFromExplevels(
    const map<string,vector<string> > & subgroup2samples)
  {
    const vector<string> * pt_samples;
    for (map<string,vector<string> >::const_iterator it = 
	   subgroup2samples.begin(); it != subgroup2samples.end(); ++it) {
      pt_samples = &(it->second);
      subgroup2explevels_.insert(
	make_pair(it->first, MapAllSamplesToTheGivenSubgroup(pt_samples)));
    }
  }

  void Samples::AddSamplesFromCovariates(
    const map<string,vector<string> > & subgroup2samples)
  {
    const vector<string> * pt_samples;
    for (map<string,vector<string> >::const_iterator it = 
	   subgroup2samples.begin(); it != subgroup2samples.end(); ++it) {
      pt_samples = &(it->second);
      subgroup2covariates_.insert(
	make_pair(it->first, MapAllSamplesToTheGivenSubgroup(pt_samples)));
    }
  }

  size_t Samples::GetIndexExplevel(const size_t & idx, const string & subgroup) const
  {
    return subgroup2explevels_.find(subgroup)->second[idx];
  }

  size_t Samples::GetIndexGenotype(const size_t & idx, const string & subgroup) const
  {
    return subgroup2genotypes_.find(subgroup)->second[idx];
  }

  size_t Samples::GetIndexCovariate(const size_t & idx, const string & subgroup) const
  {
    return subgroup2covariates_.find(subgroup)->second[idx];
  }

  void Samples::GetCommonAndUniqueIndividualsBetweenPairOfSubgroups(
    const string & subgroup1,
    const string & subgroup2,
    vector<size_t> & inds_s1s2,
    vector<size_t> & inds_s1,
    vector<size_t> & inds_s2) const
  {
    bool present_in_s1, present_in_s2;
    for (size_t idx_all = 0; idx_all < all_.size(); ++idx_all) {
      present_in_s1 = (
	subgroup2genotypes_.find(subgroup1)->second[idx_all] != string::npos
	&& subgroup2explevels_.find(subgroup1)->second[idx_all] != string::npos);
      present_in_s2 = (
	subgroup2genotypes_.find(subgroup2)->second[idx_all] != string::npos
	&& subgroup2explevels_.find(subgroup2)->second[idx_all] != string::npos);
      if (present_in_s1 && present_in_s2)
	inds_s1s2.push_back(idx_all);
      else if (present_in_s1 && ! present_in_s2)
	inds_s1.push_back(idx_all);
      else if (! present_in_s1 && present_in_s2)
	inds_s2.push_back(idx_all);
    }
  }

} // namespace quantgen
