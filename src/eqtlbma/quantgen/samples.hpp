/** \file samples.hpp
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

#ifndef QUANTGEN_SAMPLES_HPP
#define QUANTGEN_SAMPLES_HPP

#include <vector>
#include <map>
#include <string>
#include <iostream>

// #include "quantgen/gene.hpp" // do NOT uncomment, otherwise circular dependency

namespace quantgen {
  
  class Gene; // forward declaration
  
  class Samples {
  private:
    std::vector<std::string> all_;
    std::map<std::string,std::vector<bool> > subgroup2present_;
    std::map<std::string,std::vector<size_t> > subgroup2genotypes_;
    std::map<std::string,std::vector<size_t> > subgroup2explevels_;
    std::map<std::string,std::vector<size_t> > subgroup2covariates_;
    const std::vector<size_t> MapAllSamplesToTheGivenSubgroup(
      const std::vector<std::string> * pt_samples) const;
    void AddSubgroup(const std::string & subgroup,
		     const std::vector<size_t> & indices_of_all_in_subgroup);
  public:
    Samples(void);
    std::string GetSample(const size_t & idx) const;
    size_t GetTotalNbSamples(void) const { return all_.size(); };
    size_t GetTotalNbSubgroups(void) const { return subgroup2present_.size(); };
    bool IsPresent(const std::string & sample) const;
    bool IsAbsent(const std::string & sample) const;
    void AddSamplesIfNew(const std::vector<std::string> & samples);
    void AddSamplesFromData(
      const std::map<std::string,std::vector<std::string> > & subgroup2samples,
      const std::string & type_data);
    std::vector<std::string>::const_iterator begin(void) const { return all_.begin(); };
    std::vector<std::string>::const_iterator end(void) const { return all_.end(); };
    size_t GetIndexExplevel(const size_t & idx, const std::string & subgroup) const;
    size_t GetIndexGenotype(const size_t & idx, const std::string & subgroup) const;
    size_t GetIndexCovariate(const size_t & idx, const std::string & subgroup) const;
    void GetCommonAndUniqueIndividualsBetweenPairOfSubgroups(
      const std::string & subgroup1,
      const std::string & subgroup2,
      const Gene & gene,
      std::vector<size_t> & inds_s1s2,
      std::vector<size_t> & inds_s1,
      std::vector<size_t> & inds_s2) const;
    void ShowPairs(std::ostream & os) const;
    void ShowAllMappings(std::ostream & os) const;
  };

} // namespace quantgen

#endif // QUANTGEN_SAMPLES_HPP
