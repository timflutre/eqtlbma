/** \file snp.hpp
 *
 *  `Snp' is a class 
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

#ifndef QUANTGEN_SNP_HPP
#define QUANTGEN_SNP_HPP

#include <vector>
#include <map>
#include <string>
#include <iostream>

#include "quantgen/utils_io.hpp"
#include "quantgen/utils_math.hpp"

namespace quantgen {

  class Snp {
  private:
    string name_;
    string chromosome_;
    size_t pos_; // 1-based
  
    map<string,vector<double> > subgroup2genotypes_;
    map<string,double> subgroup2maf_; // minor allele frequency
  
  public:
    Snp(void);
    Snp(const string & name, const string & chr, const string & coord);
    string GetName(void) const { return name_; };
    string GetChromosome(void) const { return chromosome_; };
    size_t GetPosition(void) const { return pos_; };
    size_t GetNbSubgroups(void) const { return subgroup2genotypes_.size(); };
    bool HasGenotypesInAtLeastOneSubgroup(void) const;
    bool HasGenotypes(const string & subgroup) const;
    bool HasGenotypesInAllSubgroups(const vector<string> & subgroups) const;
    double GetMinorAllelFreq(const string & subgroup) const;
    void Show(ostream & os);
    void AddSubgroupFromImputeLine(const string & subgroup,
				   const vector<string>::const_iterator & begin,
				   const vector<string>::const_iterator & end);
    void AddSubgroupFromVcfLine(const string & subgroup,
				const vector<string>::const_iterator & begin,
				const vector<string>::const_iterator & end);
    void AddSubgroupFromDoseLine(const string & subgroup,
				 const vector<string>::const_iterator & begin,
				 const vector<string>::const_iterator & end);
    void AddSubgroup(const string & subgroup,
		     const vector<string>::const_iterator & begin,
		     const vector<string>::const_iterator & end,
		     const string & format);
    void EraseIfLowMafPerSubgroup(const double & min_maf);
    void DuplicateGenotypesFromFirstSubgroup(const string & subgroup_old,
					     const string & subgroup_new);
    int IsInCis(const size_t & start, const size_t & end, const string & anchor,
		const size_t & radius) const;
    double GetGenotype(const string & subgroup, const size_t & idx) const;
    friend class Gene;
  };

  bool operator==(const Snp& lhs, const Snp& rhs);
  bool operator!=(const Snp& lhs, const Snp& rhs);
  bool operator< (const Snp& lhs, const Snp& rhs);
  bool operator> (const Snp& lhs, const Snp& rhs);
  bool operator<=(const Snp& lhs, const Snp& rhs);
  bool operator>=(const Snp& lhs, const Snp& rhs);

  bool pt_snp_lt_pt_snp(const Snp* pt_lhs, const Snp* pt_rhs); // lt means <

} // namespace quantgen

#endif // QUANTGEN_SNP_HPP
