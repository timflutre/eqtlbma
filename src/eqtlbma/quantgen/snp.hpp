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

#include "utils/utils_io.hpp"
#include "utils/utils_math.hpp"

namespace quantgen {
  
  struct Genotypes {
    std::vector<double> values;
    double maf; // minor allele frequency
  };
  
  class Snp {
  private:
    std::string name_;
    std::string chromosome_;
    size_t pos_; // 1-based
    
    std::map<std::string,Genotypes> subgroup2genotypes_;
    
  public:
    Snp(void);
    Snp(const std::string & name, const std::string & chr, const std::string & coord);
    std::string GetName(void) const { return name_; };
    std::string GetChromosome(void) const { return chromosome_; };
    size_t GetPosition(void) const { return pos_; };
    size_t GetNbSubgroups(void) const { return subgroup2genotypes_.size(); };
    size_t GetNbSamples(const std::string & subgroup) const;
    double GetMinorAlleleFreq(const std::string & subgroup) const;
    bool HasGenotypesInAtLeastOneSubgroup(void) const;
    bool HasGenotypes(const std::string & subgroup) const;
    bool HasGenotypesInAllSubgroups(const std::vector<std::string> & subgroups) const;
    void Show(std::ostream & os);
    void AddSubgroupFromImputeLine(const std::string & subgroup,
				   const std::vector<std::string>::const_iterator & begin,
				   const std::vector<std::string>::const_iterator & end,
				   std::vector<double> & genotypes,
				   double & minor_allele_freq);
    void AddSubgroupFromVcfLine(const std::string & subgroup,
				const std::vector<std::string>::const_iterator & begin,
				const std::vector<std::string>::const_iterator & end,
				std::vector<double> & genotypes,
				double & minor_allele_freq);
    void AddSubgroupFromDoseLine(const std::string & subgroup,
				 const std::vector<std::string>::const_iterator & begin,
				 const std::vector<std::string>::const_iterator & end,
				 std::vector<double> & genotypes,
				 double & minor_allele_freq);
    void AddSubgroup(const std::string & subgroup,
		     const std::vector<std::string>::const_iterator & begin,
		     const std::vector<std::string>::const_iterator & end,
		     const std::string & format);
    void EraseIfLowMafPerSubgroup(const double & min_maf);
    void DuplicateGenotypesFromFirstSubgroup(const std::string & subgroup_old,
					     const std::string & subgroup_new);
    int IsInCis(const size_t & start, const size_t & end, const std::string & anchor,
		const size_t & radius) const;
    double GetGenotype(const std::string & subgroup, const size_t & idx) const;
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
