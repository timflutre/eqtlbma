/** \file data_loader.hpp
 *
 *  `data_loader' gathers functions useful to load gen-etics/-omics data
 *  Copyright (C) 2011-2014 Timothée Flutre
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

#ifndef QUANTGEN_DATA_LOADER_HPP
#define QUANTGEN_DATA_LOADER_HPP

#include <vector>
#include <map>
#include <set>
#include <iomanip>

#include "utils/utils_io.hpp"
#include "quantgen/samples.hpp"
#include "quantgen/gene.hpp"
#include "quantgen/snp.hpp"
#include "quantgen/covariates.hpp"

namespace quantgen {
  
  std::map<std::string,std::string> loadTwoColumnFile(
    const std::string & file,
    const int & verbose);

  void loadListsGenoExplevelAndCovarFiles(
    const std::string & file_genopaths,
    const std::string & file_exppaths,
    const std::string & file_covarpaths,
    const std::vector<std::string> & subgroups_tokeep,
    const std::string & error_model,
    const int & verbose,
    std::map<std::string,std::string> & subgroup2genofile,
    std::map<std::string,std::string> & subgroup2explevelfile,
    std::map<std::string,std::string> & subgroup2covarfile,
    std::vector<std::string> & subgroups);
  
  void loadSamplesFromMatrixEqtl(
    const std::map<std::string, std::string> & subgroup2file,
    const std::string & data_type,
    const int & verbose,
    std::vector<std::string> & samples,
    std::map<std::string,std::vector<std::string> > & subgroup2samples);
  
  void loadSamplesFromGenotypes(
    const std::map<std::string, std::string> & subgroup2genofile,
    const int & verbose,
    std::vector<std::string> & samples,
    std::map<std::string,std::vector<std::string> > & subgroup2samples_genotypes);
  
  void loadSamples(
    const std::map<std::string,std::string> & subgroup2genofile,
    const std::map<std::string,std::string> & subgroup2explevelfile,
    const std::map<std::string,std::string> & subgroup2covarfile,
    const std::string & error_model,
    const int & verbose,
    Samples & samples);
  
  void loadGeneInfo(
    const std::string & file_genecoords, const int & verbose,
    std::map<std::string,Gene> & gene2object,
    std::map<std::string,std::vector<Gene*> > & mChr2VecPtGenes);
  
  void loadExplevels(
    const std::map<std::string,std::string> & subgroup2explevelfile,
    const int & verbose,
    std::map<std::string,Gene> & gene2object);
  
  void loadSnpsToKeep(
    const std::string & file_snpstokeep,
    const int & verbose,
    std::set<std::string> & sSnpsToKeep);
  
  void loadGenosAndSnpInfoFromImpute(
    const std::map<std::string,std::string>::const_iterator & it_subgroup2genofile,
    const std::set<std::string> & sSnpsToKeep,
    const std::map<std::string, std::vector<Gene*> > & mChr2VecPtGenes,
    gzFile & genoStream,
    std::string & line,
    size_t & nb_lines,
    size_t & nb_snps_tokeep_per_subgroup,
    std::map<std::string, Snp> & snp2object,
    std::map<std::string, std::vector<Snp*> > & mChr2VecPtSnps);
  
  void loadGenosAndSnpInfoFromVcf (
    const std::map<std::string,std::string>::const_iterator & it_subgroup2genofile,
    const std::set<std::string> & sSnpsToKeep,
    const std::map<std::string, std::vector<Gene*> > & mChr2VecPtGenes,
    gzFile & genoStream,
    std::string & line,
    size_t & nb_lines,
    size_t & nb_snps_tokeep_per_subgroup,
    std::map<std::string, Snp> & snp2object,
    std::map<std::string, std::vector<Snp*> > & mChr2VecPtSnps);
  
  void duplicateGenosPerSnpInAllSubgroups(
    const std::map<std::string,std::string> & subgroup2genofile,
    const int & verbose,
    std::map<std::string, Snp> & snp2object);
  
  void loadGenosAndSnpInfo(
    const std::map<std::string, std::string> & subgroup2genofile,
    const float & min_maf,
    const std::set<std::string> & sSnpsToKeep,
    const std::map<std::string, std::vector<Gene*> > & mChr2VecPtGenes,
    const int & verbose,
    std::map<std::string, Snp> & snp2object,
    std::map<std::string, std::vector<Snp*> > & mChr2VecPtSnps);
  
  // void loadSnpInfo(const std::string & file_snpcoords,
  // 		   const std::string & file_snpcoords_idx,
  // 		   const std::set<std::string> & sSnpsToKeep,
  // 		   const std::map<std::string,Gene> gene2object,
  // 		   const std::string & anchor,
  // 		   const size_t & radius,
  // 		   const int & verbose,
  // 		   std::map<std::string, Snp> & snp2object,
  // 		   std::map<std::string, std::vector<Snp*> > & mChr2VecPtSnps);
  
  void loadSnpInfo(const std::string & snpCoordsFile,
		   const std::set<std::string> & sSnpsToKeep,
		   const int & verbose,
		   std::map<std::string, Snp> & snp2object,
		   std::map<std::string, std::vector<Snp*> > & mChr2VecPtSnps);
  
  // void loadSnpInfo(const std::string & file_snpcoords,
  // 		   const std::set<std::string> & sSnpsToKeep,
  // 		   const std::map<std::string,Gene> gene2object,
  // 		   const std::string & anchor,
  // 		   const size_t & radius,
  // 		   const int & verbose,
  // 		   std::map<std::string, Snp> & snp2object,
  // 		   std::map<std::string, std::vector<Snp*> > & mChr2VecPtSnps);
  
  void loadGenos(
    const std::map<std::string, std::string> & subgroup2genofile,
    const float & min_maf,
    const int & verbose,
    std::map<std::string, Snp> & snp2object);
  
  void loadListCovarFiles(
    const std::string & file_covarpaths,
    const std::string & sbgrpToKeep,
    const std::vector<std::string> & subgroups,
    std::map<std::string, std::string> & mCovarPaths,
    const int & verbose);
  
  void loadListCovarFiles(
    const std::string & file_covarpaths,
    const std::string & sbgrpToKeep,
    const std::vector<std::string> & subgroups,
    std::map<std::string, std::vector<std::string> > & mCovarPaths,
    const int & verbose);
  
  void loadCovariates(
    const std::map<std::string,std::string> subgroup2covarfile,
    const int & verbose,
    Covariates & covariates);
  
  // void loadRawInputData(
  //   const std::string & file_genopaths,
  //   const std::string & file_snpcoords,
  //   const std::string & file_exppaths,
  //   const std::string & file_genecoords,
  //   const std::string & anchor,
  //   const size_t & radius,
  //   const float & min_maf,
  //   const std::string & file_covarpaths,
  //   const std::string & error_model,
  //   const std::vector<std::string> & subgroups_tokeep,
  //   const std::set<std::string> & sSnpsToKeep,
  //   const int & verbose,
  //   std::vector<std::string> & subgroups,
  //   Samples & samples,
  //   std::map<std::string,Snp> & snp2object,
  //   std::map<std::string, std::vector<Snp*> > & mChr2VecPtSnps,
  //   Covariates & covariates,
  //   std::map<std::string, Gene> & gene2object);
  
  void loadListSstatsFile(
    const std::string & file_sstats,
    const int & verbose,
    std::map<std::string, std::string> & subgroup2sstatsfile);
  
  void fillGeneSnpPairsWithSstats(
    const std::map<std::string, std::string> & subgroup2sstatsfile,
    const int & verbose,
    std::map<std::string, Gene> & gene2object,
    std::map<std::string, Snp> & snp2object);
  
  void loadSummaryStats(
    const std::string & file_sstats,
    const int & verbose,
    std::vector<std::string> & subgroups,
    std::map<std::string, Gene> & gene2object,
    std::map<std::string, Snp> & snp2object);
  
} // namespace quantgen

#endif // QUANTGEN_DATA_LOADER_HPP
