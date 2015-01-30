////////////////////////////////////////////////////////////////
//
// Copyright (C) 2008 Affymetrix, Inc.
//
// This program is free software; you can redistribute it and/or modify 
// it under the terms of the GNU General Public License (version 2) as 
// published by the Free Software Foundation.
// 
// This program is distributed in the hope that it will be useful, 
// but WITHOUT ANY WARRANTY; without even the implied warranty of 
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
// General Public License for more details.
// 
// You should have received a copy of the GNU General Public License 
// along with this program;if not, write to the 
// 
// Free Software Foundation, Inc., 
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//
////////////////////////////////////////////////////////////////

#ifndef _QUANTLABELZIO_H_
#define _QUANTLABELZIO_H_

//
#include "chipstream/QuantLabelZ.h"
#include "chipstream/TsvReport.h"
//
#include "copynumber/CNLog2RatioData.h"
//
#include <cstring>
#include <string>

//
#define QUANTLABELZ__SNPTSV_DATAORDER "meanX,varX,nObsMean,nObsVar,meanY,varY,covarXY"

//
void QuantLabelZ__readSnpPriorMap(AffxArray<snp_labeled_distribution>& Priors, 
                                  const std::string& fileName);
//
void QuantLabelZ__readSnpPriorMap_v1(AffxArray<snp_labeled_distribution>& Priors, 
                                     const std::string& fileName);
void QuantLabelZ__readSnpPriorMap_v2(AffxArray<snp_labeled_distribution>& Priors, 
                                     const std::string& fileName);
//
void QuantLabelZ__readSnpPriorMap_tsv5(AffxArray<snp_labeled_distribution>& Priors,
                                       affx::File5_Tsv* tsv5);
void QuantLabelZ__readSnpPriorMap_tsv5_v1(AffxArray<snp_labeled_distribution>& Priors,
                                          affx::File5_Tsv* tsv5);
void QuantLabelZ__readSnpPriorMap_tsv5_v2(AffxArray<snp_labeled_distribution>& Priors,
                                          affx::File5_Tsv* tsv5);
void QuantLabelZ__readSnpPriorMap_tsv5_v2(int iXChromosome, int iYChromosome, CNExperiment& objExperiment, CNProbeSetArray& vProbeSets,
                                          affx::File5_Tsv* tsv5);

void QuantLabelZ__readSnpPrior_tsv5_v2(affx::File5_Tsv*tsv5, snp_distribution &cluster);
//
void QuantLabelZ__writeSnpPosteriorTsv_defC_version_2(affx::TsvReport& tsv,const std::string& prefix,int* cidx);
void QuantLabelZ__writeSnpPosteriorTsv_defD_version_2(affx::TsvReport& tsv,const std::string& prefix,int* cidx);

void QuantLabelZ__writeSnpPosteriorTsv_defC_version_3(affx::TsvReport& tsv,const std::string& prefix,int* cidx);
void QuantLabelZ__writeSnpPosteriorTsv_defD_version_3(affx::TsvReport& tsv,const std::string& prefix,int* cidx);

void QuantLabelZ__writeSnpPosteriorTsv(affx::TsvReport& tsv,
                                       const std::string& fileName,
                                       affx::TsvReport::TsvReportFmt_t format_file,
                                       int format_ver);
void QuantLabelZ__writeSnpPosteriorTsv(affx::TsvReport& tsv,
                                       const std::string& fileName,
                                       affx::TsvReport::TsvReportFmt_t format_file,
                                       int format_ver);

void QuantLabelZ__writeSnpPosteriorValue_v2_C(affx::TsvReport& tsv,
                                              const std::string& prefix,
                                              cluster_data& c);
void QuantLabelZ__writeSnpPosteriorValue_v2_D(affx::TsvReport& tsv,
                                              const std::string& prefix,
                                              snp_distribution& d);
void QuantLabelZ__writeSnpPosteriorValue_v3_C(affx::TsvReport& tsv,
                                              const std::string& prefix,
                                              cluster_data& c);
void QuantLabelZ__writeSnpPosteriorValue_v3_D(affx::TsvReport& tsv,
                                              const std::string& prefix,
                                              snp_distribution& d);
void QuantLabelZ__writeSnpPosteriorValue(affx::TsvReport& tsv,
                                         int format_ver,
                                         const std::string& probeset_id,
                                         snp_param &tsp);

std::string c_to_str(cluster_data d, const std::string& sep, int dimension);
void str_to_c(cluster_data &d,const std::string& str);
void str_to_cv(snp_distribution& Dist,const std::string& str);

void QuantLabelZ__SnpPriorFromStrings(snp_labeled_distribution &tsp,
                                      const std::string& bb, 
                                      const std::string& ab,
                                      const std::string& aa,
                                      const std::string& cv);

snp_labeled_distribution 
QuantLabelZ__CreatePriorFromStrings(const std::string& bb,
                                    const std::string& ab,
                                    const std::string& aa,
                                    const std::string& cv);
//
void 
QuantLabelZ__convert_priors_tsv(const std::string& file_in,
                                const std::string& file_out,
                                int format_ver,
                                int verbose);
void
QuantLabelZ__convert_priors_tsv5(const std::string& file_in,
                                 const std::string& tsv_name,
                                 const std::string& file_out,
                                 int format_ver,
                                 int verbose);

void QuantLabelZ__convert_priors_tsv_to_tsv5(const std::string& file_in,
                                             const std::string& tsv_name,
                                             const std::string& file_out,
                                             int format_ver,
                                             int verbose);

#endif // _QUANTLABELZIO_H_
