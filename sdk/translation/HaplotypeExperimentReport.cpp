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
/**
 * @file   HaplotypeExperimentReport.cpp
 * @author Mybrid Spalding
 * @date   Tue Mar 18 09:02:01 PDT 2008
 * @brief  DEPRECATED DMET2 Class for generating the haplotype report as a report file. 
 */

#include "translation/HaplotypeExperimentReport.h"
//
#include "util/Guid.h"

void HaplotypeExperimentReport::open(const std::string& fileName)
{
  // our init
  m_guid=affxutil::Guid::GenerateNewGuid();
  // describe our columns.
  m_tsv.defineFile("experiment"    "\t"
                   "geneId"        "\t"
                   "sample"        "\t"
                   "call"          "\t"
                   "callCount"     "\t"
                   "knownCount"    "\t"
                   "UNKExists"     "\t"
                   "basecallRate"  "\t"
                   "basecallCount" "\t"
                   "noCallCount"   "\t"
                   "possibleRareAlleleCount" "\t"
                   "notAvailableCount");
  //
  m_tsv.addHeader("dmet3-file-type","haplotype-experiment-report");
  m_tsv.addHeader("dmet3-report-guid",m_guid);
  //
  m_tsv.writeTsv_v1(fileName);
}

void HaplotypeExperimentReport::close()
{
  m_tsv.close();
}


void HaplotypeExperimentReport::report(const ExperimentGeneResults& result)
{
  // here we expect things to be in the order defined above.
  // otherwise we could say: m_tsv.set(0,"experiment",m_experiment);
  //  m_tsv.set(0, 0,result.m_experimentName   );
  //m_tsv.set(0, 1,result.m_geneName       );
  //  m_tsv.set(0, 2,result.m_sample       );
  //  m_tsv.set(0, 3,result.m_call         );
  //m_tsv.set(0, 4,result.m_callCount    );
  //m_tsv.set(0, 5,result.m_knownCount   );
  //m_tsv.set(0, 6,result.m_UNKExists    );
  // m_tsv.set(0, 7,result.m_basecallRate );
  //m_tsv.set(0, 8,result.m_basecallCount);
  //m_tsv.set(0, 9,result.m_noCallCount  );
  //m_tsv.set(0,10,result.m_possibleRareAlleleCount);
  //m_tsv.set(0,11,result.m_notAvailableCount);
  //
  m_tsv.writeLevel(0);
}

