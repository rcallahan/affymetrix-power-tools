////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
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
 * @file   QuantMethodExprReport.cpp
 * @author Chuck Sugnet
 * @date   Mon Oct 24 12:09:33 2005
 * 
 * @brief  Class for reporting results of quantification methods.
 */

//
#include "chipstream/QuantMethodExprReport.h"
//
#include "util/Fs.h"
#include "util/Util.h"
//
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <map>
#include <string>

QuantMethodExprReport::QuantMethodExprReport(int numCol) {
    init();
    m_NumCol = numCol; 
}

QuantMethodExprReport::QuantMethodExprReport() {
    init();
}

void QuantMethodExprReport::init() {
  
    m_NumCol = 0;
    m_ProbeSetsToReport = NULL;
    m_MaxNameLength = TSVREPORT_PROBESET_STRLEN;
    affx::TsvReport::init();
    // @todo: why cant we do this? -jhg
    // m_is_header_buffer=true;  
  
    m_summary_tsv_precision=5;
    m_feffects_tsv_precision=5;
    m_residuals_tsv_precision=5;

    m_DoSummary = true;
    m_DoFeatureEffects = true;
    m_DoResiduals = true;
    m_DoCompact = false;
    m_File5 = NULL;
    m_F5Group = NULL;
    m_summary_a5_group = NULL;
    m_feffects_a5_group = NULL;
    m_residuals_a5_group = NULL;
    m_WriteOldStyleFeatureEffectsFile = false;
    m_unmungePSName = false;
}

QuantMethodExprReport::~QuantMethodExprReport() {
    if (m_F5Group != NULL) {
        m_F5Group->close();
        Freez(m_F5Group);
    }
    if (m_File5 != NULL) {
        m_File5->close();
        Freez(m_File5);
    }
    affx::TsvReport::close();
}

void QuantMethodExprReport::setCompactFile5Format(int maxNameLength)  {
    std::string fileName = Fs::join(m_dir_path,m_file_prefix + ".summaries.a5");
    m_File5 = TsvReport::openA5File(fileName, true);
    m_group_name ="/";
    m_MaxNameLength = maxNameLength;
    m_F5Group = m_File5->openGroup(m_group_name,affx::FILE5_CREATE|affx::FILE5_OPEN);
    m_a5_shared_group = m_F5Group;
    m_summary_tsv.setA5SharedGroup(m_F5Group);
    m_feffects_tsv.setA5SharedGroup(m_F5Group);
    m_residuals_tsv.setA5SharedGroup(m_F5Group);
    m_file_prefix="/";
    m_format = TsvReport::FMT_A5;
    m_DoCompact = true;
    m_summary_a5_group = m_F5Group;
    m_feffects_a5_group = m_F5Group;
    m_residuals_a5_group = m_F5Group;
}

bool QuantMethodExprReport::prepare(QuantMethod &qMethod, const IntensityMart &iMart) {
    // printf("### QuantMethodExprReport: '%s'\n",m_file_name.c_str());
    std::string delim = ".";
    if (m_DoCompact) {
        delim = "";
    }

    if (m_DoSummary) {
        copyOptionsTo(m_summary_tsv);
        m_summary_tsv.setFilename(m_file_prefix + delim + "summary");
        m_summary_tsv.addHeadersFrom(m_header_buffer);
        m_summary_tsv.setPrecision(m_summary_tsv_precision);
        m_summary_tsv.setA5SharedGroup(m_summary_a5_group);
        m_summary_tsv.setFormat(m_format);
    }
    if (m_DoFeatureEffects) {
	copyOptionsTo(m_feffects_tsv);
	m_feffects_tsv.setFilename(m_file_prefix+ delim + "feature-response");
	m_feffects_tsv.addHeadersFrom(m_header_buffer);
	m_feffects_tsv.setPrecision(m_feffects_tsv_precision);
        m_feffects_tsv.setA5SharedGroup(m_feffects_a5_group);
	m_feffects_tsv.setFormat(m_format);
    }
    if (m_DoResiduals) {
        copyOptionsTo(m_residuals_tsv);
        m_residuals_tsv.setFilename(m_file_prefix+ delim + "residuals");
        m_residuals_tsv.addHeadersFrom(m_header_buffer);
        m_residuals_tsv.setPrecision(m_residuals_tsv_precision);
        m_residuals_tsv.setA5SharedGroup(m_residuals_a5_group);
        m_residuals_tsv.setFormat(m_format); 
    }

    // fix names
    std::vector<std::string> col_names=iMart.getCelFileNames();
    for (size_t i=0;i<col_names.size();i++) {
        col_names[i]=Fs::basename(col_names[i]);
    }

    // summary
    if (m_DoSummary) {
        m_summary_tsv.defineStringColumn(0,0,"probeset_id",m_MaxNameLength);
        if (m_DoCompact) {
            m_summary_tsv.defineColumns(col_names,affx::FILE5_DTYPE_FLOAT,0);
        }
        else {
            m_summary_tsv.defineColumns(col_names,affx::FILE5_DTYPE_DOUBLE,0);
        }
        //
        m_summary_tsv.writeTsv_v1();
    }

    // feffects => feature effects
    if (m_DoFeatureEffects) {
        if (m_WriteOldStyleFeatureEffectsFile) {
            if (m_feffects_tsv.m_format == affx::TsvReport::FMT_A5) {
                m_feffects_tsv.defineColumn(0,0,"probe_id",             affx::FILE5_DTYPE_INT);
                m_feffects_tsv.defineColumn(0,1,"feature_response",     affx::FILE5_DTYPE_DOUBLE);
                m_feffects_tsv.defineColumn(0,2,"probeset_id",          affx::FILE5_DTYPE_STRING, m_MaxNameLength);
            } 
            else {
                m_feffects_tsv.defineColumn(0,1,"atom_id"          ,affx::FILE5_DTYPE_INT);
                m_feffects_tsv.defineColumn(0,2,"probe_id",        affx::FILE5_DTYPE_INT);
                m_feffects_tsv.defineColumn(0,3,"x"                ,affx::FILE5_DTYPE_INT);
                m_feffects_tsv.defineColumn(0,4,"y"                ,affx::FILE5_DTYPE_INT);
                m_feffects_tsv.defineColumn(0,5,"feature_response",affx::FILE5_DTYPE_DOUBLE);
 
            }
        } 
        else { 
            if (m_feffects_tsv.m_format == affx::TsvReport::FMT_A5) {
                m_feffects_tsv.defineColumn(0,0,"probe_id",             affx::FILE5_DTYPE_INT);
                m_feffects_tsv.defineColumn(0,1,"feature_response",     affx::FILE5_DTYPE_DOUBLE);
                m_feffects_tsv.defineColumn(0,2,"probeset_id",          affx::FILE5_DTYPE_STRING, m_MaxNameLength);
                m_feffects_tsv.defineColumn(0,3,"allele_id",            affx::FILE5_DTYPE_INT);
                m_feffects_tsv.defineColumn(0,4,"context_id",           affx::FILE5_DTYPE_INT);
                m_feffects_tsv.defineColumn(0,5,"channel_id",           affx::FILE5_DTYPE_INT);
            } 
            else {
                m_feffects_tsv.defineStringColumn(0,0,"probeset_id",m_MaxNameLength);
                //m_feffects_tsv.defineColumn(0,1,"atom_id"          ,affx::FILE5_DTYPE_INT);
                m_feffects_tsv.defineColumn(0,1,"probe_id",        affx::FILE5_DTYPE_INT);
                //m_feffects_tsv.defineColumn(0,3,"x"                ,affx::FILE5_DTYPE_INT);
                //m_feffects_tsv.defineColumn(0,4,"y"                ,affx::FILE5_DTYPE_INT);
                m_feffects_tsv.defineColumn(0,2,"feature_response",affx::FILE5_DTYPE_DOUBLE);
                m_feffects_tsv.defineColumn(0,3,"allele_id",       affx::FILE5_DTYPE_INT);
                m_feffects_tsv.defineColumn(0,4,"context_id",      affx::FILE5_DTYPE_INT);
                m_feffects_tsv.defineColumn(0,5,"channel_id",      affx::FILE5_DTYPE_INT);
            }
        }
        m_feffects_tsv.writeTsv_v1();

    }



    // residuals
    if (m_DoResiduals) {
        m_residuals_tsv.defineStringColumn(0,0,"probeset_id",m_MaxNameLength);
    
        m_residuals_tsv.defineColumn(0,1,"atom_id"    ,affx::FILE5_DTYPE_INT);
        m_residuals_tsv.defineColumn(0,2,"probe_id"   ,affx::FILE5_DTYPE_INT);
        m_residuals_tsv.defineColumn(0,3,"x"          ,affx::FILE5_DTYPE_INT);
        m_residuals_tsv.defineColumn(0,4,"y"          ,affx::FILE5_DTYPE_INT);
    
        for (size_t i=0;i<col_names.size();i++) {
            m_residuals_tsv.defineColumn(0,5+i,col_names[i],affx::FILE5_DTYPE_DOUBLE);
        }
    
        m_residuals_tsv.writeTsv_v1();
    }

  
    return true;
}


bool QuantMethodExprReport::report(ProbeSetGroup &psGroup,
                                   QuantMethod &qMethod,
                                   const IntensityMart &iMart, 
                                   std::vector<ChipStream *> &iTrans,
                                   PmAdjuster &pmAdjust) {
 
 
    QuantExprMethod* qeMethod = dynamic_cast<QuantExprMethod *>(&qMethod);
    if (qeMethod == NULL) {
        Err::errAbort("Can't call QuantMethodExprReport::report() with something other than a QuantExprMethod.");
    }
    int targetCount = qeMethod->getNumTargets();

    /* We only report expression and copynumber probe sets. */
    if (psGroup.probeSets[0]->psType != ProbeSet::Expression && psGroup.probeSets[0]->psType != ProbeSet::Copynumber) {
        return false;
    }

    if (m_ProbeSetsToReport != NULL)
        if (m_ProbeSetsToReport->find(psGroup.name)==m_ProbeSetsToReport->end())
            return false;

    if (m_DoSummary) {
        // summary
		if (psGroup.displayName != NULL && strlen(psGroup.displayName) > 0)
			m_summary_tsv.set_string(0,0,psGroup.displayName);
		else
			m_summary_tsv.set_string(0,0,psGroup.name);
        for (int i = 0; i < targetCount; i++) {
            if (m_DoCompact) {
                m_summary_tsv.set_f(0,1+i,(float)qeMethod->getSignalEstimate(i));
            }
            else {
                m_summary_tsv.set_d(0,1+i,qeMethod->getSignalEstimate(i));
            }
        }
        m_summary_tsv.writeLevel(0);
    }

    // populate the map of id to names.

    std::map<int,std::string> probeToPsName;
    std::map<int,int> probeToAtom;
    std::map<int,int> probeToAllele;
    std::map<int,int> probeToContext;
    std::map<int,int> probeToChannel;

    std::map<int,std::string> probeToPsNameFE;
    std::map<int,int> probeToAtomFE;
    std::map<int,int> probeToAlleleFE;
    std::map<int,int> probeToContextFE;
    std::map<int,int> probeToChannelFE;
  
    if (m_DoResiduals || m_DoFeatureEffects) {
        for (uint32_t gIx = 0; gIx < psGroup.probeSets.size(); gIx++) {
            const ProbeSet *ps = psGroup.probeSets[gIx];
            for (uint32_t aIx = 0; aIx < ps->atoms.size(); aIx++) {
                Atom *a = ps->atoms[aIx];
                for (uint32_t pIx = 0; pIx < a->probes.size(); pIx++) {
                    Probe *p = a->probes[pIx];

                    probeToPsName[p->id] = ps->name;
                    probeToAtom[p->id] = a->id;
                    probeToAllele[p->id] = a->getAlleleCode();
                    probeToContext[p->id] = a->getContextCode();
                    probeToChannel[p->id] = a->getChannelCode();

                    //  We need to have a separate set of maps for feature effects where the apid is used
                    //  since two probes with the same p->id could be in the same probeset.    
                    int apid = p->getApid();
                    probeToPsNameFE[apid] = ps->name;
                    probeToAtomFE[apid] = a->id;
                    probeToAlleleFE[apid] = a->getAlleleCode();
                    probeToContextFE[apid] = a->getContextCode();
                    probeToChannelFE[apid] = a->getChannelCode();
                }
            }
        }
    }

    if (m_DoFeatureEffects) {
        // @todo: put code here.
        for (int featureIx = 0; featureIx < qeMethod->getNumFeatures(); featureIx++) {
            if (qeMethod->featureUsed(featureIx)) {
                const Probe *p = qeMethod->getFeature(featureIx);
                assert(p);
                if (m_WriteOldStyleFeatureEffectsFile) {
                    unsigned int x,y;
                    ChipLayout::indexToXY(p->id, m_NumCol, x, y); 

                    if (m_feffects_tsv.m_format == affx::TsvReport::FMT_A5) {
                        m_feffects_tsv.set_i(0,0,FROM_ZEROBASED(p->id));
                        m_feffects_tsv.set_d(0,1,qeMethod->getFeatureEffect(featureIx));
                        m_feffects_tsv.set_string(0,2,probeToPsName[p->id]);
                        m_feffects_tsv.writeLevel(0);
                    } 
                    else {
                        m_feffects_tsv.set_string(0,0,probeToPsName[p->id]);
                        m_feffects_tsv.set_i(0,1,probeToAtom[p->id]);
                        m_feffects_tsv.set_i(0,2,FROM_ZEROBASED(p->id));
                        m_feffects_tsv.set_i(0,3,x);
                        m_feffects_tsv.set_i(0,4,y);
                        m_feffects_tsv.set_d(0,5,qeMethod->getFeatureEffect(featureIx));
                        m_feffects_tsv.writeLevel(0);
                    }                        
                } 
                else {
                    std::string probesetName = probeToPsNameFE[p->getApid()];
                    if (m_unmungePSName && 
                        (probesetName.rfind("-A") == probesetName.length()-2 || 
                         probesetName.rfind("-B") == probesetName.length()-2)
                        ) {
                        probesetName.erase(probesetName.length()-2);
                    }
                    if (m_feffects_tsv.m_format == affx::TsvReport::FMT_A5) {     
                        m_feffects_tsv.set_i(0,0,FROM_ZEROBASED(p->id));
                        m_feffects_tsv.set_d(0,1,qeMethod->getFeatureEffect(featureIx));
                        m_feffects_tsv.set_string(0,2,probesetName);
                        m_feffects_tsv.set_i(0,3,probeToAlleleFE[p->getApid()]);
                        m_feffects_tsv.set_i(0,4,probeToContextFE[p->getApid()]);
                        m_feffects_tsv.set_i(0,5,probeToChannelFE[p->getApid()]);
                        m_feffects_tsv.writeLevel(0);
                    } 
                    else {
                        m_feffects_tsv.set_string(0,0,probesetName);
                        m_feffects_tsv.set_i(0,1,FROM_ZEROBASED(p->id));
                        m_feffects_tsv.set_d(0,2,qeMethod->getFeatureEffect(featureIx));
                        m_feffects_tsv.set_i(0,3,probeToAlleleFE[p->getApid()]);
                        m_feffects_tsv.set_i(0,4,probeToContextFE[p->getApid()]);
                        m_feffects_tsv.set_i(0,5,probeToChannelFE[p->getApid()]);
                        m_feffects_tsv.writeLevel(0);
                    }
                }
            }
        }
    }

    if (m_DoResiduals) {
        for (int featureIx = 0; featureIx < qeMethod->getNumFeatures(); featureIx++) {
            if (qeMethod->featureUsed(featureIx)) {
                const Probe *p = qeMethod->getFeature(featureIx);
                assert(p);
                unsigned int x,y;
                ChipLayout::indexToXY(p->id, m_NumCol, x, y);
        
                m_residuals_tsv.set_string(0,0,probeToPsName[p->id]);
        
                m_residuals_tsv.set_i(0,1,probeToAtom[p->id]);
                m_residuals_tsv.set_i(0,2,p->id+1);
                m_residuals_tsv.set_i(0,3,x);
                m_residuals_tsv.set_i(0,4,y);
        
                for (int chipIx = 0; chipIx < qeMethod->getNumTargets(); chipIx++) {
                    m_residuals_tsv.set_d(0,5+chipIx,qeMethod->getResidual(featureIx, chipIx));
                }
        
                m_residuals_tsv.writeLevel(0);
            }
        }
    }
  
    return true;
}

bool QuantMethodExprReport::finish(QuantMethod &qMethod) {
    //
    m_summary_tsv.close();
    if (m_DoFeatureEffects) {m_feffects_tsv.close();}
    if (m_DoResiduals) {m_residuals_tsv.close();}
    //
    close();
    //
    return true;
}
