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
 * @file   ArtifactReduction.cpp
 * @author Earl Hubble
 * @date   Fri Oct 21 18:16:03 2005
 *
 * @brief Class for doing normalization.
 */

//
#include "chipstream/ArtifactReduction.h"
//
#include "chipstream/ProbeSet.h"
//
#include "algorithm/artifact/morphology.h"
#include "file/TsvFile/TsvFile.h"
#include "file5/File5.h"
#include "stats/stats.h"
#include "util/Err.h"
#include "util/Fs.h"
#include "util/FsPath.h"
#include "util/Verbose.h"

using namespace std;
using namespace affx;

/**
 * Constructor.
 *
 */
ArtifactReduction::ArtifactReduction() {
    setUpSelfDoc(*this);
    m_Type = getDocName();
    m_SafetyZero=0.1;

    //
    m_TransformedIMart = NULL;
    m_DataCache = NULL;
    m_NumChannels = 0;
    m_a5_filename="";
    m_a5_shared_group=NULL;
    m_pobjChipLayout = NULL;
    m_FreeLayout = false;
    m_CurCount = 0;
    m_Clip = -1;
    m_Open = -1;
    m_Close = -1;
    m_Fringe = -1;
    m_ResType = -1;
    m_MapVerbose = -1;
    m_CoincidenceCount = -1;
    m_TrustCheck = false;
    m_dc_k = 0;
    m_Gradient = 0;

}

/**
 * Destructor.
 */
ArtifactReduction::~ArtifactReduction() {

    for (unsigned int i = 0; i < m_ToFree.size(); i++) {
        delete m_ToFree[i];
    }
    if (m_FreeLayout) {
        Freez(m_pobjChipLayout);
    }
    Freez(m_DataCache);
    m_ProbesetTrustTsv.clear();
}

/**
 * Fill in the information for Self documentation.
 * @param doc - Self documenter to be filled in.
 */
void ArtifactReduction::setUpSelfDoc(SelfDoc &doc) {
    doc.setDocName(ARTIFACTREDUCTIONSTR);
    doc.setDocDescription("Class for artifact reduction.");
    doc.setDocOptions(getDefaultDocOptions());
}

/**
 * @brief Default Getter method for parameters and their documentation.
 * @return map of parameters and their descriptions.
 */
std::vector<SelfDoc::Opt> ArtifactReduction::getDefaultDocOptions() {
    std::vector<SelfDoc::Opt> opts;

    SelfDoc::Opt o_Clip = {"Clip", SelfDoc::Opt::Double, "2.0","0.0","0","NA",
                           "Threshold for intensity log-ratio to be outlier"};
    opts.push_back(o_Clip);

    SelfDoc::Opt o_Open = {"Open",SelfDoc::Opt::Integer, "1","0","0","NA",
                           "Minimum size of snow (isolated residuals) to remove"};
    opts.push_back(o_Open);

    SelfDoc::Opt o_Close = {"Close",SelfDoc::Opt::Integer, "3","0","0","NA",
                            "Max gap to close between tentative artifacts"};
    opts.push_back(o_Close);

    SelfDoc::Opt o_Fringe= {"Fringe", SelfDoc::Opt::Integer, "2","0","0","NA",
                            "Add a fringe to catch border features that might be poor."};
    opts.push_back(o_Fringe);

    SelfDoc::Opt o_ResType = {"ResType",SelfDoc::Opt::Integer, "2","1","1","NA",
                              "Type=1 raw intensity, type=2 replicate comparison"};
    opts.push_back(o_ResType);

    SelfDoc::Opt o_MapVerbose = {"MapVerbose", SelfDoc::Opt::Integer,"0","0","0","NA",
                                 "0=no map, 1=map, 2=map+profiles"};
    opts.push_back(o_MapVerbose);
    SelfDoc::Opt o_CC = {"CC", SelfDoc::Opt::Integer, "1","1","1","NA",
                         "How many channels must support an artifact"};
    opts.push_back(o_CC);
    SelfDoc::Opt o_tc = {"TrustCheck",SelfDoc::Opt::Boolean,"false","false","NA","NA",
                         "Completely untrusted, blimeshed probesets in samples are no-calls"};
    opts.push_back(o_tc);
    SelfDoc::Opt o_dc = {"DC",SelfDoc::Opt::Integer,"8","0","NA","NA",
                         "Region of block to average for gradient"};
    opts.push_back(o_dc);
    SelfDoc::Opt o_Gradient = {"Gradient",SelfDoc::Opt::Integer,"0","0","NA","NA",
                               "Do gradient removal"};
    opts.push_back(o_Gradient);


    return opts;
}


/**
 * @brief This static function should be overridden by child classes
 * to return an object of the correct type initialized correctly
 * with the parameters in the string, string map. All objects
 * created this way should be deleted when finished using.
 *
 * @param param - Map of key/value pairs to initialize the object.
 *
 * @return Pointer toCreate object, this should be sub casted as necessary.
 */
SelfCreate *ArtifactReduction::newObject(std::map<std::string,std::string> &param) {
    SelfDoc doc = explainSelf();

    ArtifactReduction *qnorm = new ArtifactReduction();

    // created generic, now fill in values from param string
    fillInValue(qnorm->m_Clip, "Clip", param,doc);
    fillInValue(qnorm->m_Open, "Open", param, doc);
    fillInValue(qnorm->m_Close,"Close",param, doc);
    fillInValue(qnorm->m_ResType,"ResType",param,doc);
    fillInValue(qnorm->m_MapVerbose,"MapVerbose",param,doc);
    fillInValue(qnorm->m_Fringe,"Fringe",param,doc);
    fillInValue(qnorm->m_CoincidenceCount, "CC", param, doc);
    fillInValue(qnorm->m_TrustCheck,"TrustCheck",param,doc);
    fillInValue(qnorm->m_dc_k,"DC",param,doc);
    fillInValue(qnorm->m_Gradient,"Gradient",param,doc);

    qnorm->setOptValue("Clip", ToStr(qnorm->m_Clip));
    qnorm->setOptValue("Open", ToStr(qnorm->m_Open));
    qnorm->setOptValue("Close", ToStr(qnorm->m_Close));
    qnorm->setOptValue("ResType", ToStr(qnorm->m_ResType));
    qnorm->setOptValue("MapVerbose", ToStr(qnorm->m_MapVerbose));
    qnorm->setOptValue("Fringe", ToStr(qnorm->m_Fringe));
    qnorm->setOptValue("CC", ToStr(qnorm->m_CoincidenceCount));
    qnorm->setOptValue("TrustCheck", qnorm->m_TrustCheck);
    qnorm->setOptValue("DC", ToStr(qnorm->m_dc_k));
    qnorm->setOptValue("Gradient", ToStr(qnorm->m_Gradient));

    return qnorm;

    


}

/// Now we interface with the numerics below here


void ArtifactReduction::updateReference(std::vector<float> &data) {
    if (m_CurCount==0) {
        m_ReferenceProfile=data;
        for (int i=0; i<data.size(); i++)
            m_ReferenceProfile[i] = log(max(data[i],m_SafetyZero));
        m_CurCount = 1;
    }
    else{
        // compute stable on-line mean
        // So we don't have to read in all the data at any time
        // Or access things in weird orders
        m_CurCount++;
        for (int i=0; i<data.size(); i++) {
            m_ReferenceProfile[i] +=
                (log(max(data[i],m_SafetyZero)) -
                 m_ReferenceProfile[i]) / m_CurCount;
        }
    }
}

void ArtifactReduction::updateAllReference(std::vector<std::vector<float> > &adata) {
    if (m_CurCount==0) {
        m_AllReferenceProfile.clear();
        m_AllReferenceProfile.resize(adata.size());
        for (int cd=0; cd<adata.size(); cd++) {
            m_AllReferenceProfile[cd].assign(adata[cd].size(),0);
            for (int i=0; i<adata[cd].size(); i++)
                m_AllReferenceProfile[cd][i] = log(max(adata[cd][i],
                                                       m_SafetyZero));
        }
        m_CurCount = 1;
    }
    else{
        // compute stable on-line mean
        // So we don't have to read in all the data at any time
        // Or access things in weird orders
        m_CurCount++;
        for (int cd=0; cd<adata.size(); cd++) {
            for (int i=0; i<adata[cd].size(); i++) {
                m_AllReferenceProfile[cd][i] +=
                    (log(max(adata[cd][i],m_SafetyZero)) -
                     m_AllReferenceProfile[cd][i]) / m_CurCount;
            }
        }
    }
}


float ArtifactReduction::TypeIResidual(std::vector<float> &resid,
                                       std::vector<float> &data,
                                       std::vector<float> &ref) {
    // construct residuals to an intensity profile

    float mt;
    resid=data;
    mt = 0;
    for (int i=0; i<data.size(); i++) {
        resid[i] = log(max(data[i],m_SafetyZero))-ref[i];
        mt += (resid[i]-mt)/(i+1);
    }
    //  Verbose::out(3, "Mean " + ToStr(mt) + " is the mean of me");
    for (int i=0; i<data.size(); i++) {
        resid[i] -= mt;
    }

    return(mt);
}

void ArtifactReduction::TypeTwoResidual(std::vector<float> &resid,
                                        std::vector<float> &data) {
    // set up residual matrix to be zero everywhere
    resid.clear();
    resid.assign(data.size(),0);

    int goodcount;
    float goodmean;
    long index;
    long psetCount;
    std::vector<int> probeIds;
    // iterate over all probe sets
    //
    psetCount = m_pobjChipLayout->getProbeSetCount();
    for (unsigned int psIx = 0; psIx < psetCount; psIx++) {
        ProbeListPacked pl = m_pobjChipLayout->getProbeListAtIndex(psIx);
        // Perhaps we only want to consider genotyping/marker probeset types
        if ((pl.get_type() == ProbeSet::GenoType) ||
            (pl.get_type() == ProbeSet::Marker) ||
            (pl.get_type() == ProbeSet::MultichannelMarker)) {
            int blockCount=pl.block_cnt();
            pl.get_probeIds(probeIds);
            int pIxStart = 0;
            for (int bIx=0;bIx<blockCount;bIx++) {
                // We now assume that all of the following probes in this block
                // are identical or perhaps we only make this assumption for
                // genotype/marker probeset types
                int blockSize=pl.get_blockSize(bIx);
                goodcount = 0;
                goodmean = 0;
                for (int pIx=pIxStart; pIx < pIxStart + blockSize; pIx++) {
                    if (1 /*pl.isPm(pIx)*/) {
                        // Do something
                        //index = pl.get_probeId(pIx);
                        index = probeIds[pIx];
                        goodcount++;
                        goodmean += (log(max(data[index],m_SafetyZero)) -
                                     goodmean) / goodcount;
                    }
                } // End for loop over probes

                // now get residuals for all probes from mean
                // median might be better, but should be good enough
                for (int pIx=pIxStart; pIx < pIxStart + blockSize; pIx++) {
                    if (1 /*pl.isPm(pIx)*/) {
                        // Do something
                        //index = pl.get_probeId(pIx);
                        index = probeIds[pIx];
                        resid[index] = log(max(data[index],m_SafetyZero))-goodmean;
                    }
                } // End for loop over probes again
                pIxStart += blockSize; // update by current block to go to next block

            } // End for loop over blocks
        } // End if genotype/marker probeset type
    } // End for loop over probesets
}

void ArtifactReduction::EmptyBlemishes(std::vector<int> &tmpBlemish,
                                       std::vector<float> &tmpRes) {
    // make no blemishes here
    tmpBlemish.assign(tmpRes.size(),0);
}

long ArtifactReduction::ThresholdResiduals(std::vector<int> &tmpBlemish,
                                           std::vector<float> &tmpRes,
                                           float Clip) {
    long dumb_count=0;

    //tmpBlemish.assign(tmpRes.size(),0);

    for (int i=0; i<tmpRes.size(); i++) {
        if (tmpRes[i]>Clip) {
            tmpBlemish[i]+=1;
            dumb_count++;
        }
        if (tmpRes[i]< (-1*Clip)) {
            tmpBlemish[i]+=1;
            dumb_count++;
        }
    }
    return(dumb_count);
}

int ArtifactReduction::EraseBlemishToReference(std::vector<float> &data,
                                               std::vector<float> &ref,
                                               std::vector<int> &blemishes,
                                               float mt) {
    // erases blemishes by resetting to reference
    // this is the worst possible method and is only used for testing
    int dumb_count=0;
    for (int i=0; i<blemishes.size(); i++) {
        if (blemishes[i]>0) {
            data[i] = exp(ref[i] + mt);
            dumb_count++;
        }
    }
    return(dumb_count);
}

void ArtifactReduction::FixProbeIterationList() {
    // generate a probe iteration list over replicates
    //int replace_count=0;
    int goodcount;
    float goodmean;
    long index=0;
    long psetCount;
    // iterate over all probe sets
    //
    m_ProbeIterationList.clear();
    vector<int> tmpReplicateList;

    psetCount = m_pobjChipLayout->getProbeSetCount();
    for (unsigned int psIx = 0; psIx < psetCount; psIx++) {
        ProbeListPacked pl = m_pobjChipLayout->getProbeListAtIndex(psIx);
        // Perhaps we only want to consider genotyping/marker probeset types
        if ((pl.get_type() == ProbeSet::GenoType) ||
            (pl.get_type() == ProbeSet::Marker) ||
            (pl.get_type() == ProbeSet::MultichannelMarker)) {

            int blockCount=pl.block_cnt();
            int pIxStart = 0;
            for (int bIx=0;bIx<blockCount;bIx++) {
                // We now assume that all of the following probes in this block
                // are identical or perhaps we only make this assumption for
                // genotype/marker probeset types
                int blockSize=pl.get_blockSize(bIx);
                goodcount = 0;
                goodmean = 0;
                tmpReplicateList.clear();
                for (int pIx=pIxStart; pIx < pIxStart + blockSize; pIx++) {
                    if (1 ) {//pl.isPm(pIx)
                        // Do something
                        index = pl.get_probeId(pIx);
                        tmpReplicateList.push_back(index);
                    }
                } // End for loop over probes
                m_ProbeIterationList.push_back(tmpReplicateList);
                pIxStart += blockSize; // update by current block to go to next block

            } // End for loop over blocks
        } // End if genotype/marker probeset type
    } // End for loop over probesets
}


void ArtifactReduction::NewTypeTwoResidual(std::vector<float> &resid,
                                           std::vector<float> &data) {
    // set up residual matrix to be zero everywhere
    resid.clear();
    resid.assign(data.size(),0);

    int goodcount;
    float goodmean;
    long index;
    long psetCount;
    // iterate over all probe sets
    //
    psetCount = m_ProbeIterationList.size();
    for (unsigned int psIx = 0; psIx < psetCount; psIx++) {
        goodmean=0;
        goodcount = 0;

        for (int pIx=0; pIx < m_ProbeIterationList[psIx].size(); pIx++) {
            if (1 /*pl.isPm(pIx)*/) {
                // Do something
                index = m_ProbeIterationList[psIx][pIx];
                goodcount++;
                goodmean +=
                    (log(max(data[index],m_SafetyZero))-goodmean)/goodcount;
            }
        } // End for loop over probes

        // now get residuals for all probes from mean
        // median might be better, but should be good enough
        for (int pIx=0; pIx < m_ProbeIterationList[psIx].size(); pIx++) {
            if (1 /*pl.isPm(pIx)*/) {
                // Do something
                index = m_ProbeIterationList[psIx][pIx];
                resid[index] = log(max(data[index],m_SafetyZero))-goodmean;
            }
        } // End for loop over probes again

    } // End for loop over probesets
}

int ArtifactReduction::NewEraseBlemishToUnblemished(std::vector<float> &data,
                                                    std::vector<int> &blemishes) {

    int replace_count=0;
    int goodcount;
    float goodmean;
    long index=0;
    long psetCount;
    // iterate over all probe sets
    //
    psetCount = m_ProbeIterationList.size();
    for (unsigned int psIx = 0; psIx < psetCount; psIx++) {


        goodcount = 0;
        goodmean = 0;
        for (int pIx=0; pIx < m_ProbeIterationList[psIx].size(); pIx++) {
            if (1 ) {//pl.isPm(pIx)
                // Do something
                index = m_ProbeIterationList[psIx][pIx];
                if (blemishes[index]<1) {
                    goodcount++;
                    goodmean += (data[index]-goodmean)/goodcount;
                }
            }
        } // End for loop over probes

        // now substitute unblemished for blemished if possible
        for (int pIx=0; pIx < m_ProbeIterationList[psIx].size(); pIx++) {
            if (1 ) {//pl.isPm(pIx)
                // Do something
                index = m_ProbeIterationList[psIx][pIx];
                if (blemishes[index]>0 && goodcount>0) {
                    data[index] = goodmean;
                    replace_count++;
                }
            }
        } // End for loop over probes again


    } // End for loop over probesets
    return(replace_count);
}

int ArtifactReduction::EraseBlemishToUnblemished(std::vector<float> &data,
                                                 std::vector<int> &blemishes) {

    int replace_count=0;
    int goodcount;
    float goodmean;
    long index=0;
    long psetCount;
    std::vector<int> probeIds;

    // iterate over all probe sets
    //
    psetCount = m_pobjChipLayout->getProbeSetCount();
    for (unsigned int psIx = 0; psIx < psetCount; psIx++) {
        ProbeListPacked pl = m_pobjChipLayout->getProbeListAtIndex(psIx);
        // Perhaps we only want to consider genotyping/marker probeset types
        if ((pl.get_type() == ProbeSet::GenoType) ||
            (pl.get_type() == ProbeSet::Marker) ||
            (pl.get_type() == ProbeSet::MultichannelMarker)) {

            int blockCount=pl.block_cnt();
            pl.get_probeIds(probeIds);
            int pIxStart = 0;
            for (int bIx=0;bIx<blockCount;bIx++) {
                // We now assume that all of the following probes in this block are identical
                // or perhaps we only make this assumption for genotype/marker probeset types
                int blockSize=pl.get_blockSize(bIx);
                goodcount = 0;
                goodmean = 0;
                for (int pIx=pIxStart; pIx < pIxStart + blockSize; pIx++) {
                    if (1 ) {//pl.isPm(pIx)
                        // Do something
                        //index = pl.get_probeId(pIx);
                        index = probeIds[pIx];
                        if (blemishes[index]<1) {
                            goodcount++;
                            goodmean += (data[index]-goodmean)/goodcount;
                        }
                    }
                } // End for loop over probes

                // now substitute unblemished for blemished if possible
                for (int pIx=pIxStart; pIx < pIxStart + blockSize; pIx++) {
                    if (1 ) {//pl.isPm(pIx)
                        // Do something
                        //index = pl.get_probeId(pIx);
                        index = probeIds[pIx];
                        if (blemishes[index]>0 && goodcount>0) {
                            data[index] = goodmean;
                            replace_count++;
                        }
                    }
                } // End for loop over probes again
                pIxStart += blockSize; // update by current block to go to next block

            } // End for loop over blocks
        } // End if genotype/marker probeset type
    } // End for loop over probesets
    return(replace_count);
}

int ArtifactReduction::BlemishedByProbeset(int ChipId, std::vector<int> &blemishes) {

    int replace_count=0;
    int goodcount;
    float goodmean;
    long index=0;
    long psetCount;
    std::vector<int> probeIds;
    
    //m_ProbesetBlemishStream << ChipId;
    // iterate over all probe sets
    //
    psetCount = m_pobjChipLayout->getProbeSetCount();
    for (unsigned int psIx = 0; psIx < psetCount; psIx++) {
        ProbeListPacked pl = m_pobjChipLayout->getProbeListAtIndex(psIx);
        // Perhaps we only want to consider genotyping/marker probeset types
        if ((pl.get_type() == ProbeSet::GenoType) || (pl.get_type() == ProbeSet::Marker) || (pl.get_type() == ProbeSet::MultichannelMarker)) {
            int blockCount=pl.block_cnt();
            pl.get_probeIds(probeIds);
            int pIxStart = 0;
            int bestcount=probeIds.size();  // at most number of probes trusted probes
            for (int bIx=0;bIx<blockCount;bIx++) {
                // We now assume that all of the following probes in this block are identical
                // or perhaps we only make this assumption for genotype/marker probeset types
                int blockSize=pl.get_blockSize(bIx);
                goodcount = 0;
                goodmean = 0;
                for (int pIx=pIxStart; pIx < pIxStart + blockSize; pIx++) {
                    if (1 ) {//pl.isPm(pIx)
                        // Do something
                        //index = pl.get_probeId(pIx);
                        index = probeIds[pIx];
                        if (blemishes[index]<1){
                            goodcount++;
                        }
                    }
                } // End for loop over probes
                if (goodcount<bestcount)
                    bestcount=goodcount;
                
                
                pIxStart += blockSize; // update by current block to go to next block
                
            } // End for loop over blocks
            // now we output 'worst' of all blocks
            // since we're just trying to make trust/notrust decision
            //  only output things not to trust!
            //  If we don't find them, they're trusted!
            if (bestcount<1){
                m_ProbesetTrustTsv.set(0,"chp_id",ToStr(ChipId));
                m_ProbesetTrustTsv.set(0,"probeset_id",pl.get_name_string());
                m_ProbesetTrustTsv.set(0,"trust_count",bestcount);
                m_ProbesetTrustTsv.writeLevel(0);
            }
        } // End if genotype/marker probeset type
    } // End for loop over probesets
    
    //m_ProbesetBlemishStream << endl;
    return(replace_count);
}

int sum_map(vector<int> &myb) {
    int tot=0;
    for (int i=0; i<myb.size(); i++)
        tot += myb[i];
    return(tot);
}


void threshold_map(vector<int> &myb, int t) {
    for (int i=0; i<myb.size(); i++)
        if (myb[i]>=t)
            myb[i] =1;
        else
            myb[i] = 0;
}

void ArtifactReduction::FixBlemish(std::vector<std::vector<float> > &adata,
                                   float Clip, int chipNo) {
    vector<float> tmp;
    vector<int> myblemishes;

    float mt;
    long dumb_count,x_count, raw_count, proc_count;

    dumb_count =0;

    // dumbest possible blemish corrector
    // Clip large residuals to end
    EmptyBlemishes(myblemishes,adata[0]);

    // blemishes aggregated over all channels
    // Assume channels represent same physical object here
    for (int cd=0; cd<adata.size(); cd++) {

        if (m_ResType==1)
            mt=TypeIResidual(tmp,adata[cd], m_AllReferenceProfile[cd]);
        if (m_ResType==2)
            TypeTwoResidual(tmp,adata[cd]);

        // save actual residuals here
        if (m_MapVerbose>0)
            saveResidualsToFile("residuals."+ToStr(chipNo)+".channel" +
                                ToStr(cd) + ".txt",tmp,"ResType" +
                                ToStr(m_ResType));

        dumb_count+=ThresholdResiduals(myblemishes,tmp,Clip);
    }

    // focus on ones for which sufficient channels agree there is a blemish
    threshold_map(myblemishes,m_CoincidenceCount);

    raw_count= sum_map(myblemishes);

    int Xdim = m_pobjChipLayout->getXCount();
    int Ydim = m_pobjChipLayout->getYCount();

    //Verbose::out(1,ToStr(Xdim) + ":" + ToStr(Ydim));

    // "close" maneuver to create contiguous regions
    dilate(myblemishes, Xdim, Ydim, m_Close);
    erode(myblemishes, Xdim, Ydim,m_Close);

    // "open maneuver to eliminate 'snow'
    erode(myblemishes,Xdim,Ydim,m_Open);
    dilate(myblemishes,Xdim,Ydim,m_Open);

    // "Dilate" maneuver to add "fringe" to all contiguous regions
    dilate(myblemishes,Xdim,Ydim,m_Fringe);

    proc_count = sum_map(myblemishes);

    //x_count = EraseBlemishToReference(adata[0], m_AllReferenceProfile[0], myblemishes, mt);
    x_count=0;
    for (int cd=0; cd<adata.size(); cd++)
        x_count += EraseBlemishToUnblemished(adata[cd],myblemishes);

    // output by construction
    if (m_TrustCheck) {
        BlemishedByProbeset(chipNo,myblemishes);
    }

    // blemish map constructed
    // output for debugging purposes
    if (m_MapVerbose>0)
        saveBlemishMapToFile("artifact."+ ToStr(chipNo) + ".map.txt",
                             myblemishes, "artifact");


    Verbose::out(3, "Found " + ToStr(raw_count*100.0/(Xdim*Ydim)) + " raw residuals " + ToStr(proc_count*100.0/(Xdim*Ydim)) + " processed");


}

void fix_residuals(vector<float> &data, vector<float> &res){
    for (int i=0; i<data.size(); i++){
        data[i] *= exp(-res[i]);
    }
}

void ArtifactReduction::FixGradient(std::vector<std::vector<float> > &adata, float Clip, int chipNo){
    // fixes gradients on the image, after worst artifacts have been removed
    // fixes gradient on each channel independently, as they can have different magnitude
    vector<float> tmp;
    
    float mt;
    
    int Xdim = m_pobjChipLayout->getXCount();
    int Ydim = m_pobjChipLayout->getYCount();
    
    // Assume channels represent same physical object here
    for (int cd=0; cd<adata.size(); cd++){
        
        // recompute residuals, as artifact reduction has changed things all over the snps
        if (m_ResType==1)
            mt=TypeIResidual(tmp,adata[cd], m_AllReferenceProfile[cd]);
        if (m_ResType==2)
            TypeTwoResidual(tmp,adata[cd]);
        
        // save actual residuals here
        if (m_MapVerbose>0)
            saveResidualsToFile("Gresiduals."+ToStr(chipNo)+".channel" +ToStr(cd) + ".txt",tmp,"ResType"+ToStr(m_ResType));
        
        // actually deal with the residuals in this channel
        // if m_dc_k were 0, get back residual matrix identity
        dc_block(tmp,Xdim,Ydim,m_dc_k);
        // remove residual
        fix_residuals(adata[cd],tmp);
    }
}

void ArtifactReduction::Winsorize(std::vector<float> &data,
                                  std::vector<float> &ref,
                                  float Clip) {
    vector<float> tmp;
    vector<int> myblemishes;

    float mt;
    long dumb_count;

    // dumbest possible blemish corrector
    // Clip large residuals to end
    dumb_count=0;
    mt=TypeIResidual(tmp,data,ref);

    for (int i=0; i<data.size(); i++) {
        if (tmp[i]>Clip) {
            data[i] = exp(ref[i]+Clip+mt);
            dumb_count++;
        }
        if (tmp[i]< (-1*Clip)) {
            data[i] = exp(ref[i]-Clip+mt);
            dumb_count++;
        }
    }
    Verbose::out(3, "Found " + ToStr(dumb_count) + " residuals");
}

void ArtifactReduction::GetReferenceProfile(IntensityMart* iMart){

    int channel_count = iMart->getChannelCount();
    //std::vector<string> CelFiles = iMart->getCelFileNames();
    m_CelFiles = iMart->getCelFileNames();
    int file_count = m_CelFiles.size();

    std::vector<float> data;
    std::vector<vector <float> > adata; // data over all channels

    if (m_ReferenceProfileName.empty()) {
        Verbose::progressBegin(1, "Computing reference profile for " +
                               ToStr(file_count) + " cel files",
                               file_count, 0, file_count);

        m_ReferenceProfile.clear();
        m_CurCount = 0;
        for (int d = 0; d < file_count; d++) {
            Verbose::progressStep(1);
            adata.clear();
            for (int cd=0; cd<channel_count; cd++){
                data = iMart->getCelData(d, cd);
                adata.push_back(data);
            }
            //updateReference(data);
            updateAllReference(adata);
        }

        Verbose::progressEnd(1, "Done.");
    }
    else{
        Verbose::out(1, "Reading reference profile from: " +
                     m_ReferenceProfileName);
        readReferenceProfileFromFile(m_ReferenceProfileName);
    }
}

float ArtifactReduction::transform(int probeIx, int chipIx, float intensity, PsBoard &board) {
    int dsIx = chipIx / m_NumChannels;
    int channelIx = chipIx % m_NumChannels;
    APT_ERR_ASSERT(m_TransformedIMart == NULL, "Not expecting m_TransformedIMart to be set. Cel file: " +
                   m_TransformedIMart->getCelFileNames()[chipIx]);
    return m_DataCache->getProbeIntensity(probeIx, dsIx, channelIx);
}

float ArtifactReduction::transform(int probeIx, int chipIx, float intensity) {
    int dsIx = chipIx / m_NumChannels;
    int channelIx = chipIx % m_NumChannels;
    APT_ERR_ASSERT(m_TransformedIMart == NULL, "Not expecting m_TransformedIMart to be set. Cel file: " +
                   m_TransformedIMart->getCelFileNames()[chipIx]);
    return m_DataCache->getProbeIntensity(probeIx, dsIx, channelIx);
}

void ArtifactReduction::doReduction(IntensityMart* iMart) {
    
    // if there aren't any chipstream nodes after this, then don't store
    // all of the intensities
    if (m_Streams.empty()) {
        iMart->setStoreAllCelIntensities(false);
    }

    int channel_count = iMart->getChannelCount();

    //std::vector<string> CelFiles = iMart->getCelFileNames();
    int file_count = iMart->getCelFileCount();

    if (m_TrustCheck){
        m_ProbesetTrustTsv.defineColumn(0, 0, "chp_id");
        m_ProbesetTrustTsv.defineColumn(0, 1, "probeset_id");
        m_ProbesetTrustTsv.defineColumn(0, 2, "trust_count");
        APT_ERR_ASSERT( ! m_ProbesetTrustTmpFileName.empty(), "Empty temporay trusted file name.");
        m_ProbesetTrustTsv.writeTsv_v1( m_ProbesetTrustTmpFileName );
    }

    std::vector<float> data;

    std::vector<vector <float> > adata; // data over all channels

    if (m_ResType==1 || !m_ReferenceProfileName.empty() ||
        !m_FileOutName.empty())

        GetReferenceProfile(iMart);

    Verbose::progressBegin(1, "Applying artifact reduction to " +
                           ToStr(file_count) + " cel files", file_count, 0,
                           file_count);
    //Verbose::out(1, "Caching replicate probe list to apply to all cels");
    //FixProbeIterationList();

    for (int d = 0; d < file_count; d++) {
        Verbose::progressStep(1);
        adata.clear(); // remove entries from adata
        for (int cd=0; cd<channel_count; cd++) {
            data = iMart->getCelData(d, cd);
            adata.push_back(data); // add an element to adata
        }

        //Winsorize(data,2);
        if (m_MapVerbose==2)
            saveMultiChannelProfileToFile("origin." + ToStr(d) + ".txt",
                                          adata,"clip-profile");

        FixBlemish(adata,m_Clip,d);
        if (m_Gradient)
            FixGradient(adata,m_Clip,d);

        if (m_MapVerbose==2)
            saveMultiChannelProfileToFile("profile." + ToStr(d) + ".txt",
                                          adata,"clip-profile");

        for (int cd=0; cd<channel_count; cd++){
            iMart->setProbeIntensity(channel_count * d + cd,
                                     adata[cd]);}
    }
    Verbose::progressEnd(1, "Done.");

    if (m_TrustCheck) {
        m_ProbesetTrustTsv.clear();
    }

}

/**
 * @brief Method for being passed a new cel file worth of data.
 * @param data - Vector of vectors of cel file data.
 */
void ArtifactReduction::newDataSet(IntensityMart* iMart) {
    // set up the output Disk Mark
    m_TransformedIMart = iMart;
    doReduction(m_TransformedIMart);
    chipStreamPassNewChip(m_TransformedIMart);
}

void ArtifactReduction::newChip(std::vector<float> &data) {
    // If the data cache isn't set up yet, then get it going
    if (m_DataCache == NULL) {
        m_DataCache = new DataStore(m_DataStoreTempFile);
        vector<int> desiredOrder(data.size(), -1);
        for (int i = 0; i < data.size(); i++) {
            desiredOrder[i]  = i;
        }
        m_DataCache->setProbeOrder(desiredOrder);
        m_DataCache->initChannels(m_NumChannels);
    }
    // Just store the data for now, we'll process it in finishedChips
    int dataCount = Max(m_DataCache->getCelDataSetCount(), 0);
    m_DataCache->writeColumn(dataCount, ToStr(dataCount), data);
    if (dataCount % m_NumChannels == 0) {
        m_CacheCelNames.push_back(ToStr(m_CacheCelNames.size()));
    }
}

void ArtifactReduction::finishedChips() {
    m_DataCache->setCelFileNames(m_CacheCelNames);
    doReduction(m_DataCache);
}

/**
 * @brief Signal that no more data is coming (i.e. newChip() will not be
 * called anymore. Currently the class builds up all the sketches and then
 * does the normalization when this function is called. This makes it
 * difficult to pass through data to downstream. If a precomputed sketch
 * is supplied then the data can be passed through chip by chip.
 */
void ArtifactReduction::endDataSet() {
    /* Flush data in stream. */
    chipStreamEndChips();


    if (m_FileOutName !="") {
        saveMultiChannelProfileToFile(m_FileOutName, m_AllReferenceProfile,
                                      "mc-profile");
    }
    if (m_a5_filename!="") {
        saveReferenceProfileToFile_a5(m_a5_filename);
    }
    if (m_a5_shared_group!=NULL) {
        saveReferenceProfileToGroup_a5(m_a5_shared_group);
    }
}

/// i/o stuff past here

/**
 * Open and read a target distribution from a text file. Must have a column
 * called 'intensities' with distribution quantiles sorted from highest to
 * lowest.
 *
 * @param fileName - text file containing column of quantiles.
 */
void ArtifactReduction::readReferenceProfileFromFile(const std::string& fileName) {
    TsvFile tsv;
    double intensity;
    // file ok
    if (File5_File::isHdf5file(fileName)==1) {
        Err::errAbort("A5 implementation of reference-profile is not yet implemented.");
    }
    int rc = tsv.open(fileName);
    assert(rc==TSV_OK);
    //
    tsv.bind(0,"intensities",&intensity,affx::TSV_BIND_REQUIRED);

    m_ReferenceProfile.clear();
    while (tsv.nextLevel(0)==TSV_OK) {
        m_ReferenceProfile.push_back(intensity);
    }

    tsv.close();
}


void ArtifactReduction::readMultiChannelReferenceProfileFromFile(const std::string& fileName) {
    TsvFile tsv;
    double intensity = -1;
    // file ok
    if (File5_File::isHdf5file(fileName)==1) {
        Err::errAbort("A5 implementation of reference-profile is not yet implemented.");
    }
    int rc = tsv.open(fileName);
    assert(rc==TSV_OK);
    //
    //tsv.bind(0,"intensities",&intensity,affx::TSV_BIND_REQUIRED);
    vector<int> ch_cidx;
    ch_cidx.assign(tsv.getColumnCount(0),-1); // no data yet for channels
    for (int cd=0; cd<ch_cidx.size(); cd++) {
        ch_cidx[cd]=tsv.cname2cidx(0,"channel"+ToStr(cd));
    }

    m_AllReferenceProfile.clear();
    m_AllReferenceProfile.resize(ch_cidx.size());
    while (tsv.nextLevel(0)==TSV_OK) {
        for (int cd=0; cd<ch_cidx.size(); cd++) {
            if (ch_cidx[cd]>=0) {
                tsv.get(0,ch_cidx[cd], intensity);
            }
            m_AllReferenceProfile[cd].push_back(intensity);
	}
    }

    tsv.close();
}

void ArtifactReduction::readReferenceProfileFromFile_a5(const std::string& fileName) {
    affx::File5_File* file5=new affx::File5_File();
    file5->open(fileName,affx::FILE5_OPEN_RO);

    readReferenceProfileFromFile_a5_tsv(file5,"reference-profile");

    file5->close();
    delete file5;
}

void ArtifactReduction::readReferenceProfileFromFile_a5_tsv(affx::File5_Group* group5, const std::string& name) {
    //printf("### reading A5 profile from: '%s'\n",name.c_str());
    double intensity;
    affx::File5_Tsv* tsv5=group5->openTsv("reference-profile",
                                          affx::FILE5_OPEN_RO);

    m_ReferenceProfile.clear();
    while (tsv5->nextLine()==affx::FILE5_OK) {
        tsv5->get(0,0,&intensity);
        m_ReferenceProfile.push_back(intensity);
    }
    //
    tsv5->close();
    delete tsv5;
}

void ArtifactReduction::saveReferenceProfileToFile(const std::string &fileName) {
    saveProfileToFile(fileName,m_ReferenceProfile,"reference-profile");
}

/**
 * Save the current m_ReferenceProfile to a file.
 * @param fileName - where to write the profile.
 */
void ArtifactReduction::saveProfileToFile(const std::string& fileName,
                                          std::vector<float> &Profile,
                                          const std::string &ProfileType) {

    affx::TsvFile tsv;
    tsv.addHeader("tsv-file-type",ProfileType);
    tsv.defineFile("intensities");
    tsv.setPrecision(0,0,8);
    tsv.writeTsv_v1(fileName);
    //
    for (unsigned int i = 0; i < Profile.size(); i++) {
        tsv.set(0,0,Profile[i]);
        tsv.writeLevel(0);
    }
    tsv.close();
}

/**
 * Save the current residuals to a file.
 * @param fileName - where to write the profile.
 */
void ArtifactReduction::saveResidualsToFile(const std::string& fileName,
                                            std::vector<float> &Residuals,
                                            const std::string &ResidualType) {

    affx::TsvFile tsv;
    tsv.addHeader("tsv-file-type",ResidualType);
    tsv.defineFile("LogResiduals");
    tsv.setPrecision(0,0,8);
    tsv.writeTsv_v1(fileName);
    //
    for (unsigned int i = 0; i < Residuals.size(); i++) {
        tsv.set(0,0,Residuals[i]);
        tsv.writeLevel(0);
    }
    tsv.close();
}

/**
 * Save the current m_ReferenceProfile to a file.
 * @param fileName - where to write the profile.
 */
void ArtifactReduction::saveMultiChannelProfileToFile(const std::string& fileName, std::vector<std::vector<float> > &Profile, const std::string &ProfileType) {

    affx::TsvFile tsv;
    string t_col;

    tsv.addHeader("tsv-file-type",ProfileType);

    for (int cd=0; cd<Profile.size(); cd++) {
        t_col = "channel" + ToStr(cd);
        tsv.defineColumn(0,cd,t_col);
        tsv.setPrecision(0,cd, 6);
    }
    tsv.writeTsv_v1(fileName);
    // loop over all entries, loop over all columns
    for (unsigned int i = 0; i < Profile[0].size(); i++) {
        for (int cd=0; cd<Profile.size(); cd++) {
            tsv.set(0,cd,Profile[cd][i]);
        }
        tsv.writeLevel(0);
    }
    tsv.close();
}


/**
 * Save the current BlemishMap to a file.
 * @param fileName - where to write the profile.
 */
void ArtifactReduction::saveBlemishMapToFile(const std::string& fileName,
                                             std::vector<int> &Profile,
                                             const std::string &ProfileType) {

    affx::TsvFile tsv;
    tsv.addHeader("tsv-file-type",ProfileType);
    tsv.defineFile("artifact");
    tsv.setPrecision(0,0,8);
    tsv.writeTsv_v1(fileName);
    //
    for (unsigned int i = 0; i < Profile.size(); i++) {
        tsv.set(0,0,Profile[i]);
        tsv.writeLevel(0);
    }
    tsv.close();
}

void ArtifactReduction::saveReferenceProfileToFile_a5(const std::string& fileName) {
    try {
        affx::File5_File* file5=new affx::File5_File();
        file5->open(fileName,affx::FILE5_REPLACE);

        saveReferenceProfileToGroup_a5(file5);

        file5->close();
        delete file5;
    }
    catch(...) {
        Err::errAbort("Cannot open file: " + fileName);
    }
}

void ArtifactReduction::saveReferenceProfileToGroup_a5(affx::File5_Group* group5) {
    //
    affx::File5_Tsv* tsv5=group5->openTsv("reference-profile",affx::FILE5_CREATE);
    tsv5->defineColumn(0,0,"intensities",affx::FILE5_DTYPE_DOUBLE,0);
    //
    for (unsigned int i = 0; i < m_ReferenceProfile.size(); i++) {
        tsv5->set_d(0,0,m_ReferenceProfile[i]);
        tsv5->writeLevel(0);
    }
    //
    tsv5->close();
    delete tsv5;
}

std::string ArtifactReduction::createProbesetTrustTmpFileName( Options *opt)
{
    std::string trustFileName = Fs::basename(GlobalTmpFileFactory()->genFilename_basic(std::string(ARTIFACTREDUCTIONSTR) + ".",".trust"));
    trustFileName = Fs::join(opt->getOpt("temp-dir"), trustFileName + ".tsv" );
    return trustFileName;
}

void ArtifactReduction::setParameters(PsBoard &board) {
    string spfFile;
    // @todo refactor - Just read the chiplayout directly from blackboard for now...ChipLayout
    board.get("chiplayout", &m_pobjChipLayout);
    m_NumChannels = m_pobjChipLayout->numChannels();
    Options *o = board.getOptions();
    setreadReferenceProfile(o->getOpt("reference-profile"));
    // @todo refactor - get the write profile working by setting the output directory
    // if (!o->getOptBool("write-profile")) {
    //     string outFile = m_WriteProfileDir + description + s->getType() + ToStr(".generic-target.txt");
    //     saveReferenceProfile(outFile);
    // }
    
    m_ProbesetTrustTmpFileName = createProbesetTrustTmpFileName(o);
    
    board.set(getProbesetTrustFileKey(), m_ProbesetTrustTmpFileName);
    m_DataStoreTempFile = Fs::join(o->getOpt("temp-dir"),"artifact-reduction-temp.a5");
}
