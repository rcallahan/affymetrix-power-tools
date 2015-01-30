////////////////////////////////////////////////////////////////
//
// Copyright (C) 1989, 1991 Free Software Foundation, Inc.
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

#include "copynumber/CNAnalysisMethodCN.h" 
#include "copynumber/CNAnalysisMethodFactory.h"
#include "copynumber/CNAnalysisMethodGaussianSmooth.h"
#include "copynumber/CNAnalysisMethodLOH.h"
#include "copynumber/CNWorkflowEngine.h"
#include "copynumber/CPPTest/Setup.h"
//
#include "util/Util.h"
//
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/HelperMacros.h> 
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
//
#include <cstring>
#include <iostream>
#include <map>
#include <string>
//

using namespace std;
/**
 * @class CNAnalysisMethodFactoryTest
 * @brief cppunit class for testing CNAnalysisMethodFactory functions.
 * last change by vliber on 03/18/09
 */

class CNAnalysisMethodFactoryTest : public CppUnit::TestFixture   
{
  CPPUNIT_TEST_SUITE(CNAnalysisMethodFactoryTest);
  CPPUNIT_TEST(functionsTest);  
  CPPUNIT_TEST_SUITE_END();

public:  
  void functionsTest();
  
};
CPPUNIT_TEST_SUITE_REGISTRATION(CNAnalysisMethodFactoryTest);



void CNAnalysisMethodFactoryTest::functionsTest()
{

/*   Disabling this CNAnalysisMethodFactoryTest while code base is in flux for SNP7 analysis module.
	cout<<endl;
	Verbose::out(1, "****CNAnalysisMethodFactoryTest::functionsTest****");
	// constructor
	CNAnalysisMethodFactory cnAF;
	std::vector<SelfDoc> myDocs=cnAF.getDocs();
	CPPUNIT_ASSERT(myDocs.size()==25);
	CPPUNIT_ASSERT(myDocs[0].getState()=="pdnn-reference-method");
	CPPUNIT_ASSERT(myDocs[1].getState()=="wave-correction-reference-method.trim=2.0.percentile=0.75.wave-count=6.demean=false");
	CPPUNIT_ASSERT(myDocs[2].getState()=="additional-waves-reference-method.trim=2.0.percentile=0.75.additional-wave-count=1.demean=false.cn-qc-cutoff=0.27.snp-qc-cutoff=1.1.waviness-seg-count-cutoff=100.use-high-waviness-seg-count=true.force=false.keep-temp-data=false");
	CPPUNIT_ASSERT(myDocs[3].getState()=="pdnn-intensity-adjustment-method.predicted-intensity-bin-count=20.gc-bin-count=20.residual-trim=2.0");
	CPPUNIT_ASSERT(myDocs[4].getState()=="high-pass-filter-intensity-adjustment-method.data-block-rows=320.data-block-cols=2015.mini-block-rows=8.mini-block-cols=8.global-smooth-weight=256.0.local-smooth-weight=64.0.converged=0.0001");
	CPPUNIT_ASSERT(myDocs[5].getState()=="wave-correction-log2ratio-adjustment-method.bandwidth=101.bin-count=25.wave-count=-1.wave-smooth=true");
	CPPUNIT_ASSERT(myDocs[6].getState()=="cn-state.hmmCN_state=0,1,2,3,4.hmmCN_prior_prob=0.2,0.2,0.2,0.2,0.2.hmmCN_mu=-2,-0.533,0,0.363,0.567.hmmCN_sigma=0.2,0.2,0.2,0.2,0.2.hmmCN_TransitionDecay=1e+9.hmmCN_StateEstimationMethod=EM.hmmCN_EMIterations=1.hmmCN_EMConvergenceThreshold=0.0001.hmmCN_NormalState=2.hmmCN_ForwardOnly=0.hmmCN_NormalStateMinObservations=2.hmmCN_SmoothOutliers=1.hmmCN_TransTypeStat=0.PostCNFitMaxOutlierRemoveRunSize=1");
    CPPUNIT_ASSERT(myDocs[7].getState()=="cn-cyto2.hmmCN_state=0,1,2,3,4,5.hmmCN_mu=-1.63,-0.58,0,0.45,0.72,0.93.hmmCN_sigma=0.3,0.3,0.3,0.3,0.3,0.3.diagonal-weight=0.995.mapd-weight=0.22.min-segment-size=5.hmm-confidence-weight=0.6");
	CPPUNIT_ASSERT(myDocs[8].getState()=="cn-snp7.hmmCN_state=0,1,2,3,4,5.hmmCN_mu=-1.63,-0.58,0,0.45,0.72,0.93.hmmCN_sigma=0.3,0.3,0.3,0.3,0.3,0.3.diagonal-weight=0.995.min-segment-size=5.hmm-confidence-weight=0.6");
        
	CPPUNIT_ASSERT(myDocs[9].getState()=="log2-ratio.gc-correction=true.median-autosome-median-normalization=true.median-smooth-marker-count=5");
    CPPUNIT_ASSERT(myDocs[10].getState()=="log2-ratio-cyto2.gc-correction=true.median-autosome-median-normalization=true.median-smooth-marker-count=5.trim-high=2.0.trim-low=-2.5");
    CPPUNIT_ASSERT(myDocs[11].getState()=="allelic-difference.outlier-trim=3.0");
    CPPUNIT_ASSERT(myDocs[12].getState()=="gaussian-smooth.expSmoothSignal=1.smooth_sigma_multiplier=2.smooth_bw=50000");
    CPPUNIT_ASSERT(myDocs[13].getState()=="kernel-smooth.sigma_span=50.0");
    CPPUNIT_ASSERT(myDocs[14].getState()=="loh.lohCN_errorrate=0.05.lohCN_beta=0.001.lohCN_alpha=0.01.lohCN_separation=1000000.lohCN_nMinMarkers=10.lohCN_NoCallThreshold=0.05.lohCN_minGenomicSpan=1000000");
    CPPUNIT_ASSERT(myDocs[15].getState()=="loh-cyto2.lohCNSegSeparation=1000000.minInformation=100.lambdaCritical=8.0");
    CPPUNIT_ASSERT(myDocs[16].getState()=="cn-neutral-loh");
    CPPUNIT_ASSERT(myDocs[17].getState()=="normal-diploid");
    CPPUNIT_ASSERT(myDocs[18].getState()=="mosaicism.gains-boundries=0.08764945,0.15380349,0.21465931,0.27100300.losses-boundries=-0.08293345,-0.17551812,-0.28048196,-0.40165383.marker-bandwidth=6000.confidence-window=251");
    CPPUNIT_ASSERT(myDocs[19].getState()=="cn-gender.male-chrX-lower-threshold=0.8.male-chrX-upper-threshold=1.3.male-chrY-lower-threshold=0.8.male-chrY-upper-threshold=1.2.female-chrX-lower-threshold=1.9.female-chrX-upper-threshold=2.1.female-chrY-lower-threshold=0.female-chrY-upper-threshold=0.4.mapd-threshold=0.5");
    CPPUNIT_ASSERT(myDocs[20].getState()=="cn-cyto2-gender.cutoff=0.5");
    CPPUNIT_ASSERT(myDocs[21].getState()=="cn-segment");
    CPPUNIT_ASSERT(myDocs[22].getState()=="loh-segment");
    CPPUNIT_ASSERT(myDocs[23].getState()=="allele-peaks.step=50.window=250.point-count=128.bandwidth=0.45.cutoff=0.05.clean-threshold=0.25.symmetry=true");
    CPPUNIT_ASSERT(myDocs[24].getState()=="chipstream.NDBandwidth=10.0");
*/

}
