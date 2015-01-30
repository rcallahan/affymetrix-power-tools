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
#include "chipstream/AnalysisInfo.h"
#include "chipstream/TsvReport.h"
#include "dmet/DmetCopyNumberData.h"
#include "dmet/DmetCopyNumberEngine.h"
#include "file5/File5.h"
#include "file5/File5_Tsv.h"
#include "stats/stats-distributions.h"
#include "util/Err.h"
#include "util/Fs.h"
#include "util/Guid.h"
#include "util/Util.h"
//
#include <cmath>
#include <iostream>
#include <vector>
//

using namespace affx;
using namespace std;

DmetCopyNumberEngine::Reg DmetCopyNumberEngine::reg;

DmetCopyNumberEngine * DmetCopyNumberEngine::FromBase(BaseEngine *engine)
{
	if (engine != NULL && engine->getEngineName() == DmetCopyNumberEngine::EngineName())
		return (DmetCopyNumberEngine *)engine;
	return NULL;
}

double DmetCopyNumberEngine::minDouble(double a, double b) { return (a<=b) ? a : b; }

double DmetCopyNumberEngine::maxDouble(double a, double b) { return (a>b) ? a : b; }

double DmetCopyNumberEngine::Log2(double input) { return ( log(input) / log(2.0) ) ; }

void DmetCopyNumberEngine::parseRegionDataFile(string regionDataFileName, DmetCopyNumberData &copyNumberData ){

  affx::TsvFile tsv;
  if(tsv.open(regionDataFileName) != affx::TSV_OK)
      Err::errAbort("Unable to open region file '"+regionDataFileName+"'");
  while(tsv.nextLevel(0)==TSV_OK){
    string inputString;
    double inputDouble;

    if(tsv.get(0,"region",inputString)==TSV_OK)
      copyNumberData.addRegion(inputString);
    else
      Err::errAbort("Unable to parse region file: missing column: region.");

    if(tsv.get(0, "nocall_prob", inputDouble)==TSV_OK)
      copyNumberData.setNoCallProbability(inputString, inputDouble); 
    else
      Err::errAbort("Unable to parse region file: missing column: nocall_prob.");

    if(tsv.get(0, "outlier_prob", inputDouble)==TSV_OK)
      copyNumberData.setOutlierProbability(inputString, inputDouble); 
    else
      Err::errAbort("Unable to parse region file: missing column: outlier_prob.");

    if(tsv.get(0, "prior_0", inputDouble)==TSV_OK)
      copyNumberData.setPrior_0(inputString, inputDouble); 
    else
      Err::errAbort("Unable to parse region file: missing column: prior_0.");

    if(tsv.get(0, "prior_1", inputDouble)==TSV_OK)
      copyNumberData.setPrior_1(inputString, inputDouble); 
    else
      Err::errAbort("Unable to parse region file: missing column: prior_1.");

    if(tsv.get(0, "prior_2", inputDouble)==TSV_OK)
      copyNumberData.setPrior_2(inputString, inputDouble); 
    else
      Err::errAbort("Unable to parse region file: missing column: prior_2.");

    if(tsv.get(0, "s_01_0", inputDouble)==TSV_OK)
      copyNumberData.setS_01_0(inputString, inputDouble); 
    else
      Err::errAbort("Unable to parse region file: missing column: s_01_0.");

    if(tsv.get(0, "m_01_0", inputDouble)==TSV_OK)
      copyNumberData.setM_01_0(inputString, inputDouble); 
    else
      Err::errAbort("Unable to parse region file: missing column: m_01_0.");

    if(tsv.get(0, "s_01_1", inputDouble)==TSV_OK)
      copyNumberData.setS_01_1(inputString, inputDouble); 
    else
      Err::errAbort("Unable to parse region file: missing column: s_01_1.");

    if(tsv.get(0, "m_01_1", inputDouble)==TSV_OK)
      copyNumberData.setM_01_1(inputString, inputDouble); 
    else
      Err::errAbort("Unable to parse region file: missing column: m_01_1.");

    if(tsv.get(0, "s_12_1", inputDouble)==TSV_OK)
      copyNumberData.setS_12_1(inputString, inputDouble); 
    else
      Err::errAbort("Unable to parse region file: missing column: s_12_1.");

    if(tsv.get(0, "m_12_1", inputDouble)==TSV_OK)
      copyNumberData.setM_12_1(inputString, inputDouble); 
    else
      Err::errAbort("Unable to parse region file: missing column: m_12_1.");

    if(tsv.get(0, "s_12_2", inputDouble)==TSV_OK)
      copyNumberData.setS_12_2(inputString, inputDouble); 
    else
      Err::errAbort("Unable to parse region file: missing column: s_12_2.");

    if(tsv.get(0, "m_12_2", inputDouble)==TSV_OK)
      copyNumberData.setM_12_2(inputString, inputDouble); 
    else
      Err::errAbort("Unable to parse region file: missing column: m_12_2.");

    if(tsv.get(0, "n_0", inputDouble)==TSV_OK)
      copyNumberData.setN_0(inputString, inputDouble); 
    else
      Err::errAbort("Unable to parse region file: missing column: n_0.");

    if(tsv.get(0, "n_1", inputDouble)==TSV_OK)
      copyNumberData.setN_1(inputString, inputDouble); 
    else
      Err::errAbort("Unable to parse region file: missing column: n_1.");

    if(tsv.get(0, "n_2", inputDouble)==TSV_OK)
      copyNumberData.setN_2(inputString, inputDouble); 
    else
      Err::errAbort("Unable to parse region file: missing column: n_2.");
    
  }
}

void DmetCopyNumberEngine::parseProbeSetDataFile(string probeSetDataFileName, DmetCopyNumberData &copyNumberData ){

  affx::TsvFile tsv;
  if(tsv.open(probeSetDataFileName) != affx::TSV_OK)
      Err::errAbort("Unable to open region file '"+probeSetDataFileName+"'");
  while(tsv.nextLevel(0)==TSV_OK){
    string inputString1, inputString2;
    double inputDouble;

    if(tsv.get(0,"probeset_id",inputString1)==TSV_OK){

      if(tsv.get(0, "region", inputString2)==TSV_OK)
        copyNumberData.addPredictionProbeSet(inputString2, inputString1); 
      else
        Err::errAbort("Unable to parse probeset model file: missing column: retion.");

      if(tsv.get(0, "m_0", inputDouble)==TSV_OK)
        copyNumberData.setMean_0(inputString1, inputDouble); 
      else
        Err::errAbort("Unable to parse probeset model file: missing column: m_0.");

      if(tsv.get(0, "m_1", inputDouble)==TSV_OK)
        copyNumberData.setMean_1(inputString1, inputDouble); 
      else
        Err::errAbort("Unable to parse probeset model file: missing column: m_1.");

      if(tsv.get(0, "m_2", inputDouble)==TSV_OK)
        copyNumberData.setMean_2(inputString1, inputDouble); 
      else
        Err::errAbort("Unable to parse probeset model file: missing column: m_2.");

      if(tsv.get(0, "w01", inputDouble)==TSV_OK)
         copyNumberData.setW_01(inputString1, inputDouble); 
       else
        Err::errAbort("Unable to parse probeset model file: missing column: w01.");

      if(tsv.get(0, "w12", inputDouble)==TSV_OK)
        copyNumberData.setW_12(inputString1, inputDouble); 
      else
        Err::errAbort("Unable to parse probeset model file: missing column: w12.");
    } else // This else is for the failure to get the first inputString2. 
        Err::errAbort("Unable to parse probeset model file: missing column: probeset_id.");
  }// end loop for nextLevel.
  
    
}


void DmetCopyNumberEngine::parseCNRegionGenotypeProbesetFile(string cn_region_gt_marker_filename, DmetCopyNumberData &copyNumberData) {
  affx::TsvFile cn_region_gt_markers;
  if ( cn_region_gt_markers.open(cn_region_gt_marker_filename) != TSV_OK ) {
    Err::errAbort(cn_region_gt_marker_filename + ": failed to open CN region genotype probesets file.");
  }
  else {
    while (cn_region_gt_markers.nextLevel(0) == TSV_OK) {
      string probeset;
      string region;
      if(cn_region_gt_markers.get(0, "probeset_id", probeset) != TSV_OK){
        Err::errAbort(cn_region_gt_marker_filename + ": cannot read probeset field.");
      }
      if(cn_region_gt_markers.get(0, "region", region) != TSV_OK){
        Err::errAbort(cn_region_gt_marker_filename + ": cannot read CN region field.");
      }

      copyNumberData.addGenotypeProbeSet(region, probeset);
    }
  }
}

//  Calculates the copyNumber value by interpolating probeSet intensity values between known CN0 and CN1 intensity
//  values.  Then take a weighted sum over all the probeSets. 
 
double DmetCopyNumberEngine::calculate01CopyNumberEstimate(string region, 
            set<string> probeSetNames, string sampleName,
            const DmetCopyNumberData &copyNumberData){
 
  double slope=0.0; 
  double cumulativeSampleSum=0.0;

  cumulativeSampleSum=0.0;
 
  set<string>::iterator beginProbeSet=probeSetNames.begin();
  set<string>::iterator endProbeSet=probeSetNames.end();

  for( ; beginProbeSet!=endProbeSet; beginProbeSet++){
    string probeSetName = *beginProbeSet;

    slope=maxDouble(0.000001, copyNumberData.getMean_1(probeSetName) - copyNumberData.getMean_0(probeSetName));            
    cumulativeSampleSum+= ( ( (copyNumberData.getSampleProbeSetIntensity(sampleName, probeSetName) - 
				copyNumberData.getMean_0(probeSetName)) / slope) * 
				copyNumberData.getW_01(probeSetName));

/* Debugging
  double mean0 = copyNumberData.getMean_0(probeSetName);
  double mean1 = copyNumberData.getMean_1(probeSetName);
  double intensity = copyNumberData.getSampleProbeSetIntensity(sampleName, probeSetName);
  double w_01 =  copyNumberData.getW_01(probeSetName);
*/
	 
   } 
    
  return cumulativeSampleSum;

}


double DmetCopyNumberEngine::calculate12CopyNumberEstimate(string region, 
    set<string> probeSetNames, string sampleName,
    const DmetCopyNumberData &copyNumberData){
 
  double slope=0.0; 
  double cumulativeSampleSum=0.0;

  cumulativeSampleSum=0.0;
 
  set<string>::iterator beginProbeSet=probeSetNames.begin();
  set<string>::iterator endProbeSet=probeSetNames.end();
  for( ; beginProbeSet!=endProbeSet; beginProbeSet++){

    slope=maxDouble(0.000001, copyNumberData.getMean_2(*beginProbeSet) - copyNumberData.getMean_1(*beginProbeSet));            
    cumulativeSampleSum+= ( ( 1.0 + ( (copyNumberData.getSampleProbeSetIntensity(sampleName, *beginProbeSet) - 
        copyNumberData.getMean_1(*beginProbeSet)) / slope) ) * 
        copyNumberData.getW_12(*beginProbeSet));
    } 
    
  return cumulativeSampleSum;

}


double DmetCopyNumberEngine::calculateCNDensity( double inputValue,
        double mean, double stDeviation, double degFreedom){

  if(degFreedom < 3)
    Err::errAbort("DmetCopyNumberengine::calculateCNDensity: deg freedom < 3");
  degFreedom = maxDouble(3.0, (degFreedom-1.0));
  double aa = pow((degFreedom/(degFreedom-2.0)), 0.5)/stDeviation;
  double ww = aa*(inputValue - mean);
  double output = (aa * affxstat::tDistribution(ww, degFreedom));
  return output;
 
}
  
int DmetCopyNumberEngine::findMaxIndex(const vector<double> inputVector){
  if(inputVector.size() < 1)
      Err::errAbort("Empty input vector in DmetCopyNumberEngine::findMaxIndex.");
  double presentMax = -1e-25;
  int maxIndex = 0;
  int index = 0;
  std::vector<double>::const_iterator begin = inputVector.begin(); 
  std::vector<double>::const_iterator  end = inputVector.end();
  for( ; begin!=end; begin++){
    if(inputVector[index] > presentMax){ 
      presentMax = inputVector[index];
      maxIndex = index;
    }
    index++; 
  }
  return maxIndex;
}

int DmetCopyNumberEngine::tailQuantileLower(const vector<double> inputVector, double pValue){
  if(inputVector.size() < 1)
      Err::errAbort("Empty input vector in DmetCopyNumberEngine::findMaxIndex.");
  if(pValue < 0 || pValue > 1)
      Err::errAbort("Invalid pvalue '"+ToStr(pValue)+"'");
  
  double cumulativeSum=0.0;
  for(size_t i=0; i<inputVector.size(); i++) {
    cumulativeSum += inputVector[i];
    if(cumulativeSum >= pValue)
      return i;
  }
  //  If the function makes it to here the sum of all the elements in the inputVector
  //  is still less than pValue and we return the index of the last element in the vector.
  return inputVector.size() - 1;
} 
    

int DmetCopyNumberEngine::tailQuantileUpper(const vector<double> inputVector, double pValue){
  if(inputVector.size() < 1)
      Err::errAbort("Empty input vector in DmetCopyNumberEngine::findMaxIndex.");
  if(pValue < 0 || pValue > 1)
      Err::errAbort("Invalid pvalue '"+ToStr(pValue)+"'");


  int returnIndex = inputVector.size()-1; 
  double cumulativeSum=0.0;

  vector<double>::const_reverse_iterator begin= inputVector.rbegin();
  vector<double>::const_reverse_iterator end= inputVector.rend();
  for( ; begin!=end; begin++){
    cumulativeSum += inputVector[returnIndex];
    if(cumulativeSum >= pValue)
      return returnIndex;
    returnIndex--;
  }
  //  If the function makes it to here the sum of all the elements in the inputVector
  //  is still less than pVal and we return the index of the first element in the vector.
  //  ie. 0.
  return 0; 
} 
 
void DmetCopyNumberEngine::parseSummaryDataFile5(string summaryDataFileName, DmetCopyNumberData &copyNumberData){
    string str;
    vector<string> sampleNames;
    affx::File5_File file5;
    file5.open(summaryDataFileName, affx::FILE5_OPEN_RO);
    affx::File5_Tsv* tsv5 = NULL;
    tsv5 = file5.openTsv("dmet.summary", affx::FILE5_OPEN_RO);
    if(tsv5 == NULL)
        Err::errAbort("Unable to get summary data from a5 file");

    // get the experiment names
    unsigned int colCount = tsv5->getColumnCount(0);
    for (unsigned int colIndex = 1; (colIndex < colCount); colIndex++)
    {
        tsv5->getColumnName(0, colIndex, &str);
        sampleNames.push_back(str);
        copyNumberData.addSampleName(str);
        //Verbose::out(1,"Found sample '"+str+"'in summary file");
    }

    // Load up signal estimates
    double estimate = 0;
    string probeSetName;
    while (tsv5->nextLine() == affx::FILE5_OK) {
            tsv5->get(0, 0, &probeSetName);
            //Verbose::out(1,"Found probeset '"+probeSetName+"'in summary file");
            for (unsigned int colIndex = 1; colIndex < colCount; colIndex++) {
                tsv5->get(0, colIndex, &estimate);
                copyNumberData.setSampleProbeSetIntensity(sampleNames[colIndex-1], probeSetName, Log2(estimate));
                //Verbose::out(1,"Setting sample data: " + sampleNames[colIndex-1] + ", " + probeSetName + ", " + ToStr(estimate));
            }
    }

    // Close/free resources
    tsv5->close();
    delete tsv5;
    file5.close();

}

void DmetCopyNumberEngine::parseSummaryDataFile(string summaryDataFileName, DmetCopyNumberData &copyNumberData){

  affx::TsvFile tsv;
  if(tsv.open(summaryDataFileName) != affx::TSV_OK)
      Err::errAbort("Unable to open region file '"+summaryDataFileName+"'");

  //  First we add the sample names to out copyNumberData data structure.
  int numberOfColumns = tsv.getColumnCount(0);
  vector<string> sampleNames;
  for(int i=1; i<numberOfColumns; i++){
    string columnName = tsv.getColumnName(0,i);
    copyNumberData.addSampleName(columnName);
    // We save the sample names for later use.
    sampleNames.push_back(columnName); 
  }    

  // We now loop over the probeset collecting and storing the intensity data for the 
  // probeset/sample values. 
  int columnIndex=1; 
  double estimate;
  while(tsv.nextLevel(0)==TSV_OK){
    string probeSetName;
    if(tsv.get(0,0,probeSetName)){

      // Note that we do not add ProbeSet names to our DmetCopyNumberData structure at this point.
      // That is done on reading the ProbeSet specific parameter file.  We do a check when summary
      // data is accessed to verify that intensity data for the probes in the ProbeSet specific 
      // does exist in the structure that is filled at this point. 

      columnIndex=1;
      while(columnIndex < numberOfColumns) {
        tsv.get(0, columnIndex, estimate);
        copyNumberData.setSampleProbeSetIntensity(sampleNames[columnIndex-1], probeSetName, Log2(estimate));
        columnIndex++;
      } // end of while loop over samples
    } // end of if
    else  
        Err::errAbort("Unable to read in summary file.");
  }
}

void DmetCopyNumberEngine::configureReport(TsvReport *tsv, affx::TsvReport::TsvReportFmt_t format, AnalysisInfo &info, const string &guid) {

  tsv->setDirPath(getOpt("out-dir"));
  tsv->setFileprefix(getOpt("set-analysis-name"));
  tsv->setFilename(tsv->getFileprefix()+".copynumber");

  tsv->setFormat(format);

  tsv->defineStringColumn(0,0,"cel",         TSVREPORT_CELFILE_STRLEN);
  tsv->defineStringColumn(0,1,"region",      TSVREPORT_PROBESET_STRLEN);
  tsv->defineColumn(      0,2,"cn_call",     affx::FILE5_DTYPE_INT);
  tsv->defineColumn(      0,3,"cn_force",    affx::FILE5_DTYPE_INT);
  tsv->defineColumn(      0,4,"cn_estimate", affx::FILE5_DTYPE_DOUBLE);
  tsv->defineColumn(      0,5,"cn_lower",    affx::FILE5_DTYPE_INT);
  tsv->defineColumn(      0,6,"cn_upper",    affx::FILE5_DTYPE_INT);
  tsv->defineColumn(      0,7,"pvalue",      affx::FILE5_DTYPE_DOUBLE);

  if (tsv->writeTsv_v1()!=affx::TSV_OK) {
    Err::errAbort("Cant open tsv report.");
  }

   tsv->addHeader("guid",guid);
   std::vector<std::string>::const_iterator keyIx, paramIx;
   for(keyIx = info.m_ParamNames.begin(), paramIx = info.m_ParamValues.begin();
       keyIx != info.m_ParamNames.end() && paramIx != info.m_ParamValues.end();
       ++keyIx, ++paramIx) {
          tsv->addHeader("affymetrix-algorithm-param-" + *keyIx, *paramIx);
   }
}

void DmetCopyNumberEngine::runImp(){

  vector<TsvReport *> reports;
  vector<affx::TsvReport*> gt_probesets_rpt;
  try {

    const string analysisGuid = affxutil::Guid::GenerateNewGuid();

    string regionDataFileName = getOpt("region-model");
    string probeSetDataFileName = getOpt("probeset-model");
    
    // Meta Info For Output File
    AnalysisInfo info;

    //  Output File
    if(getOptBool("text-output")) {
        TsvReport *report = new TsvReport;
        configureReport(report, affx::TsvReport::FMT_TSV, info, analysisGuid);
        reports.push_back(report);
    }
    
    if(getOptBool("a5-output")) {
        TsvReport *report = new TsvReport;
        configureReport(report, affx::TsvReport::FMT_A5, info, analysisGuid);
        reports.push_back(report);
    }
    
    /*
    affx::TsvFile tsv;
    string outputFileName = Fs::join(getOpt("out-dir"),"dmet-copynumber.txt");
    tsv.defineFile("cel\tregion\tcn_call\tcn_force\tcn_estimate\tcn_lower\tcn_upper\tpvalue");
    tsv.writeTsv_v1(outputFileName);
    */
    
    double copyNumberEstimate01=0.0;
    double copyNumberEstimate12=0.0;
    
    DmetCopyNumberData copyNumberData;
    set<string> probeSetNames;
    
   
    if(getOpt("summaries")!="")
        parseSummaryDataFile(getOpt("summaries"), copyNumberData);
    else if(getOpt("a5-summaries")!="")
        parseSummaryDataFile5(getOpt("a5-summaries"), copyNumberData);
    parseRegionDataFile(regionDataFileName, copyNumberData );
    parseProbeSetDataFile(probeSetDataFileName, copyNumberData );
    

    //  These are the samples and regions we will be examining. 
    vector<string> sampleNames = copyNumberData.getSampleNames();
    vector<string> regions = copyNumberData.getRegionNames();

    string cn_region_gt_marker_filename = getOpt("cn-region-gt-probeset-file");
    if (!cn_region_gt_marker_filename.empty()) {
        parseCNRegionGenotypeProbesetFile(cn_region_gt_marker_filename, copyNumberData);
    
	//  Output File Format
	if(getOptBool("a5-output")) {
	    affx::TsvReport *a5 = new TsvReport();
	    a5->setDirPath(getOpt("out-dir"));
	    a5->setFilename(getOpt("set-analysis-name") + ".gt.markers.cn.call");	    a5->setFormat(affx::TsvReport::FMT_A5);
	    gt_probesets_rpt.push_back(a5);
	}

	if(getOptBool("text-output")) {
	    affx::TsvReport *tsv = new TsvReport();
	    tsv->setDirPath(getOpt("out-dir"));
	    tsv->setFilename(getOpt("set-analysis-name") + ".gt.markers.cn.call");
	    tsv->setFormat(affx::TsvReport::FMT_TSV);
	    gt_probesets_rpt.push_back(tsv);
	}	

	for (int it = 0; it < gt_probesets_rpt.size(); it++) {
	    gt_probesets_rpt[it]->defineColumn(0, 0, "sample_id", affx::FILE5_DTYPE_STRING, 100);
	    
	    // make column headers using probeset ids
	    int cidx = 1;
	    for (int i = 0; i < regions.size(); i++) {
		vector<string> gt_probesets = copyNumberData.getGenotypeProbeSets(regions[i]);
		for (int j = 0; j < gt_probesets.size(); j++) {
		    gt_probesets_rpt[it]->defineColumn(0, cidx++, gt_probesets[j], affx::FILE5_DTYPE_INT);
		}
	    }
	    
	    // -1 is error code for TsvReport.  This will probably change,
	    // compelling a change here as well.
	    if (gt_probesets_rpt[it]->writeTsv_v1() != affx::TSV_OK) {
        Err::errAbort(gt_probesets_rpt[it]->getFilePath() + ": Unable to open for writing.");
	    }
	}
	
    }
    
    //  We iterate over all samples i.e. CELL files being examined.    
    vector<string>::iterator sampleNamesBegin= sampleNames.begin();
    vector<string>::iterator sampleNamesEnd= sampleNames.end();
    for( ; sampleNamesBegin!=sampleNamesEnd; sampleNamesBegin++){


        string sampleName = *sampleNamesBegin; 
        if (!cn_region_gt_marker_filename.empty()) {
	    for (int it = 0; it < gt_probesets_rpt.size(); it++) {
		gt_probesets_rpt[it]->set_string(0, 0, sampleName);
	    }
	}
	// index used for filling in the rest of the level for gt_probesets_rpt
	int gt_probesets_rpt_cidx = 1;

        // We now iterate over all CN regions.      
        VerboseErrHandler* errHandler=new VerboseErrHandler(true);
        vector<string>::iterator regionNamesBegin= regions.begin();
        vector<string>::iterator regionNamesEnd= regions.end();
        for( ; regionNamesBegin!=regionNamesEnd; regionNamesBegin++){

            string regionName = *regionNamesBegin; 
        
            bool error = false;
            string errorMsg = "";

            Verbose::setOutput(false);
            Err::pushHandler(errHandler);

            int copyNumberCall=-1;

            try {
    
                // Some variables used in calculating CN for each region.
                double sum=0.0;
                double sumWithOutliers=0.0;
                double copyNumberEstimate = 0.0; 
                vector<double> posteriors;
                vector<double> posteriorsWithOutliers;
    
                // These are the output values needed for each Copy Number regions. 
                int copyNumberForce=-1;
                double pValue = 0.0; 
                int copyNumberLower=-1;
                int copyNumberUpper=-1;
    
                // Set the priors from the parameter files
                double prior_0 = copyNumberData.getPrior_0(regionName);
                double prior_1 = copyNumberData.getPrior_1(regionName);
                double prior_2 = copyNumberData.getPrior_2(regionName);
    
                // Some local variables for calculation. 
                double d_01_0 = 0.0;
                double d_01_1 = 0.0;
                double d_12_1 = 0.0;
                double d_12_2 = 0.0;
                double copyNumberEstimate_b = 0.0;
    
    
                // We split the copy number determination into two sub processes.  In the first we assume that
                // there is CN=2 training data and do two subdeterminatnions and choose the "best".  In the second
                // case we assume that there is no CN=2 training data, so the choice is only between 0 and 1.
                // Note that the variables "sum" and "posteriors" are computed in both sub processes and used
                // to calculate other output.      
            
                probeSetNames=copyNumberData.getPredictionProbeSetNamesForRegion(regionName);  

                if( copyNumberData.getN_2(regionName) != 0 ){
                    copyNumberEstimate01= calculate01CopyNumberEstimate(regionName,
                                    probeSetNames, sampleName, copyNumberData);
    
                    copyNumberEstimate01 = minDouble(1.0, maxDouble(0.0, copyNumberEstimate01));
    
                    copyNumberEstimate12= calculate12CopyNumberEstimate(regionName,
                                    probeSetNames, sampleName, copyNumberData);
    
                    double copyNumberEstimate12_b = minDouble(2.5, copyNumberEstimate12);
    
                    d_01_0 =  calculateCNDensity(copyNumberEstimate01,
                        copyNumberData.getM_01_0(regionName),
                        copyNumberData.getS_01_0(regionName),
                        copyNumberData.getN_0(regionName));
    
                    d_01_1 =  calculateCNDensity(copyNumberEstimate01,
                        copyNumberData.getM_01_1(regionName),
                        copyNumberData.getS_01_1(regionName),
                        copyNumberData.getN_1(regionName));
    
                    d_12_1 =  calculateCNDensity(copyNumberEstimate12_b,
                        copyNumberData.getM_12_1(regionName),
                        copyNumberData.getS_12_1(regionName),
                        copyNumberData.getN_1(regionName));
    
                    d_12_2 =  calculateCNDensity(copyNumberEstimate12_b,
                        copyNumberData.getM_12_2(regionName),
                        copyNumberData.getS_12_2(regionName),
                        copyNumberData.getN_2(regionName));
    
    
                    double sum_01 = (prior_0*d_01_0) + ((1-prior_0)*d_01_1);
                    double post_gt0 = ((1-prior_0)*d_01_1)/sum_01;
    
                    sum = (prior_0*d_01_0) + (post_gt0*( (prior_1*d_12_1) + (prior_2*d_12_2) ) );
                    sumWithOutliers = sum + copyNumberData.getOutlierProbability(regionName);
    
                    posteriors.push_back(prior_0*d_01_0/sum);
                    posteriors.push_back((post_gt0*prior_1*d_12_1)/sum);
                    posteriors.push_back((post_gt0*prior_2*d_12_2)/sum);
            
                    posteriorsWithOutliers.push_back(prior_0*d_01_0/sumWithOutliers);
                    posteriorsWithOutliers.push_back((post_gt0*prior_1*d_12_1)/sumWithOutliers);
                    posteriorsWithOutliers.push_back((post_gt0*prior_2*d_12_2)/sumWithOutliers);
            
                    if( post_gt0 <= 0.5 || copyNumberEstimate12 < 1.0)
                        copyNumberEstimate = copyNumberEstimate01 ;
                    else
                        copyNumberEstimate = copyNumberEstimate12;
            
                } else {
        
                    copyNumberEstimate = calculate01CopyNumberEstimate(regionName,
                            probeSetNames, sampleName, copyNumberData);
                    copyNumberEstimate =  maxDouble(0.0, copyNumberEstimate);
                    copyNumberEstimate_b = minDouble(1.0, copyNumberEstimate);
    
                    d_01_0 =  calculateCNDensity(copyNumberEstimate_b,
                        copyNumberData.getM_01_0(regionName),
                        copyNumberData.getS_01_0(regionName),
                        copyNumberData.getN_0(regionName));
    
                    d_01_1 =  calculateCNDensity(copyNumberEstimate_b,
                        copyNumberData.getM_01_1(regionName),
                        copyNumberData.getS_01_1(regionName),
                        copyNumberData.getN_1(regionName));
    
                    sum=prior_0*d_01_0 + prior_1*d_01_1;
                    sumWithOutliers = sum + copyNumberData.getOutlierProbability(regionName);
    
                    posteriors.push_back(prior_0*d_01_0/sum);
                    posteriors.push_back(prior_1*d_01_1/sum);
    
                    posteriorsWithOutliers.push_back(prior_0*d_01_0/sumWithOutliers);
                    posteriorsWithOutliers.push_back(prior_1*d_01_1/sumWithOutliers);
    
                }// end else  

                copyNumberForce = findMaxIndex(posteriors);
                copyNumberCall = copyNumberForce;

                pValue = 1-posteriorsWithOutliers[copyNumberForce]; 

                ///@todo remove hack. Setting of calls should be determined using "label" field in region model file
                if(copyNumberCall >= 1) {
                    copyNumberCall = -3;
                }
    
                copyNumberLower = tailQuantileLower(posteriors, copyNumberData.getNoCallProbability(regionName));
                copyNumberUpper = tailQuantileUpper(posteriors, copyNumberData.getNoCallProbability(regionName)); 
    
                if(copyNumberData.getNoCallProbability(regionName) < pValue)
                    copyNumberCall = -1;
    
                posteriors.erase(posteriors.begin(), posteriors.end());
                posteriorsWithOutliers.erase(posteriorsWithOutliers.begin(), posteriorsWithOutliers.end());
    
                // Write out the results
                for(int i=0; i<reports.size(); i++) {
                    reports[i]->set_string(0, 0, sampleName);
                    reports[i]->set_string(0, 1, regionName);
                    reports[i]->set_i(     0, 2, copyNumberCall); 
                    reports[i]->set_i(     0, 3, copyNumberForce); 
                    reports[i]->set_d(     0, 4, copyNumberEstimate); 
                    reports[i]->set_i(     0, 5, copyNumberLower); 
                    reports[i]->set_i(     0, 6, copyNumberUpper); 
                    reports[i]->set_d(     0, 7, pValue); 
                    reports[i]->writeLevel(0);
                }

            }// end try-catch
            catch(Except &e) {
                error = true;
                errorMsg = e.what();
            }
            catch(...){
                error = true;
                errorMsg = "";
            }
            Err::popHandler();
            Verbose::setOutput(true);

            if(error) {
                copyNumberCall = -1;
                if(sampleNamesBegin == sampleNames.begin()) {
                    Verbose::out(1,"Unable to process region " + regionName + ". Most likely summary data not available/consented.");
                    if ( errorMsg.length() && (errorMsg.find("FATAL ERROR: ") != std::string::npos ) )
                      errorMsg = errorMsg.substr(errorMsg.find("FATAL ERROR: ")+13);
                    if(errorMsg != "")
                        Verbose::out(3,errorMsg);
                }
            }
            if (!cn_region_gt_marker_filename.empty()) {
                vector<string> gt_probesets = copyNumberData.getGenotypeProbeSets(regionName);
                for (int i = 0; i < gt_probesets.size(); i++) {
                    // Tsvreport does not (yet) let you access
                    // columns by name, so a column index is used.
		    for (int it = 0; it < gt_probesets_rpt.size(); it++) {
			gt_probesets_rpt[it]->set_i(0, gt_probesets_rpt_cidx, copyNumberCall);
		    }
		    gt_probesets_rpt_cidx++;
                }
            }

            // GUI Hack -- provide a chance to cancel
            Verbose::out(1,"",false);

            error = false;
        }// end loop over Regions 

        if (!cn_region_gt_marker_filename.empty()) {
	    for (int it = 0; it < gt_probesets_rpt.size(); it++) {
		gt_probesets_rpt[it]->writeLevel(0);
	    }
	}
    }// end loop over Samples.
  }
  catch(...) {
      for(int i=0; i<reports.size(); i++) {
          reports[i]->close();
          Freez(reports[i]);
      }
      reports.resize(0);
      for (int it = 0; it < gt_probesets_rpt.size(); it++) {	
	  gt_probesets_rpt[it]->close();
	  Freez(gt_probesets_rpt[it]);
      }
      gt_probesets_rpt.resize(0);
      throw;
  }
  for(int i=0; i<reports.size(); i++) {
      reports[i]->close();
      Freez(reports[i]);
  }
  for (int it = 0; it < gt_probesets_rpt.size(); it++) {	
      gt_probesets_rpt[it]->close();
      Freez(gt_probesets_rpt[it]);
  }
  reports.resize(0);
  gt_probesets_rpt.resize(0);
}// end mainFunction 
 
DmetCopyNumberEngine::DmetCopyNumberEngine() {
    defineOptions();
}

DmetCopyNumberEngine::~DmetCopyNumberEngine() {
}

void DmetCopyNumberEngine::defineOptions() {

  defineOptionSection("Input Options");
  defineOption("", "summaries", PgOpt::STRING_OPT,
                     "Probeset summary values. ",
                     "");
  defineOption("", "a5-summaries", PgOpt::STRING_OPT,
                     "Probeset summary values in a5 format. ",
                     "");
  defineOption("", "region-model", PgOpt::STRING_OPT,
                     "Regions model parameters. ",
                     "");
  defineOption("", "probeset-model", PgOpt::STRING_OPT,
                     "Probeset model parameters. ",
                     "");
  defineOption("", "cn-region-gt-probeset-file", PgOpt::STRING_OPT, 
                     "File of genotype probesets contained in copy number regions. ", "");

  defineOptionSection("Output Options");
  defineOption("", "text-output", PgOpt::BOOL_OPT,
                     "Generate text output file. ",
                     "false");
  defineOption("", "a5-output", PgOpt::BOOL_OPT,
                     "Generate a5 output file. ",
                     "true");

  defineOptionSection("Analysis Options");
  defineOption("", "set-analysis-name", PgOpt::STRING_OPT,
                     "Explicitly set the analysis name. "
                     "This affects output file names (ie prefix) and various meta info.",
                     "dmet");
}

void DmetCopyNumberEngine::defineStates() {
  defineOption("", "chip-type", PgOpt::STRING_OPT,
                     "The chip-type reported for the run.",
                     "");
}

/**
 * Make sure that our options are sane. Call Err::errAbort if not.
 */
void DmetCopyNumberEngine::checkOptionsImp() {

    defineStates();

    setLibFileOpt("region-model");
    setLibFileOpt("probeset-model");
    setLibFileOpt("cn-region-gt-probeset-file");

    string outDir = getOpt("out-dir");

    // setup output folder
    if(!Fs::isWriteableDir(outDir)) {
        if(Fs::mkdirPath(outDir, false) != APT_OK) {
            Err::errAbort("Can't make or write to directory: " + ToStr(outDir));
        }
    }

    ///@todo does the input summary file chip type match model files chip type
    ///@todo are the model files from the same set
    ///@todo was the input summary file summarized correctly -- ie analysis spec match that in model files

}

