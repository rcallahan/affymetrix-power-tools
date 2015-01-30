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

//
#include "chipstream/MultiChannelNonOverlapCelListener.h"
//
#include "chipstream/GenoUtility.h"
#include "chipstream/apt-geno-qc/GenoQC.h"
#include "stats/stats.h"
#include "util/Fs.h"

//////////

//MultiChannelNonOverlapCelListener::MultiChannelNonOverlapCelListener() {
//  declareMetrics();
//  m_layout=NULL;
//  m_Label = "UNSET";
//}

MultiChannelNonOverlapCelListener::~MultiChannelNonOverlapCelListener()
{
  // no op
}

// Process another cel files worth of data.
// void MultiChannelNonOverlapCelListener::newChip(affymetrix_fusion_io::FusionCELData *cel){
//
//   std::vector<ChipSummary::Metric> metrics;
//   m_SummaryStats.push_back(metrics);
//   m_CelNames.push_back(cel->GetFileName());
//
//   int chipIdx=getNextChipIdx();
//
//   std::string medianDiff_label = m_Label+"_median_diff_over_iqr";
//   std::string at_median_label = m_Label+"_AT_median";
//   std::string gc_median_label = m_Label+"_GC_median";
//
//   // Populate the Summary Metrics
//   setMetric(chipIdx,m_Label,          1.0);
//   setMetric(chipIdx,medianDiff_label, 1.0);
//   setMetric(chipIdx,at_median_label,  1.0);
//   setMetric(chipIdx,gc_median_label,  1.0);
//
//   setValid(true);
// }

void MultiChannelNonOverlapCelListener::declareMetrics()
{
  // These were declared "uninformative" and should no longer be reported.
  //declareMetric(m_Label+"_AT_median",ChipSummary::Metric::Double);
  //declareMetric(m_Label+"_GC_median",ChipSummary::Metric::Double);
  //declareMetric(m_Label+"_median_diff_over_iqr",ChipSummary::Metric::Double);
  declareMetric(m_Label+"_DQC",ChipSummary::Metric::Double);
  declareMetric(m_Label+"_log_diff_qc",ChipSummary::Metric::Double);
  declareMetric("saturation_GC",ChipSummary::Metric::Double);
  declareMetric("saturation_AT",ChipSummary::Metric::Double);
}

MultiChannelNonOverlapCelListener::MultiChannelNonOverlapCelListener(std::vector<ProbeListPacked>& probeSets,
                                                                     QCProbesetOptions& psOpts,
                                                                     ChipLayout* layout,
                                                                     std::string label) :
  m_ProbeSets(probeSets),
  m_psOpts(psOpts),
  m_layout(layout),
  m_Label(label),
  m_cel (NULL)
{
  m_K = 2.0;
  declareMetrics();
  Verbose::out(2, "Initializing Dish QC Okay!");
}

/**
 * Process another cel files worth of data.
 */
void
MultiChannelNonOverlapCelListener::newChip(affymetrix_fusion_io::FusionCELData* cel)
{
  // remember the pointer for later.
  m_cel=cel;
  //
  m_CelNames.push_back(cel->GetFileName());

  // Get Contrast Values
  vector<double> ATcontrastValues;
  vector<double> GCcontrastValues;
  // Get LogDiff values
  vector<double> ATLogDiffValues;
  vector<double> GCLogDiffValues;

  vector<double> ATSaturationValues;
  vector<double> GCSaturationValues;

  // Vincent Bressler 11/15/2012
  //vector<float> intensities;
  //double sat_gc = 0;
  //double sat_at = 0;

  //{
	 // int n = cel->GetNumCells ();
	 // intensities.resize (n);
	 // int check = cel->GetIntensities (0,intensities); // GC is channel 0
	 // int nsat = 0;
	 // for (int i=0; i<n; i++)
	 // {
		//  if (intensities[i] >= 3800)
		//	  nsat++;
	 // }

	 // sat_gc = (double)nsat / (double)n;
  //}

  //{
	 // int n = cel->GetNumCells ();
	 // cel->GetIntensities (1,intensities); // AT is channel 1
	 // int nsat = 0;
	 // for (int i=0; i<n; i++)
	 // {
		//  if (intensities[i] >= 3800)
		//	  nsat++;
	 // }

	 // sat_at = (double)nsat / (double)n;
  //}
  
  //need to edit here to change the contrast
  GimmeContrast(m_ProbeSets,
                ATcontrastValues, GCcontrastValues,
                ATLogDiffValues,  GCLogDiffValues,
				ATSaturationValues, GCSaturationValues,
                m_K);

  Verbose::out(4, "AT nonpolymorphic probeset size is: " + ToStr(ATcontrastValues.size()));
  Verbose::out(4, "GC nonpolymorphic probeset size is: " + ToStr(GCcontrastValues.size()));

  Verbose::out(4, "calculate new log diff score");
  double atLogDiffScore, gcLogDiffScore, logDiffQC ;
  calcLogDiff(ATLogDiffValues, atLogDiffScore);
  calcLogDiff(GCLogDiffValues, gcLogDiffScore);
  logDiffQC = atLogDiffScore + gcLogDiffScore;

  CumulativeStats<double> GCstats;
  CumulativeStats<double> ATstats_full;
  CumulativeStats<double> GCstats_full;

  double nonOverlap = calcNonOverlap(ATcontrastValues, GCcontrastValues, GCstats);
  
  // difference between median / avg IQR
  //double medianDiff = calcMedianDiff(ATcontrastValues, GCcontrastValues, GCstats_full, ATstats_full);
  int chipIdx=getNextChipIdx();
  setMetric(chipIdx,m_Label+"_DQC"                 , nonOverlap);
  //setMetric(chipIdx,m_Label+"_median_diff_over_iqr", medianDiff);
  //setMetric(chipIdx,m_Label+"_AT_median"           , ATstats_full.getMedian());
  //setMetric(chipIdx,m_Label+"_GC_median"           , GCstats_full.getMedian());
  setMetric(chipIdx,m_Label+"_log_diff_qc"         , logDiffQC);
  setMetric(chipIdx,"saturation_GC"       , GCSaturationValues[0]);
  setMetric(chipIdx,"saturation_AT"       , ATSaturationValues[0]);
  setValid(true);
}

///** clear a probeset structure */
//void MultiChannelNonOverlapCelListener::clearProbeSet(ProbeSet &ps) {
//  //check ownMem is true first
//  delete [] ps.name;
//  ps.name = NULL;
//  //ps.atoms.clear();
//};

/**
 * Go thru GC contrast values and find the mean and stdev
 * go thru AT contrast values and count how many are within 2stdev of GC mean
 * calculate the fraction of non-overlap
 * @param ATcontrast Values
 * @param GCcontrast Values
 * @param GCstats vector to fill in with stats
 */
double
MultiChannelNonOverlapCelListener::calcNonOverlap(const std::vector<double> &ATcontrastValues,
                                                  const std::vector<double> &GCcontrastValues,
                                                  CumulativeStats<double>  &GCstats)
{
  for(int gcIx = 0; gcIx < GCcontrastValues.size(); gcIx++) {
    float intensity = GCcontrastValues[gcIx];
    GCstats.addData(intensity);
  }
  float gc_mean = GCstats.getMean();
  float gc_stdev = GCstats.getStdev();

  float gcthreshold = gc_mean + 2*gc_stdev;
  int inGC = 0;

  for (int atIx = 0; atIx< ATcontrastValues.size(); atIx++){
    if (ATcontrastValues[atIx] < gcthreshold){
      inGC++;
    }
  }

  double result = (double) inGC / int(ATcontrastValues.size());
  Verbose::out(4, "in GC fraction is " + ToStr(result));
  result = 1.0-result;
  return result;
  
}


/**
   * Go thru one channel's contrast values and calculate the fraction 
   * of probesets above or under given threshold.
   * @param contrastValues for a given channel
   * @param threshold to decide boundary
   * @para greater: count number of probeset greater than threshold or not
   */
double MultiChannelNonOverlapCelListener::calcCrossThreshold(const std::vector<double> &contrastValues,
double threshold, bool greater ) {
    int cnt = 0;		
    if(greater){
    	for (int ix = 0; ix< contrastValues.size(); ix++){
    	     if (contrastValues[ix] > threshold){
      		cnt++;
     	    }
        }
  }
   else{
        for (int ix = 0; ix< contrastValues.size(); ix++){
           if (contrastValues[ix] < threshold){
      	      cnt++;
    	   }
        }
  }
  
  double result = (double) cnt/ contrastValues.size();
  Verbose::out(4, "cross threshold fraction is " + ToStr(result));
  return result;
}

void MultiChannelNonOverlapCelListener::calcLogDiff(std::vector<double> &logdiffValues, double &avg, double &score, double &IQRS){
  CumulativeStats<double> values;
  values.setFull();
  if(logdiffValues.size() == 0){
	Err::errAbort("Empty logdiffvalues.");
  }
  for(int idx = 0; idx < logdiffValues.size(); idx++) {
    double logdiff = logdiffValues[idx];
    values.addData(logdiff);
  }
  values.setUpFullMethod();
  avg = values.getMean();
  double diff_std = values.getStdev();
  double IQR = values.getIQR();
  if(diff_std == 0 || IQR == 0){
    Err::errAbort("No variance or inter quantile rangefor sample.");
  }
  score = avg/diff_std;
  IQRS = avg/IQR;
  
}

void MultiChannelNonOverlapCelListener::calcLogDiff(std::vector<double> &logdiffValues, double &score){
  CumulativeStats<double> values;
  values.setFull();
  if(logdiffValues.size() == 0){
	Err::errAbort("Empty logdiffvalues.");
  }
  for(int idx = 0; idx < logdiffValues.size(); idx++) {
    double logdiff = logdiffValues[idx];
    values.addData(logdiff);
  }
  values.setUpFullMethod();
  double avg = values.getMean();
  double diff_std = values.getStdev();
  if(diff_std == 0){
    Verbose::out(4, "No variance range for sample. Something wrong, please check.");
    score=-100;
  }
  else{
    score = avg/diff_std;
  }
}


/**
   * Go thru GC contrast values and find the median and IQR
   * go thru AT contrast values and find the median and IQR
 * calculate the abs(GC_median - AT_median) / average_IQR
 * @param ATcontrast Values
 * @param GCcontrast Values
 *
 */
double
MultiChannelNonOverlapCelListener::calcMedianDiff(const std::vector<double> &ATcontrastValues,
                                                  const std::vector<double> &GCcontrastValues,
                                                  CumulativeStats<double> &ATstats,
                                                  CumulativeStats<double> &GCstats)
{
  GCstats.setFull();
  ATstats.setFull(); //use full method

  for(int gcIx = 0; gcIx < GCcontrastValues.size(); gcIx++) {
    float intensity = GCcontrastValues[gcIx];
    GCstats.addData(intensity);
  }
  for(int atIx = 0; atIx < ATcontrastValues.size(); atIx++) {
    float intensity = ATcontrastValues[atIx];
    ATstats.addData(intensity);
  }
  GCstats.setUpFullMethod();
  ATstats.setUpFullMethod();

  double gc_median = GCstats.getMedian();
  double gc_IQR = GCstats.getIQR();
  double at_median = ATstats.getMedian();
  double at_IQR = ATstats.getIQR();
  double result = abs (gc_median - at_median) / ((gc_IQR + at_IQR)/2);
  //double gc_mean = GCstats.getMean();
  /*Verbose::out(4, "gc_median is " + ToStr(gc_median));
    Verbose::out(4, "at_median is " + ToStr(at_median));
    Verbose::out(4, "gc_IQR is " + ToStr(gc_IQR));
    Verbose::out(4, "at_IQR is " + ToStr(at_IQR));
    Verbose::out(4, "gc_mean is " + ToStr(gc_mean));

    Verbose::out(4, "medianDiff is " + ToStr(result));
  */
  return result;


}


/**
 * Loop through the probesets provided and calculate a contrast value
 * for each one using the median of PM probes for A allele and B
 * allele.
 */
///@todo refactor? simmilar code in EmGenderCelListener
void MultiChannelNonOverlapCelListener::fillInContrastValues(int celIx,
                                                             std::vector<ProbeListPacked>& probeSets,
                                                             std::vector<double> &contrastValues,
                                                             double k)
{
  contrastValues.clear();
  contrastValues.reserve(probeSets.size());
  /// @todo why is this done?
}

/**
 * Add a new probeset subset via a mask associated with identifier name.
 *
 * @param name - string identifier associated with this mask (i.e. "all", "pm")
 * @param mask - bit mask with probesets to calculate stats for set to true.
 */
void MultiChannelNonOverlapCelListener::addProbeMask(const std::string &name, const std::vector<bool> &mask)
{
  m_MaskNames.push_back(name);
  m_ProbesetMasks.push_back(mask);
}

void MultiChannelNonOverlapCelListener::readSampleAlleleMap(const std::string& filename,
                                                            std::map<std::string,std::vector<int> >& sampleAlleleMap)
{
  affx::TsvFile tsv;
  Verbose::out(1, "reading channel-file: "+filename);

  std::string probeset;
  int alleleA, alleleB;

  tsv.bind(0, "probeset_name",    &probeset, affx::TSV_BIND_REQUIRED);
  tsv.bind(0, "allele_A_channel", &alleleA,  affx::TSV_BIND_REQUIRED);
  tsv.bind(0, "allele_B_channel", &alleleB,  affx::TSV_BIND_REQUIRED);

  if (tsv.open(filename)!=affx::TSV_OK) {
    APT_ERR_ABORT("Could not open channel-file: '" + filename + "'");
  }

  while (tsv.nextLevel(0)==affx::TSV_OK) {
    //printf("channel: '%s' %d %d\n",probeset.c_str(),alleleA,alleleB);
    sampleAlleleMap[probeset].push_back(alleleA);
    sampleAlleleMap[probeset].push_back(alleleB);
  }
  tsv.close();
  printf("channel: read %d items\n",int(sampleAlleleMap.size()));
}

void
MultiChannelNonOverlapCelListener::GimmeContrast(std::vector<ProbeListPacked>& probeSets,
                                                 std::vector<double>& ATcontrastValues,
                                                 std::vector<double>& GCcontrastValues,
                                                 std::vector<double> &ATLogDiffValues, 
                                                 std::vector<double> &GCLogDiffValues,
												 std::vector<double> &ATSaturationValues,
												 std::vector<double> &GCSaturationValues,
                                                 double k)
{
  std::string fileName=m_cel->GetFileName();
  std::string fileNameRoot=Fs::basename(fileName);

  Verbose::out(1, "computing contrast for '"+fileName+"' in MultiChannelHomHiLoCelListener::GimmeContrast");

  ATcontrastValues.clear();
  GCcontrastValues.clear();
  ATLogDiffValues.clear();
  GCLogDiffValues.clear();

  // Construct the analysis streams for each factory
  AnalysisStreamFactory asFactory;
  AnalysisStream*       asStream;

  ///
  std::map<std::string,std::vector<int> > sampleAlleleMap;
  if (m_psOpts.m_channelFile!="") {
    readSampleAlleleMap(m_psOpts.m_channelFile,sampleAlleleMap);
  }
  else {
    //Verbose::out(1,"No channel file. Assuming allele A on channel A and allele B on channel B.");
    Verbose::out(1,"No channel file. Taking channel info from probesets.");
  }

  if (m_psOpts.m_sketchInFile!="") {
    asFactory.readTargetSketchFromFile(m_psOpts.m_sketchInFile);
  }
  asStream=asFactory.constructAnalysisStream(m_psOpts.m_exprAnalysisString,*m_layout);

  std::vector<std::string> celFiles;
  celFiles.push_back(fileName);

  SparseMart* sMart=NULL;
  sMart=new SparseMart(m_layout->getProbeLayoutOrder(),
                              celFiles,
                              true);

  CelReader cReader;
  cReader.registerStream(asStream->getChipStreamHead());
  cReader.registerIntensityMart(sMart);

  // Setup the MultiChannelQuantMethod
  MultiQuantMethodListener multiQuant;
  QuantMethodReportListener* listener = new QuantMethodReportListener();
  listener->registerListener(&multiQuant);
  asStream->addReporter(listener);

  // Dummies needed by the reporter
  PmOnlyAdjust pmAdjust;
  //ChipStream * iTrans;

  //populate Chipstreams and pmAdjusters by taking advantage of analysisStream object
  std::vector<ChipStream*>* cStream;
  PmAdjuster* pmAdjuster;
  cStream=asStream->getChipStream();
  pmAdjuster=asStream->getPmAdjuster();

  // Read in cel files
  // we could just call readFiles() on each cReader, but if there are
  // cross dependencies between the channels, this will not work. For example
  // if there is a chip stream node in one channel which corrects for cross
  // talk and thus needs values from an earlier step in all the other channels */
  cReader.setFiles(celFiles);
  cReader.readFiles();

  GCSaturationValues.resize (celFiles.size ());
  ATSaturationValues.resize (celFiles.size ());
  for (int ifile=0; ifile<celFiles.size (); ifile++)
  {
	  GCSaturationValues[ifile] = cReader.GetSaturationValue (ifile*2);
	  ATSaturationValues[ifile] = cReader.GetSaturationValue (ifile*2+1);
  }

  // Loop through all of our probe sets and do an expression quantification for each one.
  // @todo: perhaps here is a place to kick in the QCprobesetOptions

  // Verbose::progressBegin(1,"Processing probesets",
  // 40, (int)(layout.getProbeSetCount()/40), layout.getProbeSetCount());

  ProbeSet m_Aallele;
  ProbeSet m_Ballele;

  //int jhg_cnt=0;

  for (int psIx = 0; psIx < probeSets.size(); psIx++) {
    // Verbose::progressStep(1);
    ProbeSet* ps = ProbeListFactory::asProbeSet(probeSets[psIx]);
    ProbeSetGroup group(ps);     // psGroup should delete the memory for ps...

    // Verbose::out(4, "processing probeset: " + ToStr(ps->name));
    asStream->doAnalysis(group, *sMart, true);

    fillInAlleleProbeSets(*ps, m_Aallele, m_Ballele);

    ///@todo much of this should probably be in a MultiQuantLabelZ
    ///      rather than using QuantLabelZ directly. Accessing
    ///      m_Results directly is ugly as well.

    // Figure out the channels for this SNP
    // @todo default should be global (ie if we have a channel file,
    // then all probesets must be found in it. only default to 0/1 when there is no channel file)
    //int alleleAindex = 0;
    //int alleleBindex = 1;
//    if (sampleAlleleMap.find(ps->name) != sampleAlleleMap.end()) {
//      alleleAindex = sampleAlleleMap[ps->name][0];
//      alleleBindex = sampleAlleleMap[ps->name][1];
//    }
//    else if () {
//    }
//    else {
//      alleleAindex = 0;
//      alleleBindex = 1;
//      Verbose::out(1,"Unable to find probeset '" + ToStr(ps->name) + "'"
//                   " in channel file. Assuming channels 1 and 2.");
//    }

    bool success = true;
    std::vector<double> Asignals;
    std::vector<double> Bsignals;

    if ((asStream->getQuantMethod()) == NULL) {
      Err::errAbort("missing quantMethod...");
    }

    // A
    success=summarizeAllele(&m_Aallele,
                            Asignals,
                            *sMart,
                            cStream,
                            *pmAdjuster,
                            asStream->getQuantMethod(),
                            true);
    if (!success) {
      Err::errAbort("problem with at MultiChannelNonOverlapLoCelListener::summarizeAllele. "
                    "Cannot succesfully summarize value for A allele.");
    }
    // B
    success = summarizeAllele(&m_Ballele,
                              Bsignals,
                              *sMart,
                              cStream,
                              *pmAdjuster,
                              asStream->getQuantMethod(),
                              true);
    if (!success) {
      Err::errAbort("problem with at MultiChannelNonOverlapCelListener::summarizeAllele. "
                    "Cannot succesfully summarize value for B allele.");
    }

    //
    double contrast = -2, strength = -2;
    ContrastExtremesStretch(k, Asignals[0], Bsignals[0], contrast, strength);
    // iterate thru the number of masks, usually just AT and GC
    for (int m = 0; m<m_ProbesetMasks.size(); m++){
      string probesetclass = m_MaskNames[m];
      double logDiff=0.0;
      if (m_ProbesetMasks[m][psIx] == true && probesetclass.compare("AT") == 0){
        ATcontrastValues.push_back(contrast);
        logDiff=log(Asignals[0]/Bsignals[0]); 
        ATLogDiffValues.push_back(logDiff);
      }
      else if (m_ProbesetMasks[m][psIx] == true && probesetclass.compare("GC") == 0){
        GCcontrastValues.push_back(contrast);
        logDiff =  log(Bsignals[0]/Asignals[0]);                                                              
        GCLogDiffValues.push_back(logDiff);                                                                   
      }
   }
    //Verbose::out(1,"SIGNALS:"+ToStr(ps->name)+":"+
    //ToStr(Asignals[0])+":"+ToStr(Bsignals[0])+":"+ToStr(contrast));
    //
    // clear out results for the current probeset
    m_Aallele.clear();
    m_Ballele.clear();
    multiQuant.m_Results.clear();

    //
    //if (jhg_cnt++>100) {
    //  break;
    //}
  }
  //  Verbose::progressEnd(1, "Done.");
  
  //
  delete asStream;
  delete sMart;
}

bool
MultiChannelNonOverlapCelListener::summarizeAllele(ProbeSet* pSet,
                                                   vector<double>& summaryValues,
                                                   SparseMart& sMart,
                                                   std::vector<ChipStream *>* iTrans,
                                                   PmAdjuster& pmAdjust,
                                                   QuantMethod* quantMethod,
                                                   bool lowPrecision)
{
  bool success = true;
  QuantExprMethod *qeMethod = dynamic_cast<QuantExprMethod *>(quantMethod);
  summaryValues.clear();
  /* This is a little interesting as we are spoofing a probeSetGroup for calling
     functions but don't want the destructor called while that group thinks it
     owns those probesets. */
  ProbeSetGroup group(pSet);
  //int probesize = group.countPmProbes();
  //int probesetsize = group.probeSets.size();
  //Verbose::out(4, "countPmProbes() is " + ToStr(probesize));
  //Verbose::out(4, "group.probeSets.size() is " + ToStr(probesetsize));
  if(!qeMethod->setUp(group, sMart, *iTrans, pmAdjust))
    success = false;
  else {
    qeMethod->computeEstimate();

    for(unsigned int i = 0; i < qeMethod->getNumTargets(); i++) { //for each metaCelFile
      double val = qeMethod->getTargetEffect(i); //get signal intensities of the pSet
      if(lowPrecision) {
        val = Util::round(val);
      }
      summaryValues.push_back(val);
    }
    success = true;
  }

  /* clear out probesets to prevent deleting them. */
  group.probeSets.clear();
  //clearProbeSet(*pSet);
  pSet->clear();
  return success;
}

/* stolen from QuantLabelZ, QuantGTyoeMethod */
bool
MultiChannelNonOverlapCelListener::fillInAlleleProbeSets(const ProbeSet &gtPs,
                                                         ProbeSet &aAllele,
                                                         ProbeSet &bAllele)
{
  aAllele.clear();
  bAllele.clear();
    
  // First set up the probe sets to do the summarization of a and b alleles.
  aAllele.psType = ProbeSet::Expression;
  bAllele.psType = ProbeSet::Expression;

  // Set the names for our spoofed probesets with suffix A and B
  string ps_name = gtPs.name;
  aAllele.name = Util::cloneString(ps_name.c_str());
  bAllele.name = Util::cloneString(ps_name.c_str());

  //int groupOffset = 0; // atoms are stored in blocks for each allele, keep track of which block we're in.
  /*  Verbose::out(4, "ProbeSet::Aforward Areverse Bforward Breverse are "
      + ToStr(ProbeSet::Aforward) + ", " + ToStr(ProbeSet::Areverse) + ", " + ToStr(ProbeSet::Bforward)
      + ", " +ToStr(ProbeSet::Breverse));
  */
  //Verbose::out(4, "numGroups for " + ToStr(gtPs.name) + " is " + ToStr(gtPs.numGroups));
  //for(unsigned int groupIx = 0; groupIx < gtPs.numGroups; groupIx++) {
  //Verbose::out(4, "atomsPerGroup for " + ToStr(groupIx) + " is " + ToStr(gtPs.atomsPerGroup[groupIx]));
  //assuming one one group

  // single atom from an old spf file.
  if (gtPs.atoms.size()!=1) {
    APT_ERR_ABORT("The probeset '"+ps_name+"' has "+ToStr(gtPs.atoms.size())+" blocks. (should be 1)");
  }
  for(unsigned int atomIx = 0; atomIx < gtPs.atoms.size(); atomIx++) {
    // Duplicate the single block onto the two alleles
    // use Atom::deepCopy to avoid messing with the memory management of ProbeSets.
    Atom* a;
    //
    a=Atom::deepCopy(*gtPs.atoms[atomIx]);
    aAllele.atoms.push_back(a);
    aAllele.atoms[atomIx]->channel_code=1;
    a=NULL;
    //
    a=Atom::deepCopy(*gtPs.atoms[atomIx]);
    bAllele.atoms.push_back(a);
    bAllele.atoms[atomIx]->channel_code=0;
    a=NULL;
  }
  
  return true;
}
//  // dual atom from a newer spf file...
//  else if (gtPs.atoms.size()==2) {
//    // walk the atoms, putting them on the correct alleles.
//    for(unsigned int atomIx = 0; atomIx < gtPs.atoms.size(); atomIx++) {
//      // split the blocks onto seperate probesets.
//      Atom* a=gtPs.atoms[atomIx];
//      if (a->channel_code==1) {
//        aAllele.atoms.push_back(Atom::deepCopy(*a));
//      }
//      else if (a->channel_code==0) {
//        bAllele.atoms.push_back(Atom::deepCopy(*a));
//      }
//      else {
//        APT_ERR_ABORT("Bad channel code: '"+ToStr(a->channel_code)+"'  Should be '0' or '1'.");
//      }
//    }
//  }
//  // whoops!
//  else {
//    APT_ERR_ABORT("Bad number of blocks: '"+ToStr(gtPs.atoms.size())+"'  Should be '1' or '2'.");
//  }
//    
//  return true;
//}
