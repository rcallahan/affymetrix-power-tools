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
#include "chipstream/MultiChannelHomHiLoCelListener.h"
//
#include "chipstream/GenoUtility.h"
#include "chipstream/apt-geno-qc/GenoQC.h"
#include "stats/stats.h"
#include "util/Fs.h"

//using namespace std;
using namespace affx;

//MultiChannelHomHiLoCelListener::MultiChannelHomHiLoCelListener() {
//  declareMetrics();
//  m_layout=NULL;
//}

MultiChannelHomHiLoCelListener::MultiChannelHomHiLoCelListener(std::vector<ProbeListPacked>& probeSets,
                                                               QCProbesetOptions& psOpts,
                                                               ChipLayout* layout,
                                                               double k,
                                                               double emThresh,
                                                               double binSize,
                                                               std::string label) :
  m_ProbeSets(probeSets),
  m_psOpts(psOpts),
  m_layout(layout),
  m_K(k),
  m_EmThresh(emThresh),
  m_BinSize(binSize),
  m_Label(label)
{
  declareMetrics();
  Verbose::out(2, "Initializing Multichannel HomHiLoCelListener Okay!");
  m_cel = NULL;
}

MultiChannelHomHiLoCelListener::~MultiChannelHomHiLoCelListener()
{
  // nothing.
}

void MultiChannelHomHiLoCelListener::declareMetrics()
{
  // the reference version has "cqc" as part of the name.
  // what does "cqc" mean?
  declareMetric(m_Label+"-minhilo",ChipSummary::Metric::Double);
  declareMetric(m_Label+"-maxhilo",ChipSummary::Metric::Double);
  declareMetric(m_Label+"-aa-diff",ChipSummary::Metric::Double);
  declareMetric(m_Label+"-bb-diff",ChipSummary::Metric::Double);
  declareMetric(m_Label+"-bb-aa-diff",ChipSummary::Metric::Double);
}

/** clear a probeset structure */
void MultiChannelHomHiLoCelListener::clearProbeSet(ProbeSet &ps) {
  //check ownMem is true first
  delete [] ps.name;
  ps.name = NULL;
  ps.atoms.clear();
};

/**
 * Process another cel files worth of data.
 */
void MultiChannelHomHiLoCelListener::newChip(affymetrix_fusion_io::FusionCELData* cel) 
{
  std::string fileRoot=Fs::basename(cel->GetFileName());

  // Get Contrast Values
  vector<double> contrastValues;

  // grab copy of the pointer for tmp use.
  m_cel=cel;
  // need to edit here to change the contrast
  GimmeContrast(m_ProbeSets, contrastValues, m_K);

  Verbose::out(4, "contrastValues size is: " + ToStr(contrastValues.size()));
  // Find the 3 peaks in contrast values (AA, AB, BB)
  double peak1, peak2, peak3;
  computeClusterPeaks(contrastValues,peak1,peak2,peak3,m_EmThresh);
  Verbose::out(4,"MultiChannelHomHiLo - " + fileRoot + " -  " + m_Label +
               " - Peaks: " + ToStr(peak1) + ", " + ToStr(peak2) + ", " + ToStr(peak3));

  // Find the valleys
  double valley1 = computeClusterValley(contrastValues,peak1,peak2,m_BinSize);
  double valley2 = computeClusterValley(contrastValues,peak2,peak3,m_BinSize);
  Verbose::out(4,"MultiChannelHomHiLo - " + fileRoot + " -  " + m_Label +
               " - Valleys: " + ToStr(valley1) + ", " + ToStr(valley2));

  // How many contrast values under each peak
  double densityPeak1 = getContrastDensity(contrastValues,peak1,m_BinSize);
  double densityPeak3 = getContrastDensity(contrastValues,peak3,m_BinSize);
  Verbose::out(4,"MultiChannelHomHiLo - " + fileRoot + " -  " + m_Label +
               " - Peak Density: " + ToStr(densityPeak1) + ", " + ToStr(densityPeak3));

  // How many contrast values under each valley
  double densityValley1 = getContrastDensity(contrastValues,valley1,m_BinSize);
  double densityValley2 = getContrastDensity(contrastValues,valley2,m_BinSize);
  Verbose::out(4,"MultiChannelHomHiLo - " + fileRoot + " -  " + m_Label +
               " - Valley Density: " + ToStr(densityValley1) + ", " + ToStr(densityValley2));

  double minHomHiLo = 0.1;
  double maxHomHiLo;
  double aa_diff, bb_diff, bb_aa_diff;

  string minhilo_label    = m_Label+"-minhilo";
  string maxhilo_label    = m_Label+"-maxhilo";
  string aa_diff_label    = m_Label+"-aa-diff";
  string bb_diff_label    = m_Label+"-bb-diff";
  string bb_aa_diff_label = m_Label+"-bb-aa-diff";

  double mn1 = densityPeak1 - densityValley1;
  double mn2 = densityPeak3 - densityValley2;

  // Get Min
  if (mn1 < mn2)
    minHomHiLo = mn1;
  else
    minHomHiLo = mn2;

  Verbose::out(4,"MultichannelHomHiLo - " + fileRoot + " -  " + m_Label +
               " - mn1,mn2,min: " + ToStr(mn1) + ", " + ToStr(mn2) + ", " + ToStr(minHomHiLo));

  int chipIdx=getNextChipIdx();

  if (minHomHiLo == mn1)
    maxHomHiLo = mn2;
  else
    maxHomHiLo = mn1;

  setMetric(chipIdx,minhilo_label, minHomHiLo);
  setMetric(chipIdx,maxhilo_label, maxHomHiLo);

  if ((peak1 < peak2) && (peak2 < peak3)) {
    bb_diff = densityPeak1 - densityValley1;
    aa_diff = densityPeak3 - densityValley2;
    bb_aa_diff = bb_diff - aa_diff;
    setMetric(chipIdx,aa_diff_label, aa_diff);
    setMetric(chipIdx,bb_diff_label, bb_diff);
    setMetric(chipIdx,bb_aa_diff_label, bb_aa_diff);
  }
  else {
    Err::errAbort("MultiChannelHomHiLoCelListener::newMultiChannelChip() peaks are not in order: " +
                  ToStr(peak1) + " " + ToStr(peak2) + " " + ToStr(peak3));
  }
  //
  setValid(true);
}

/**
 * Loop through the probesets provided and calculate a contrast value
 * for each one using the median of PM probes for A allele and B
 * allele.
 */
///@todo refactor? simmilar code in EmGenderCelListener
void MultiChannelHomHiLoCelListener::fillInContrastValues(int celIx,
                                                          std::vector<ProbeList> &probeSets,
                                                          std::vector<double> &contrastValues, double k) {
  contrastValues.clear();
  contrastValues.reserve(probeSets.size());

  //for(int psIx = 0; psIx < probeSets.size(); psIx++) {
  //const ProbeSet *ps =  ProbeListFactory::asProbeSet(probeSets[psIx]);
  //contrastValues.push_back(CalculateEMContrast(celA, celB, ps, k));
  //delete ps;
  //}
}

/**
 * Calculate contrast values
 */
///@todo refactor? simmilar code in EmGenderCelListener
/*
  double HomHiLoCelListener::CalculateEMContrast(FusionCELData* celA, FusionCELData* celB, const ProbeSet* ps, const double k) {
  float Amedian = -1, Bmedian = -1;
  bool medianOk = summarize(o, celA, celB, ps, Amedian, Bmedian);
  if(!medianOk || Amedian <= 0 || Bmedian <= 0) {
  Err::errAbort("MultiChannelHomHiLoCelListener::fillInContrastValues() - Warning. alleleMedian() failed for probeset: " +
  ToStr(ps->name));
  }
  double contrast = -2, strength = -2;
  ContrastExtremesStretch(k, Amedian, Bmedian, contrast, strength);

  if(fabs(contrast) > 1) {
  Err::errAbort("MultiChannelHomHiLoCelListener::fillInContrastValues() - Can't have abs(contrast) > 1 for probeset: " +
  ToStr(ps->name));
  }
  return contrast;
  }
*/
/**
 * Compute the cluster peaks using EM
 */
void MultiChannelHomHiLoCelListener::computeClusterPeaks(const std::vector<double>& contrastValues,
                                                         double &peak1,
                                                         double &peak2,
                                                         double &peak3,
                                                         double emThresh)
{
  CEMSeed seed;
  seed.setMu(-0.66f, 0.0f, 0.66f);
  seed.setSigma(0.1f, 0.1f, 0.1f);
  seed.setWeight(0.33f, 0.34f, 0.33f);
  seed.setMinMu(-2.0f, -.05f, .25f);
  seed.setMinSigma(0.02f, 0.02f, 0.02f);
  seed.setMaxMu(-.25f, .05f, 2.0f);
  seed.setMaxSigma(.3f, .3f, .3f);

  CPrimeEM em;
  em.setData(contrastValues);
  em.setThreshold(emThresh);
  em.EMEstimate(seed);
  CEMEst* pEst = em.getEMEstimates();

  peak1 = pEst->m_mu[0];
  peak2 = pEst->m_mu[1];
  peak3 = pEst->m_mu[2];
}

/**
 * Compute the Valley Between Two Points
 */
double MultiChannelHomHiLoCelListener::computeClusterValley(const std::vector<double>& contrast,
                                                            double x1,
                                                            double x2,
                                                            double bin){
  double valleyX = x1;
  double valleyDensity = getContrastDensity(contrast, valleyX, bin);
  double x = x1 + 2.0 * bin;
  while (x < x2) {
    double density = getContrastDensity(contrast, x, bin);
    if (density < valleyDensity){
      valleyX = x;
      valleyDensity = density;
    }
    x += 2.0 * bin;
  }
  return valleyX;
}

/**
 * Compute the density of values for a given bin
 */
double MultiChannelHomHiLoCelListener::getContrastDensity(const std::vector<double>& contrast, double x, double bin){
  double lb = x - bin;
  double ub = x + bin;
  size_t n = contrast.size();
  int count = 0;
  for (size_t i=0; i < n; i++) {
    if ( (contrast[i]> lb) && (contrast[i] <= ub))
      count++;
  }
  return 100.0*(double(count)/double(n));
}

void MultiChannelHomHiLoCelListener::GimmeContrast(std::vector<ProbeListPacked>& probeSets,
                                                   std::vector<double>& contrastValues,
                                                   double k)
{
  // "m_cel" is a pointer to the celfile we are working on.
  // this assumes that the AT channel is 1 and GC is 0

  std::string fileName=m_cel->GetFileName();
  std::string fileNameRoot=Fs::basename(fileName);

  Verbose::out(1, "computing contrast for '"+fileNameRoot+"' in MultiChannelHomHiLoCelListener::GimmeContrast");
  contrastValues.clear();
  contrastValues.reserve(probeSets.size());

  // Construct the analysis streams for each factory
  AnalysisStreamFactory asFactory;
  AnalysisStream* asStream;

  Verbose::out(1,"Allele estimates computed using: " + ToStr(m_psOpts.m_exprAnalysisString));

  //
  if(m_psOpts.m_sketchInFile!="") {
    asFactory.readTargetSketchFromFile(m_psOpts.m_sketchInFile);
  }
  asStream=asFactory.constructAnalysisStream(m_psOpts.m_exprAnalysisString,*m_layout);
  
  // Setup Intensity Mart for each channel
  std::vector<std::string> celFiles;
  celFiles.push_back(fileName);
  
  SparseMart* sMart=NULL;
  sMart = new SparseMart(m_layout->getProbeLayoutOrder(),
                                celFiles, 
                                true);
  
  // Setup CelReaders for each Channel
  CelReader cReader;
  cReader.registerStream(asStream->getChipStreamHead());
  cReader.registerIntensityMart(sMart);

  // Setup the MultiChannelQuantMethod
  MultiQuantMethodListener multiQuant;
  QuantMethodReportListener *listener = new QuantMethodReportListener();
  listener->registerListener(&multiQuant);
  asStream->addReporter(listener);

  // Dummies needed by the reporter
  PmOnlyAdjust pmAdjust;
  //ChipStream* iTran;

  //populate Chipstreams and pmAdjusters by taking advantage of analysisStream object
  std::vector<ChipStream*>* cStream; //now vector of vector of chip streams
  PmAdjuster* pmAdjuster; //now a vector of PmAdjuster
  cStream=asStream->getChipStream();
  pmAdjuster=asStream->getPmAdjuster();

  // register the cel files to open.
  Verbose::out(1,"Registering CEL files for each channel.");

  // Read in cel files
  // we could just call readFiles() on each cReader, but if there are
  // cross dependencies between the channels, this will not work. For example
  // if there is a chip stream node in one channel which corrects for cross
  // talk and thus needs values from an earlier step in all the other channels
  cReader.setFiles(celFiles);
  cReader.readFiles();

  //
  // Loop through all of our probe sets and do an expression quantification for each one.
  // @todo: perhaps here is a place to kick in the QCprobesetOptions

  for(unsigned int psIx = 0; psIx < probeSets.size(); psIx++) {
    // Verbose::progressStep(1);
    //ProbeList pList = layout.getProbeListAtIndex(psIx);
    ProbeSet* ps = ProbeListFactory::asProbeSet(probeSets[psIx]);
    ProbeSetGroup group(ps);     // psGroup should delete the memory for ps...

    // Verbose::out(4, "processing probeset: " + ToStr(ps->name));
    asStream->doAnalysis(group, *sMart, true);

    ProbeSet m_Aallele;
    ProbeSet m_Ballele;
    fillInAlleleProbeSets(*ps, m_Aallele, m_Ballele);

    //fillInAlleleProbeSets(*ps, m_Aallele, m_Ballele);

    ///@todo much of this should probably be in a MultiQuantLabelZ
    ///      rather than using QuantLabelZ directly. Accessing
    ///      m_Results directly is ugly as well.

    // Figure out the channels for this SNP
    ///@todo default should be global (ie if we have a channel file, then all probesets must be found in it. only default to 0/1 when there is no channel file)
    //int alleleAindex = 0;
    //int alleleBindex = 1;

    //m_Aallele.clear();
    //m_Ballele.clear();

//    if(sampleAlleleMap.find(ps->name) != sampleAlleleMap.end()) {
//      alleleAindex = sampleAlleleMap[ps->name][0];
//      alleleBindex = sampleAlleleMap[ps->name][1];
//    } else {
//      Verbose::out(1,"Unable to find probeset " + ToStr(ps->name) + " in channel file. Assuming channels 1 and 2.");
//    }
    bool success = true;
    std::vector<double> Asignals;
    std::vector<double> Bsignals;

    if((asStream->getQuantMethod()) == NULL ) {
      Err::errAbort("missing quantMethod...");
    }
    // Verbose::out(4, "quantMethod is A:" + ToStr(streams[alleleAindex]->getQuantMethod()));
    // Verbose::out(4, "quantMethod is A:" + ToStr(streams[alleleBindex]->getQuantMethod()));

    //success = summarizeAllele (&m_Aallele, Asignals, *layout,
    // *(iMarts[alleleAindex]), cStreamVec[alleleAindex], pmAdjust, (streams[alleleAindex]->getQuantMethod()), true);
    success = summarizeAllele(&m_Aallele,
                              Asignals,
                              *m_layout,
                              *sMart,
                              cStream,
                              *pmAdjuster,
                              asStream->getQuantMethod(),
                              true);
    // Verbose::out(4, "Asignal is " + ToStr(Asignals[0]));
    if (!success) {
      Err::errAbort("problem with at MultiChannelHomHiLoCelListener::summarizeAllele. "
                    "Cannot succesfully summarize value for A allele.");
    }

    //success = summarizeAllele (&m_Ballele, Bsignals, *layout,
    // *(iMarts[alleleBindex]), cStreamVec[alleleBindex], pmAdjust, (streams[alleleBindex]->getQuantMethod()),true);
    success = summarizeAllele(&m_Ballele,
                              Bsignals,
                              *m_layout,
                              *sMart,
                              cStream,
                              *pmAdjuster,
                              asStream->getQuantMethod(),
                              true);
    if (!success) {
      Err::errAbort("problem with at MultiChannelHomHiLoCelListener::summarizeAllele. "
                    "Cannot succesfully summarize value for B allele.");
    }

    double contrast = -2, strength = -2;
    ContrastExtremesStretch(k, Asignals[0], Bsignals[0], contrast, strength);

    //Verbose::out(1, ToStr(ps->name) + ":" +ToStr(Asignals[0]) + ":" + ToStr(Bsignals[0]) + ":" + ToStr(contrast));
    //Verbose::out(4, ToStr(ps->name) + "-B:" + ToStr(Bsignals[0]));
    //Verbose::out(4, ToStr(ps->name) + "-contrast:" + ToStr(contrast));
    contrastValues.push_back(contrast);
    //Verbose::out(4, ToStr(ps->name)+": "+ToStr(contrastValues[psIx]));

    //
    m_Aallele.clear();
    m_Ballele.clear();

    // clear out results for the current probeset
    multiQuant.m_Results.clear();
  }
  //Verbose::progressEnd(1, "Done.");

  // Free up memory
  delete asStream;
  delete sMart;
}

/** summarize an allele using expression summary methods specified*/
/*no report given, hacked from QuantGTypedMethod.cpp summarizedAllele*/
//returns the median summary of each probeset
//Here summaryValues will be the signal intensity for each sample

bool MultiChannelHomHiLoCelListener::summarizeAllele(ProbeSet* pSet,
                                                     vector<double>& summaryValues,
                                                     ChipLayout& layout,
                                                     SparseMart& sMart,
                                                     std::vector<ChipStream*>* iTrans,
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
  if (!qeMethod->setUp(group, sMart, *iTrans, pmAdjust)) {
    success = false;
  }
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
  clearProbeSet(*pSet);
  return success;
}
 
/* stolen from QuantLabelZ, QuantGTyoeMethod */
bool MultiChannelHomHiLoCelListener::fillInAlleleProbeSets(const ProbeSet &gtPs, ProbeSet &aAllele, ProbeSet &bAllele) {
  /* First set up the probe sets to do the summarization of a and b alleles. */
  aAllele.psType = ProbeSet::Expression;
  bAllele.psType = ProbeSet::Expression;
  // Set the names for our spoofed probesets with suffix A and B
  string name = gtPs.name;
  //name += "-A";
  aAllele.name = Util::cloneString(name.c_str());
  name = gtPs.name;
  //name += "-B";
  bAllele.name = Util::cloneString(name.c_str());
  //int groupOffset = 0; // atoms are stored in blocks for each allele, keep track of which block we're in.
  /*  Verbose::out(4, "ProbeSet::Aforward Areverse Bforward Breverse are "
      + ToStr(ProbeSet::Aforward) + ", " + ToStr(ProbeSet::Areverse) + ", " + ToStr(ProbeSet::Bforward)
      + ", " +ToStr(ProbeSet::Breverse));
  */
  //Verbose::out(4, "numGroups for " + ToStr(gtPs.name) + " is " + ToStr(gtPs.numGroups));
  //for(unsigned int groupIx = 0; groupIx < gtPs.numGroups; groupIx++) {
  //Verbose::out(4, "atomsPerGroup for " + ToStr(groupIx) + " is " + ToStr(gtPs.atomsPerGroup[groupIx]));
  //assuming one one group

  for (unsigned int atomIx = 0; atomIx < gtPs.atoms.size(); atomIx++) {
    Atom* a;
    //aAllele.atoms.push_back(gtPs.atoms[atomIx]);
    a=Atom::deepCopy(*gtPs.atoms[atomIx]);
    a->channel_code=0;
    aAllele.atoms.push_back(a);
    //bAllele.atoms.push_back(gtPs.atoms[atomIx]);
    a=Atom::deepCopy(*gtPs.atoms[atomIx]);
    a->channel_code=1;
    bAllele.atoms.push_back(a);
  }
  
  return true;
}
