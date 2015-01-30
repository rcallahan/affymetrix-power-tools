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
#include "chipstream/apt-geno-qc/GenoQC.h"
//
#include "chipstream/BioTypes.h"
#include "chipstream/CelStatListener.h"
#include "chipstream/ChipLayout.h"
#include "chipstream/CnProbeGenderCelListener.h"
#include "chipstream/DmListener.h"
#include "chipstream/EmGenderCelListener.h"
#include "chipstream/EngineUtil.h"
#include "chipstream/KitAOCelListener.h"
#include "chipstream/ProbeListFactory.h"
#include "chipstream/QuantMethodRunReport.h"
//
#include "broadutil/BroadException.h"
#include "calvin_files/exception/src/ExceptionBase.h"
#include "calvin_files/utils/src/StringUtils.h"
#include "dm/DM.h"
#include "file/TsvFile/TsvFile.h"
#include "portability/affy-base-types.h"
#include "stats/stats.h"
#include "util/Err.h"
#include "util/Fs.h"
#include "util/Util.h"
#include "util/Verbose.h"
//
#include <cassert>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>

//using namespace affymetrix_calvin_exceptions;
//using namespace std;
//using namespace affx;

GenoQC::Reg GenoQC::reg;

GenoQC * GenoQC::FromBase(BaseEngine *engine)
{
	if (engine != NULL && engine->getEngineName() == GenoQC::EngineName())
		return (GenoQC *)engine;
	return NULL;
}

// constructor
GenoQC::GenoQC(ChipLayout *chipLayout) {
  m_kitao_db=NULL;

  defineOptions();

  if(chipLayout != NULL) {
      m_FreeChipLayout = false;
      m_ChipLayout = chipLayout;
  } else {
      m_FreeChipLayout = true;
      m_ChipLayout = NULL;
  }
  m_Report = NULL;
}

GenoQC::~GenoQC() {
  clear();
}

void GenoQC::clear()
{
  for(int cl = 0; cl < (int)m_CelListeners.size(); cl ++)
    Freez(m_CelListeners[cl]);
  m_CelListeners.clear();
  Freez(m_Report);
  //
  if (m_kitao_db!=NULL) {
    delete m_kitao_db;
    m_kitao_db=NULL;
  }
}

void GenoQC::defineOptions() {

  defineOptionSection("Input Options");

  defineOption("c", "cdf-file", PgOpt::STRING_OPT,
                     "File defining probe sets.",
                     "");
  defineOption("", "spf-file", PgOpt::STRING_OPT,
                     "SPF File defining probe sets.",
                     "");
  defineOption("q", "qcc-file", PgOpt::STRING_OPT,
                     "File defining QC probesets.",
                     "");
  defineOption("a", "qca-file", PgOpt::STRING_OPT,
                     "File defining QC analysis methods.",
                     "");
  defineOption("", "cel-files", PgOpt::STRING_OPT,
                     "Text file specifying cel files to process, one per line with the first line being 'cel_files'.",
                     "");
  defineOption("", "chrX-probes", PgOpt::STRING_OPT,
                     "File containing probe_id (1-based) of probes on chrX. "
                     "Used for copy number probe chrX/Y ratio gender calling. "
                     "[Experimental]",
                     "");
  defineOption("", "chrY-probes", PgOpt::STRING_OPT,
                     "File containing probe_id (1-based) of probes on chrY. "
                     "Used for copy number probe chrX/Y ratio gender calling. "
                     "[Experimental]",
                     "");
  defineOption("", "chrZ-probes", PgOpt::STRING_OPT,
                     "File containing probe_id (1-based) of probes on chrZ. "
                     "Used for copy number probe chrZ/W ratio avian gender calling. "
                     "[Experimental]",
                     "");
  defineOption("", "chrW-probes", PgOpt::STRING_OPT,
                     "File containing probe_id (1-based) of probes on chrW. "
                     "Used for copy number probe chrZ/W ratio avian gender calling. "
                     "[Experimental]",
                     "");
  defineOption("", "probe-class-file", PgOpt::STRING_OPT,
                    "File containing probe_id (1-based) of "
                          "probes and a 'class' designation. Used to  "
                          "compute mean probe intensity by class for "
                          "report file.",
                    "");
  defineOption("","target-sketch", PgOpt::STRING_OPT,
                         "Sketch file.",
                         "");
  defineOption("","channel-file", PgOpt::STRING_OPT,
               "Channel file.",
               "");
  //
  defineOption("","reagent-kit-discriminator",PgOpt::STRING_OPT,
               "list of probeset names with pc1s and means to use for "
               "classifying the reagent kits.",
               "");
//  defineOption("","reagent-kit-discriminator",PgOpt::STRING_OPT,
//               "list of probeset names with pc1s and means to use for "
//               "classifying the reagent kits.",
//               "");

  defineOptionSection("Output Options");

  defineOption("", "out-file", PgOpt::STRING_OPT,
               "Name to use for the output file.",
               "apt-geno-qc.report.txt");
  defineOption("", "dm-out", PgOpt::STRING_OPT,
               "Folder to use for DM output. Enables DM "
               "output. One per CEL file. [experimental]",
               "");

  defineOptionSection("Analysis Options");

  defineOption("", "dm-het-mult", PgOpt::DOUBLE_OPT,
                     "DM Het Mult to use for DM output. [default '1.25']",
                     "1.25");
  defineOption("", "dm-thresh", PgOpt::DOUBLE_OPT,
                     "DM threshold to use for making no calls. [default '0.33']",
                     "0.33");
  defineOption("", "female-thresh", PgOpt::DOUBLE_OPT,
                    "Threshold for calling females when using cn-probe-chrXY-ratio or cn-probe-chrZW-ratio method.",
                    "0.48");
  defineOption("", "male-thresh", PgOpt::DOUBLE_OPT,
                    "Threshold for calling females when using cn-probe-chrXY-ratio or cn-probe-chrZW-ratio method.",
                    "0.71");

  defineOptionSection("Engine Options (Not used on command line)");

  defOptMult("", "cels", PgOpt::STRING_OPT,
             "Cel files to process.",
             "");
}

void GenoQC::defineStates() { }

  /**
   * @brief Populate the input qcAnalysisOptions from qca file
   * @param qcaFile - the file to read
   */
void GenoQC::readQCAFile(std::string qcaFile) {
  APT_ERR_ASSERT(qcaFile!="","qcaFile cannot be blank.");
  QCAnalysisOptions::readQCAFile(qcaFile, qcAnalysisOpts);
}

  /**
   * @brief Populate the input qcProbesetOptions from qcc file
   * @param qccFile - the file to read
   */
void GenoQC::readQCCFile(std::string qccFile) {
  APT_ERR_ASSERT(qccFile!="","qccFile cannot be blank.");
  QCProbesetOptions::readQCCFile(qccFile, qcProbesetOpts);
  // copy the temp-dir to the options as it will be passed down into GimmieContrast
  qcProbesetOpts.m_tempdir=getOpt("temp-dir");
  //printf("### JHG: readQCCFile(): temp-dir='%s'\n",qcProbesetOpts.m_tempdir.c_str());
}

/**
 * Make sure that our options are sane. Call Err::errAbort if not.
 */
void GenoQC::checkOptionsImp()
{
    defineStates();

    setLibFileOpt("cdf-file");
    setLibFileOpt("spf-file");
    setLibFileOpt("qcc-file");
    setLibFileOpt("qca-file");
    setLibFileOpt("target-sketch");
    setLibFileOpt("channel-file");
    setLibFileOpt("chrX-probes");
    setLibFileOpt("chrY-probes");
    setLibFileOpt("chrZ-probes");
    setLibFileOpt("chrW-probes");
    setLibFileOpt("probe-class-file");
    setLibFileOpt("reagent-kit-discriminator");

    if(getOpt("out-file") == "") {
        Err::errAbort("Must provide an output file using --out-file option.");
    }
    
    /* Read in cel file list from other file if specified. */
    vector<string> celFiles;
    EngineUtil::getCelFiles(celFiles, this);
    if (celFiles.size() == 0) {
      Err::errAbort("No cel files specified.");
    }
    setOpt("cels",celFiles);

    //std::string tmpstr;
    //tmpstr=getOpt("temp-dir");
    //printf("### JHG: GenoQC::checkOptionsImp(): temp-dir-1: '%s'\n",tmpstr.c_str());

    /// @todo this should be in BaseEngine?
    if (getOpt("temp-dir")=="") {
      std::string temp_dir=Fs::dirname(getOpt("out-file"));
      if (temp_dir=="") {
        temp_dir=".";
      }
      setOpt("temp-dir",temp_dir);
    }

    //tmpstr=getOpt("temp-dir");
    //printf("### JHG: GenoQC::checkOptionsImp(): temp-dir-2: '%s'\n",tmpstr.c_str());

    readQCAFile(getOpt("qca-file"));
    readQCCFile(getOpt("qcc-file"));

  // no setLibFileOpt call as we dont want it them be expanded.
  qcProbesetOpts.m_sketchInFile=getOpt("target-sketch");
  qcProbesetOpts.m_channelFile=getOpt("channel-file");

  ///@todo allow computing of DM calls when no qcc/qca file
    int nMethods = qcAnalysisOpts.methods.size();
    if((nMethods < 0) && getOpt("dm-out") == "")
        Err::errAbort("There is at least one method need to be specified");

    // check that SNP lists are available
    for(int i = 0; i < nMethods; i ++) {
        string key = qcAnalysisOpts.methods[i].groupName;

        // snp list specified by a map
        if(key != "") 
            if(qcProbesetOpts.probesetGroups.find(key) == 
                qcProbesetOpts.probesetGroups.end())
                Err::errAbort("SNP list:" + qcAnalysisOpts.methods[i].groupName + " does not exist");

    }

    // This is not great in that it causes a copy of the 
    // probeset lists - which can be quite large. On the
    // other hand keeping just a pointer seems a bit
    // dangerous as one could initialize with one set of options
    // and then modify them before computing the stats. 
    if(getOptBool("force"))
        Verbose::out(1, "WARNING: Chip type force option was chosen, no chip type validation will be done!");
    else
        CheckChipType();

    // create output folder for DM calls
	std::string strDmOut = getOpt("dm-out");
    if(strDmOut != "") {
        if(!Fs::isWriteableDir(strDmOut)) {
            if(Fs::mkdirPath(strDmOut, false) != APT_OK) {
                Err::errAbort("Can't make or write to directory: " + strDmOut);
            }
        }
    }

}

// Compute the various metrics over all chips.
void GenoQC::runImp() {
  //
  if(m_ChipLayout == NULL)
    m_ChipLayout = new ChipLayout();
  //
  map<string,bool> psMap;
  map<string,bool>::iterator psMapIter;
  std::set<const char *, Util::ltstr> probeSetsToLoad;

 try { // outer try for handling exception
  try { // inner try for memory cleanup
    Verbose::out(1, "Initializing QC Process.");
    Verbose::out(1, "MAJOR PROGRESS UPDATE: Running GenoQC.");

    int nChips = (int)getOptVector("cels").size();
    int nMethods = qcAnalysisOpts.methods.size();
    
    // KitAO - part 1
    if (getOpt("reagent-kit-discriminator")!="") {
      m_kitao_db=new KitAODb();
      APT_ERR_ASSERT(m_kitao_db!=NULL,"internal error.");
      m_kitao_db->readTsv(getOpt("reagent-kit-discriminator"));
      // add these names to the probesets to load in case we end up loading
      // m_ChipLayout on our own.
      for (KitAODb::name2entry_iter_t entry_i=m_kitao_db->m_name2entry.begin();
           entry_i!=m_kitao_db->m_name2entry.end();
           entry_i++) {
        probeSetsToLoad.insert(entry_i->first.c_str());
      }
    }


    // Figure out the union of all qc probeset lists
    // if the user wants dm calls, then we leave the set empty so that
    // everything is loaded.
    if (getOpt("dm-out") == "") {
      Verbose::out(1,"Figuring out what probesets to load...");
      //
      std::string ps_name;
      std::string group_name;
      vector<pair<string,string> > *vPtr;

      for(int m = 0; m < nMethods; m ++) {
        group_name=qcAnalysisOpts.methods[m].groupName;
        Verbose::out(1,"Loading group '"+group_name+"'");
        //vector<string> *vPtr = &qcProbesetOpts.probesetGroups[qcAnalysisOpts.methods[m].groupName];
        //
        vPtr=&qcProbesetOpts.probesetGroups[group_name];
        for (int i=0; i < vPtr->size(); i++) {
          psMap[vPtr->at(i).first] = true;
        }
        // multiple probeset lists are required by this method 
        // and we need to get the lists from the options
        if (qcAnalysisOpts.methods[m].analysis == QCMethod::ADJHOMHILO) {
          vPtr = &qcProbesetOpts.probesetGroups[qcAnalysisOpts.methods[m].options["ps.rand"]];
          for(int i=0; i < vPtr->size(); i++) {
            psMap[vPtr->at(i).first] = true;
          }
          vPtr = &qcProbesetOpts.probesetGroups[qcAnalysisOpts.methods[m].options["ps.set1"]];
          for(int i=0; i < vPtr->size(); i++) {
            psMap[vPtr->at(i).first] = true;
          }
          vPtr = &qcProbesetOpts.probesetGroups[qcAnalysisOpts.methods[m].options["ps.set2"]];
          for (int i=0; i < vPtr->size(); i++) {
            psMap[vPtr->at(i).first] = true;
          }
        }
      }
      //
      for(psMapIter=psMap.begin(); psMapIter != psMap.end(); psMapIter++) {
        probeSetsToLoad.insert(psMapIter->first.c_str());
        //printf("probeSetsToLoad: '%s'\n",psMapIter->first.c_str());
      }
    }

    // get the chip layout if one was not passed in to the constructor
    if (m_FreeChipLayout) {
        std::set<affxcdf::GeneChipProbeSetType> psTypesToLoadEmpty;
        std::vector<bool> probeSubsetEmpty;
        std::string stringEmpty;
        std::string strCdfFile = getOpt("cdf-file");
        std::string strSpfFile = getOpt("spf-file");
        if (strCdfFile != "") {
          Verbose::out(1,"Opening cdf file: " + strCdfFile);
          m_ChipLayout->openCdf(strCdfFile,probeSetsToLoad, probeSubsetEmpty, stringEmpty);
        } 
        else if(strSpfFile != "") {
          Verbose::out(1,"Opening spf file: " + strSpfFile);
          m_ChipLayout->openSpf(strSpfFile,
                                probeSetsToLoad,
                                NULL,
                                probeSubsetEmpty,
                                stringEmpty,
                                false,
                                psTypesToLoadEmpty
                                );
        } 
        else {
          Err::errAbort("Must provide a CDF or SPF file.");
        }
    }
    
    // create our reporter
    std::vector<std::string> cels = getOptVector("cels");
    m_Report = new QuantMethodRunReport(cels);
    m_Report->setFormat(affx::TsvReport::FMT_TSV);
    m_Report->setUseDefaultSuffix(0); // avoid default report.txt suffix
    m_Report->setFilename(getOpt("out-file"));

    // Go through each method and setup a CEL listener, attach to reporter
    for(int m = 0; m < nMethods; m ++) {
        m_CelListeners.push_back(
            qcAnalysisOpts.methods[m].createCelListener(
                m_ChipLayout,cels,qcProbesetOpts));
        ChipSummary *chipSummary = dynamic_cast<ChipSummary *>(m_CelListeners[m]);
        if(chipSummary != NULL)
            m_Report->registerChipSummary(chipSummary);

        // save some memory
        DmListener *dmListener = dynamic_cast<DmListener *>(m_CelListeners[m]);
        if(dmListener != NULL)
            dmListener->setKeepGenotypes(false);
    }

    // Create other CelListeners not specified in QCA file
    // CN probe based gender calling
    bool zw_gender_calling = false;
    if (getOpt("chrZ-probes") != "") {
        setOpt("chrX-probes", getOpt("chrZ-probes"));
        setOpt("chrY-probes", getOpt("chrW-probes"));
        zw_gender_calling = true;
    }
    if(getOpt("chrX-probes") != "") {
        CnProbeGenderCelListener *cnProbeCelListener = 
            new CnProbeGenderCelListener(getOpt("chrX-probes"), getOpt("chrY-probes"), getOptDouble("female-thresh"), getOptDouble("male-thresh"), zw_gender_calling);
        m_CelListeners.push_back(cnProbeCelListener);

        ChipSummary *chipSummary = dynamic_cast<ChipSummary *>(cnProbeCelListener);
        if(chipSummary != NULL)
            m_Report->registerChipSummary(chipSummary);
    }

    // Cel Stat Listener
    if(getOpt("probe-class-file") != "") {
        map<string, vector<bool> > masks;
    
        vector<bool> pmMask = m_ChipLayout->getPmProbeMask();
        masks["pm"] = pmMask;
    
        vector<bool> mmMask = m_ChipLayout->getMmProbeMask();
        int mmCount = 0;
        for(int i=0; i<mmMask.size(); i++)
            if(mmMask[i])
                mmCount++;
        if(mmCount > 5000) {
            masks["mm"] = mmMask;
        }

        map<string, vector<bool> >::iterator maskIter;
        EngineUtil::readProbeClassFile(getOpt("probe-class-file"), m_ChipLayout->getProbeCount(), masks);
        for(maskIter = masks.begin(); maskIter != masks.end(); maskIter++) {
            int count = 0;
            for(vector<bool>::iterator mIx = maskIter->second.begin(); mIx != maskIter->second.end(); mIx++)
                if(*mIx) 
                    count++;
            Verbose::out(2,"Number of " + maskIter->first + " probes is: " + ToStr(count));
        }

        CelStatListener *celStats = new CelStatListener(masks);
        m_CelListeners.push_back(celStats);

        ChipSummary *chipSummary = dynamic_cast<ChipSummary *>(celStats);
        if(chipSummary != NULL)
            m_Report->registerChipSummary(chipSummary);
    }
    /// @todo  DM chrX het rate gender calling?

    // setup dmProbesets
    vector<ProbeListPacked> dmProbesets;
    if(getOpt("dm-out") != "") {
        for(uint32_t i = 0; i < m_ChipLayout->getProbeSetCount(); i++) {
            ProbeListPacked pl = m_ChipLayout->getProbeListAtIndex(i);
            ProbeSet *ps = ProbeListFactory::asProbeSet(pl);
            if(DmListener::okForDm(ps))
                dmProbesets.push_back(pl);
            delete ps;
        }
    }

    // KitAO - part 2
    if (getOpt("reagent-kit-discriminator")!="") {
      APT_ERR_ASSERT(m_kitao_db!=NULL,"internal error.");
      //
      KitAOCelListener* kitaocellistener=new KitAOCelListener();
      kitaocellistener->setDb(m_kitao_db);
      kitaocellistener->setLayout(m_ChipLayout);

      kitaocellistener->m_sketchFileName=getOpt("target-sketch");
      APT_ERR_ASSERT(kitaocellistener->m_sketchFileName!="",
                     "a target sketch is required.");
      //
      m_CelListeners.push_back(kitaocellistener);
      ChipSummary *chipSummary = dynamic_cast<ChipSummary *>(kitaocellistener);
      m_Report->registerChipSummary(chipSummary);
    }

    //// now the work is done over each cel file.
    Verbose::out(1, "Processing CEL files for genotype QC analysis:");
    for(int celIx = 0; celIx < nChips; celIx ++)
    {
        affymetrix_fusion_io::FusionCELData CEL;
        CEL.SetFileName(cels[celIx].c_str());

        if(!CEL.Exists()) 
            Err::errAbort("CEL file " + cels[celIx] + " does not exist");
    
        Verbose::out(1, "Processing " + ToStr(celIx + 1) + " of " + ToStr(nChips) + ": " + CEL.GetFileName());
    
        // read in the cel file
        if(!CEL.Read())
            Err::errAbort("Can not read CEL file " + cels[celIx]);

        // Pass to each CEL listener
        for(int cl = 0; cl < m_CelListeners.size(); cl ++) 
          m_CelListeners[cl]->newChip(&CEL);
        
        // Generate DM call output if requested
        if(getOpt("dm-out") != "") {
            affx::TsvFile tsv;
            tsv.defineColumn(0,0,"probeset_id");
            tsv.defineColumn(0,1,"call");
            tsv.defineColumn(0,2,"confidence");
            tsv.addHeaderComment("Calls: -1=NN, 0=AA, 1=AB, 2=BB");
            tsv.addHeader("cdf-file",getOpt("cdf-file"));
            tsv.addHeader("spf-file",getOpt("spf-file"));
            tsv.addHeader("cel-file",cels[celIx]);
            tsv.addHeader("number-SNPs",(int)dmProbesets.size());
            tsv.addHeader("dm-thresh",getOpt("dm-thresh"));
            tsv.addHeader("dm-het-mult",getOpt("dm-het-mult"));
            string filename = cels[celIx];
            filename = Fs::noextname1(filename)+".dm.txt";
            filename = Fs::join(getOpt("dm-out"),Fs::basename(filename));
            tsv.writeTsv_v1(filename);

            for(int psIx = 0; psIx < dmProbesets.size(); psIx++){
                ProbeSet *ps = ProbeListFactory::asProbeSet(dmProbesets[psIx]);
                vector<CQuartet> qVec;
                DmListener::fillInQuartet(qVec,ps,&CEL);
                std::pair<float,int> call;
                string name = ps->name;
                tsv.set(0,0,name);
                tsv.set(0,1,(int)-1);
                tsv.set(0,2,(float)1.0);
                if(CallDM(qVec,call,getOptDouble("dm-het-mult"))) {
                    if(call.first < getOptDouble("dm-thresh"))
                        tsv.set(0,1,call.second);
                    else
                        tsv.set(0,1,(int)-1);
                    tsv.set(0,2,call.first);
                }else{
                    Verbose::out(1,"DM Call on probeset " + name + " failed.");
                }
                tsv.writeLevel(0);
                delete ps;
            }
        }

        CEL.Close();
    } 

    // Flush metrics to report file
    if(!m_Report->finish()) {
      Err::errAbort("Problem writing the report file " + getOpt("out-file"));
    }
    //
    probeSetsToLoad.clear();
    psMap.clear();

  } // end inner try
  catch (...) {
    // free up memory
   
    if(m_FreeChipLayout)
        Freez(m_ChipLayout);
    clear();

    //re-throw
    throw;
  }
 } // end outer try
 /* When things go wrong see if we can die gracefully here. */
 catch(Except &e) {
   Verbose::out(0,"");
   Err::errAbort(e.what());
 }
 catch(const std::bad_alloc &e) {
   Verbose::out(0,"");
   Err::errAbort("Ran out of memory. "
                 "Try using quitting other applications.");
 }
 catch(affymetrix_calvin_exceptions::CalvinException &ce) {
   Verbose::out(0,"");
   /// @todo The error message is null for most calvin files.
   ///       But we should mark this as a calvin exception so users know.
   Err::errAbort("Affymetrix GeneChip Command Console library has thrown an exception. "
                 "Bad Calvin File: "
                 "Description: '" +StringUtils::ConvertWCSToMBS(ce.Description()) + "'");
 }
 catch(BroadException &e) {
   Verbose::out(0,"");
   Err::errAbort("Exception. msg: '" + ToStr(e.m_msg) + "' source file: '" + ToStr(e.m_sourcefile) +
                 ":" + ToStr(e.m_sourceline) + "'");
 }
 catch(const std::exception &e) {
   Verbose::out(0,"");
   Err::errAbort("Exception caught. "
                 "Message is: " + ToStr(e.what()));
 }
 catch(...) {
   Verbose::out(0,"");
   Err::errAbort("Unknown exception caught.");
 }

 // free up memory
 if(m_FreeChipLayout)
    Freez(m_ChipLayout);
 clear();
}

// Check the chip type
void GenoQC::CheckChipType() {
    Verbose::out(1, "Checking chip type on all cel files.");
    colrow_t numRows, numCols;
    int probeCount;
    int probeSetCount;
    vector<string> chipTypes;

    std::string strCdfFile = getOpt("cdf-file");
    std::string strSpfFile = getOpt("spf-file");
    std::vector<std::string> cels = getOptVector("cels");

    if(strCdfFile != "") {
      EngineUtil::getCdfChipType(chipTypes, numRows, numCols, probeCount, probeSetCount, strCdfFile);
    }
    else {
      EngineUtil::getSpfChipType(chipTypes, numRows, numCols, probeCount, probeSetCount, strSpfFile);
    }

    EngineUtil::checkCelChipTypes(chipTypes, probeCount, cels, numRows, numCols);
    if ((getOpt("chrX-probes") != "") && !getOptBool("force")) {
      EngineUtil::checkTsvFileChipType(getOpt("chrX-probes"), chipTypes);
    }
    if ((getOpt("chrY-probes") != "") && !getOptBool("force")) {
      EngineUtil::checkTsvFileChipType(getOpt("chrY-probes"), chipTypes);
    }
}

