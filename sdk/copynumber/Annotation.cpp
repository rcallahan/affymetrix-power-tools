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
 * @file   Annotation.cpp
 * @brief  Wrapper for the SQLite Annotation database.
 *
 */

#include "copynumber/Annotation.h"
#include "copynumber/CNAnalysisMethodCovariateParams.h"
#include "copynumber/CNAnalysisMethodCovariateSignalAdjuster.h"
//
#include "copynumber/CNAnalysisMethod.h"
#include "copynumber/CNLog2RatioData.h"
//
#include "file/TsvFile/TsvFile.h"
#include "file5/File5.h"
#include "file5/File5_File.h"
#include "file5/File5_Tsv.h"

#include "util/Util.h"
#include "util/Verbose.h"
//

std::vector<affymetrix_calvin_parameter::ParameterNameValueType> Annotation::m_vParams;

// A test function.
void Annotation::test()
{
    AffxString strFileNameIn = "c:/annotation/apt-cyto-annot-chipc.txt";
    std::string strFileName = "c:/annotation/cyto2_HG49m_C.annot.a5";
    affx::File5_File file5;
    file5.open(strFileName, affx::FILE5_CREATE | affx::FILE5_REPLACE);
    affx::File5_Group* group5 = file5.openGroup("Annotation", affx::FILE5_REPLACE);
    affx::File5_Tsv* tsv5 = group5->openTsv("Annotation", affx::FILE5_REPLACE);
    tsv5->defineColumn(0, 0, "probeset_id", affx::FILE5_DTYPE_STRING, 40);
    tsv5->defineColumn(0, 1, "Chromosome", affx::FILE5_DTYPE_INT);
    tsv5->defineColumn(0, 2, "Position", affx::FILE5_DTYPE_INT);
    tsv5->defineColumn(0, 3, "GC", affx::FILE5_DTYPE_DOUBLE);
    affx::TsvFile tsv;
    tsv.m_optAutoTrim = true;
    AffxByteArray baLine;
    AffxByteArray baColumn;

    AffxString strCol1, strCol2;
    int iCol3 = 0;
    double dCol4 = 0.0;
    
    tsv.open(strFileNameIn);
    tsv.bind(0,0,&strCol2,  affx::TSV_BIND_REQUIRED);
    tsv.bind(0,1,&strCol2,  affx::TSV_BIND_REQUIRED);
    tsv.bind(0,2,&iCol3,  affx::TSV_BIND_REQUIRED);
    tsv.bind(0,3,&dCol4,  affx::TSV_BIND_REQUIRED);
    
    while (tsv.nextLevel(0) == affx::TSV_OK) {
      tsv5->set_string(0, 0, strCol1);
      strCol2 = strCol2.substring(3);
      if (strCol2 == "X") {tsv5->set_i(0, 1, 24);}
      else if (strCol2 == "Y") {tsv5->set_i(0, 1, 25);}
      else {tsv5->set_i(0, 1, ::getInt(strCol2));}
      tsv5->set_i(0, 2, iCol3);
      tsv5->set_d(0, 3, dCol4);
      tsv5->writeLevel(0);
    }
    tsv.clear();
    tsv5->close();
    delete tsv5;
    group5->close();
    delete group5;
    file5.close();
}

/**
 * @brief Constructor
 */
Annotation::Annotation() : SQLiteDatabase()
{
}

/**
 * @brief Destructor
 */
Annotation::~Annotation()
{
}

/**
 * @brief Setup a display string from the run time.
 * @param const AffxString& - The name of the action we are timing.
 * @param time_t - The run time in seconds.
 * @return AffxString - The resulting string
 */
AffxString Annotation::runTime(const AffxString& strAction, time_t runTime)
{
    AffxString str = "";
    double dRunTime = (double)runTime;

    if (dRunTime < 60)
    {
        str = strAction + "  Run time = ";
        str += ::getDouble(dRunTime, 3, true);
        str += " seconds";
    }
    else if (dRunTime < (60 * 60))
    {
        str = strAction + "  Run time = ";
        str += ::getDouble((dRunTime) / 60, 3, true);
        str += " minutes";
    }
    else
    {
        str = strAction + "  Run time = ";
        str += ::getDouble((dRunTime) / (60 * 60), 3, true);
        str += " hours";
    }
    return str;
}

void Annotation::loadAnnotationForCopynumber(   BaseEngine& engine,
                                                bool bNewProbeSets,
                                                AffxSplitArray<CNProbeSet> &vProbeSets,
                                                AffxArray<CNProbeSet>& arProbeSets)
{
        time_t startTime = time(NULL);
    if ((engine.isOptDefined("annotation-file")) && (engine.getOpt("annotation-file") != ""))
    {
        if (bNewProbeSets) {Annotation().newCNProbeSetsSQLite(engine, vProbeSets);}
        else {Annotation().loadAnnotationSQLite(engine, arProbeSets);}
    }
    else
    {
        Err::errAbort("annotation-file parameter has not been specified.");
    }
        time_t endTime = time(NULL);
    Verbose::out(1, runTime("Load NetAffx", endTime - startTime));
}

bool Annotation::isCytoScanHD(const std::string& strFileName)
{
	bool bCytoScanHD = false;   
    AffxString strSQL = "select probeSet_ID, chr_id, start, Percent_GC, chrx_par, Chromosome.shortname, enzyme_fragment, process_flag, probe_count, probeset_type, probe_gc_a from Annotations, Chromosome where Annotations.chr_id = Chromosome.key and process_flag > 0";
    SQLiteDatabase db;
	try 
	{
		db.open(strFileName);
		SQLiteRecordset rset(db);
		try {rset.open(strSQL); rset.close(); bCytoScanHD = true;} catch(...) {}
		db.close();
    }
    catch (...) {db.close();}
	return bCytoScanHD;
}

void Annotation::newCNProbeSetsSQLite(BaseEngine& engine, AffxSplitArray<CNProbeSet> &vProbeSets)
{
    Verbose::out(3, "Annotation::newCNProbeSets");
    int iXChromosome = engine.getOptInt("xChromosome");
    int iYChromosome = engine.getOptInt("yChromosome");
    AffxString strFileName = engine.getOpt("annotation-file");
    if (affx::File5_File::isHdf5file(strFileName)) {return newCNProbeSetsHDF5(engine, vProbeSets);}
    bool bEnzymes = !((engine.isOptDefined("cyto2")) && (engine.getOptBool("cyto2")));
    AffxArray<AffxString> vStyAdapters;
    AffxArray<AffxString> vNspAdapters;
    vProbeSets.clear();
    m_vParams.clear();

    AffxArray<AffxString> arRestrictList;
    bool bRestrictList = (engine.getOpt("probeset-ids") != "");
    if (bRestrictList)
    {
        loadProbeSetNamesFromRestrictList(engine, arRestrictList);
    }

    AffxByteArray ba;
    AffxString str;
    AffxString strProbeSetName;
    bool bSnpInterference = false;
    int iSearchIndex = 0;
    affymetrix_calvin_parameter::ParameterNameValueType param;
    unsigned int uiProbeSetIndex = 0;
    AffxString strSQL1;
    AffxString strSQL2;
    AffxString strSQL3;
    unsigned int uiRowCount = 0;
    try {
    open(strFileName);
    SQLiteRecordset rset(*this);

    strSQL2 = "select key, value from information";
    rset.open(strSQL2);
    int iKeyCount = 0;
    while (rset.fetch())
    {
        ba.assign(rset.getString(0));
        ba.replace(" ", "-");
        ba.replace("_", "-");
        ba.replace("netaffx-annotation-netaffx-build", "netaffx-build");
        if (ba.equals("genome-version-ucsc")) {iKeyCount++;}
        if (ba.equals("dbsnp-version")) {iKeyCount++;}
        if (ba.equals("netaffx-build")) {iKeyCount++;}
        param.SetName(StringUtils::ConvertMBSToWCS("affymetrix-algorithm-param-" + ba.toString()));
        param.SetValueAscii(rset.getString(1));
        m_vParams.push_back(param);
    }
    rset.close();
    if (iKeyCount < 3) {throw(Except("The annotation file must contain all three of these parameters: genome-version-ucsc, dbsnp-version, and netaffx-build"));}

    bool bCyto2 = false;	
	bool bCytoScanHD = false;
    strSQL1 = "select count(*) from Annotations, Chromosome where Annotations.chr_id = Chromosome.key and process_flag > 0";
    strSQL2 = "select count(*) from Annotations, Chromosome where Annotations.chr_id = Chromosome.key and in_final_list = 1";
    try {
        rset.open(strSQL1);
        bCyto2 = true;
    } catch (...) {
        rset.open(strSQL2);
        bCyto2 = false;
    }
    if (rset.fetch())
    {
        uiRowCount = rset.getInteger(0);
    }
    rset.close();

    strSQL1 = "select probeSet_ID, chr_id, start, Percent_GC, chrx_par, Chromosome.shortname, process_flag, probe_count, ProbeSet_Type  from Annotations, Chromosome where Annotations.chr_id = Chromosome.key and process_flag > 0";
    strSQL2 = "select probeSet_ID, chr_id, start, Percent_GC, chrx_par, Chromosome.shortname, enzyme_fragment, snp_interference from Annotations, Chromosome where Annotations.chr_id = Chromosome.key and in_final_list = 1";
    strSQL3 = "select probeSet_ID, chr_id, start, Percent_GC, chrx_par, Chromosome.shortname, enzyme_fragment, process_flag, probe_count, probeset_type, probe_gc_a, fragment_gc from Annotations, Chromosome where Annotations.chr_id = Chromosome.key and process_flag > 0";
    if (uiRowCount > 1)
    {
        Verbose::progressBegin(1, "Annotation::newCNProbeSetsSQLite(...) ", (uiRowCount / 100000), 1, (uiRowCount / 100000));
        vProbeSets.allocate(uiRowCount);

		try {
            rset.open(strSQL3);
            bCyto2 = false;
            bCytoScanHD = true;
        }
		catch(...) 
		{
			try {
                rset.open(strSQL1);
                bCyto2 = true;
            }
			catch(...) 
			{
				rset.open(strSQL2);
                bCyto2 = false;
			}
		}
        
        CovariateParams::determineCovariateMap(engine);
        int fragmentAdapterTypeIndex = CovariateParams::mapCovariateNameToIndex("fragment-adapter-type");
        int fragmentLengthIndex = CovariateParams::mapCovariateNameToIndex("fragment-length");
        int fragmentGCIndex = CovariateParams::mapCovariateNameToIndex("fragment-gc");
        int probeGCIndex = CovariateParams::mapCovariateNameToIndex("probe-gc");
        int localGCIndex = CovariateParams::mapCovariateNameToIndex("local-gc");
        int markerClassIndex = CovariateParams::mapCovariateNameToIndex("marker-class");


        AffxString strChromosome;
        char cChromosome = 0;
        int i = 0;
        bool bInWorkflowEngine = (!(engine.isOptDefined("cyto2") && engine.getOptBool("cyto2")) && !engine.getOptBool("cytoscan-hd"));
        while (rset.fetch())
        {
            if ((i % 100000) == 0) {Verbose::progressStep(1);}
            i++;
            strChromosome = rset.getString(5);
            if (strChromosome.length() > 0) {cChromosome = strChromosome.getAt(0);}
            strProbeSetName = rset.getString(0);
            strProbeSetName.removeSurroundingQuotes();
            CNProbeSet* p = vProbeSets.getAt(uiProbeSetIndex);
            if (arRestrictList.getCount() > 0)
            {
                iSearchIndex = arRestrictList.binarySearch(strProbeSetName, 0);
                if (iSearchIndex == -1) {p->setProcess(false);}
                else {p->setProcess(true);}
            }

            p->setProbeSetName(strProbeSetName);
            std::string strProbeSetType;
			if (bCyto2)
			{
				p->setProcessFlag(rset.getInteger(6));
				strProbeSetType = rset.getString(8);
				if(strProbeSetType == "SNP")
				{
					p->setTrulySNP(true);
				}
				if(strProbeSetType == "CN")
				{
					p->setTrulyCN(true);
				}
			}
			else
			{
				if (bCytoScanHD)
				{
					p->setProcessFlag(rset.getInteger(7));
					strProbeSetType = rset.getString(9);
					if(strProbeSetType == "SNP")
					{
						p->setTrulySNP(true);
					}
					if(strProbeSetType == "CN")
					{
						p->setTrulyCN(true);
					}
				}
				else
				{
					if (strProbeSetName.startsWith("SNP"))
					{
						p->setProcessFlag(3);
						p->setTrulySNP(true);
					}
					else
					{
						p->setProcessFlag(1);
						p->setTrulyCN(true);
					}
				}
			}

            p->setChromosome(rset.getInteger(1));
            p->setPosition(rset.getInteger(2));
            double gcContent = rset.getDouble(3);
            if (bInWorkflowEngine && (gcContent != gcContent)) { p->setGCContent(0.0f); }
            else  { p->setGCContent((float)gcContent); }
            p->setPseudoAutosomalRegion(rset.getInteger(4) > 0);
			if (bCytoScanHD) {p->setReplicateCount(rset.getInteger(8));}
			else {p->setReplicateCount(rset.getInteger(7));}
            float probe_gc_a = 0.0f;
            if (!bCyto2)
            {
                AffxString strFragment = rset.getString(6);
                strFragment.removeSurroundingQuotes();
                char nspIcode = parseAdapterCode(strFragment, "NspI", vNspAdapters);
                p->setStyAdapterCode(parseAdapterCode(strFragment, "StyI", vStyAdapters));
                p->setNspAdapterCode(nspIcode);
                if (!bCytoScanHD)
                {
                    bSnpInterference = (rset.getInteger(7) == 1);
                }
                else
                {
                    bSnpInterference = false;

                    probe_gc_a = (float)rset.getDouble(10);
                    double fragmentGC = rset.getDouble(11);

                    // Allocate space and assign initial NaN value for all covariates.
                    if (CovariateParams::m_allCovariateMap.size() > 0)
                    {
                        p->setCovariateValue(CovariateParams::m_allCovariateMap.size()-1, std::numeric_limits<float>::quiet_NaN());
                    }
                    // Add covariate information.
                    if (fragmentAdapterTypeIndex > -1)
                    {
                        p->setCovariateValue(fragmentAdapterTypeIndex, 
                            ((nspIcode == -1) ? std::numeric_limits<float>::quiet_NaN() : nspIcode));
                    }
                    if (fragmentLengthIndex > -1) 
                    {
                        int fragmentLength = parseFragmentLength(strFragment);
                        p->setCovariateValue(fragmentLengthIndex, 
                            ((fragmentLength == -1) ? std::numeric_limits<float>::quiet_NaN() : fragmentLength));
                    }
                    if (probeGCIndex > -1) p->setCovariateValue(probeGCIndex, probe_gc_a);
                    if (localGCIndex > -1) p->setCovariateValue(localGCIndex, gcContent);
                    if (fragmentGCIndex > -1) p->setCovariateValue(fragmentGCIndex, fragmentGC);
                    if (markerClassIndex > -1) p->setCovariateValue(markerClassIndex, (p->processAsSNPNormalize() ? 3.f :
                                                                                          (p->processAsCNNormalize() ? 1.f :
                                                                                              std::numeric_limits<float>::quiet_NaN())));
                }
            }
            else
            {
                p->setStyAdapterCode(-1);
                p->setNspAdapterCode(-1);
                bSnpInterference = false;
            }
            bool bBadAnnotation = false;
            if (cChromosome == 'X') {p->setChromosome(iXChromosome);}
            else if (cChromosome == 'Y') {p->setChromosome(iYChromosome);}
            if (!bCyto2)
            {
                if ((cChromosome == 'M') || (cChromosome == 'm')) {p->setChromosome(0);}
                if (p->getPosition() == 0) {bBadAnnotation = true;}
                if (p->getGCContent() == 0) {bBadAnnotation = true;}                    // SNP6
                if (p->getGCContent() != p->getGCContent()) {bBadAnnotation = true;}    // CytoScanHD
                if ((bEnzymes) && (!p->isSty()) && (!p->isNsp())) {bBadAnnotation = true;}
                if (bSnpInterference) {bBadAnnotation = true;}
                if (probe_gc_a != probe_gc_a) {bBadAnnotation = true;}
                if (bCytoScanHD && strProbeSetType == "") {bBadAnnotation = true;}
            }
            if (p->getChromosome() == 0) {bBadAnnotation = true;}
            if (!bBadAnnotation)
            {
                uiProbeSetIndex++;
            }
            else
            {
                p->clear();
            }
        }
        rset.close();
        Verbose::progressEnd(1, "Done.");
    }

    close();
    }
    catch (Except&) {
        vStyAdapters.deleteAll();
        vNspAdapters.deleteAll();
        arRestrictList.deleteAll();
        close();
        throw;
    }
    catch (SQLiteException& e)
    {
        vStyAdapters.deleteAll();
        vNspAdapters.deleteAll();
        arRestrictList.deleteAll();
        close();
        Verbose::out(1, "*"); Verbose::out(1, strSQL2);  Verbose::out(1, "*"); throw(Except(e.getMessage()));
    }
    catch (...)
    {
        vStyAdapters.deleteAll();
        vNspAdapters.deleteAll();
        arRestrictList.deleteAll();
        close();
        throw(Except("Unknown exception during Annotation processing. ProbeSetIndex = " + ToStr(uiProbeSetIndex) + ", RowCount = " + ToStr(uiRowCount)));
    }

    vStyAdapters.deleteAll();
    vNspAdapters.deleteAll();
    arRestrictList.deleteAll();
}

void Annotation::newCNProbeSetsHDF5(BaseEngine& engine, AffxSplitArray<CNProbeSet> &vProbeSets)
{
    AffxString strFileName = engine.getOpt("annotation-file");
    vProbeSets.clear();
    m_vParams.clear();

    affymetrix_calvin_parameter::ParameterNameValueType param;
    param.SetName(StringUtils::ConvertMBSToWCS("affymetrix-algorithm-param-genome-version-ucsc"));
    param.SetValueAscii("hg18");
    m_vParams.push_back(param);
    param.SetName(StringUtils::ConvertMBSToWCS("affymetrix-algorithm-param-dbsnp-version"));
    param.SetValueAscii("128");
    m_vParams.push_back(param);
    param.SetName(StringUtils::ConvertMBSToWCS("affymetrix-algorithm-param-netaffx-build"));
    param.SetValueAscii("26.1");
    m_vParams.push_back(param);

    AffxArray<AffxString> arRestrictList;
    bool bRestrictList = (engine.getOpt("probeset-ids") != "");
    if (bRestrictList)
    {
        loadProbeSetNamesFromRestrictList(engine, arRestrictList);
    }

    AffxByteArray ba;
    AffxString str;
    AffxString strProbeSetName;
    bool bSnpInterference = false;
    int iSearchIndex = 0;
    unsigned int uiProbeSetIndex = 0;

    affx::File5_File file5;
    affx::File5_Group* group5 = NULL;
    affx::File5_Tsv* tsv5 = NULL;

    file5.open(strFileName, affx::FILE5_OPEN_RO);
    group5 = file5.openGroup("Annotation", affx::FILE5_OPEN);
    tsv5 = group5->openTsv("Annotation", affx::FILE5_OPEN);
    unsigned int uiRowCount = 0;
    while (tsv5->nextLine() == affx::FILE5_OK)
    {
        uiRowCount++;
    }
    tsv5->close();
    delete tsv5;
    int i = 0;
    double d = 0;
    if (uiRowCount > 0)
    {
        vProbeSets.allocate(uiRowCount);
        tsv5 = group5->openTsv("Annotation", affx::FILE5_OPEN);
        while (tsv5->nextLine() == affx::FILE5_OK)
        {
            tsv5->get(0, 0, &str);

            CNProbeSet* p = vProbeSets.getAt(uiProbeSetIndex);
            if (arRestrictList.getCount() > 0)
            {
                iSearchIndex = arRestrictList.binarySearch(str, 0);
                if (iSearchIndex == -1) {p->setProcess(false);}
                else {p->setProcess(true);}
            }
            p->setProbeSetName(str);
            p->setProcessFlag(((str.startsWith("SNP")) || (str.startsWith("rs")) || (str.startsWith("S-")) ? 3 : 1));
            tsv5->get(0, 1, &i); p->setChromosome(i);
            tsv5->get(0, 2, &i); p->setPosition(i);
            tsv5->get(0, 3, &d); p->setGCContent((float)d);
            if (p->getGCContent() > 1) {p->setGCContent(p->getGCContent() / 49.0);}
            tsv5->get(0, 4, &i);
            p->setPseudoAutosomalRegion(i == 1);
            bSnpInterference = false;
            p->setStyAdapterCode(-1);
            p->setNspAdapterCode(-1);

            uiProbeSetIndex++;
        }
        tsv5->close();
        delete tsv5;
    }

    group5->close();
    delete group5;
    file5.close();

    arRestrictList.deleteAll();
}

void Annotation::loadAnnotationSQLite(BaseEngine& engine, AffxArray<CNProbeSet>& arProbeSets)
{
    int iXChromosome = engine.getOptInt("xChromosome");
    int iYChromosome = engine.getOptInt("yChromosome");
    AffxArray<AffxString> vStyAdapters;
    AffxArray<AffxString> vNspAdapters;

    AffxString strFileName = engine.getOpt("annotation-file");
    if (affx::File5_File::isHdf5file(strFileName)) {return loadAnnotationHDF5(engine, arProbeSets);}
    affymetrix_calvin_parameter::ParameterNameValueType param;
    CNProbeSet objSearch;
    int iSearchIndex = -1;
    AffxString str;
    AffxByteArray ba;
    m_vParams.clear();

    AffxString strSQL1;
    AffxString strSQL2;
    AffxString strSQL3;
    try    {
    open(strFileName);
    SQLiteRecordset rset(*this);

    strSQL2 = "select key, value from information";
    rset.open(strSQL2);
    int iKeyCount = 0;
    while (rset.fetch())
    {
        ba.assign(rset.getString(0));
        ba.replace(" ", "-");
        ba.replace("_", "-");
        ba.replace("netaffx-annotation-netaffx-build", "netaffx-build");
        if (ba.equals("genome-version-ucsc")) {iKeyCount++;}
        if (ba.equals("dbsnp-version")) {iKeyCount++;}
        if (ba.equals("netaffx-build")) {iKeyCount++;}
        param.SetName(StringUtils::ConvertMBSToWCS("affymetrix-algorithm-param-" + ba.toString()));
        param.SetValueAscii(rset.getString(1));
        m_vParams.push_back(param);
    }
    rset.close();
    if (iKeyCount < 3) {throw(Except("The annotation file must contain all three of these parameters: genome-version-ucsc, dbsnp-version, and netaffx-build"));}

    strSQL1 = "select probeSet_ID, chr_id, start, Percent_GC, chrx_par, Chromosome.shortname, process_flag, probe_count,ProbeSet_Type from Annotations, Chromosome where Annotations.chr_id = Chromosome.key ";
    strSQL2 = "select probeSet_ID, chr_id, start, Percent_GC, chrx_par, Chromosome.shortname, enzyme_fragment from Annotations, Chromosome where Annotations.chr_id = Chromosome.key and in_final_list = 1";
    strSQL3 = "select probeSet_ID, chr_id, start, Percent_GC, chrx_par, Chromosome.shortname, enzyme_fragment, process_flag, probe_count, probeset_type, probe_gc_a from Annotations, Chromosome where Annotations.chr_id = Chromosome.key";
    bool bCyto2 = false;
	bool bCytoScanHD = false;
	try {rset.open(strSQL3); bCyto2 = false; bCytoScanHD = true;} 
	catch(...) 
	{
		try {rset.open(strSQL1); bCyto2 = true;} 
		catch(...) 
		{
			rset.open(strSQL2); bCyto2 = false;
		}
	}
    AffxString strChromosome;
    char cChromosome = 0;
    Verbose::out(1, "MAJOR PROGRESS UPDATE: Begin Loading Annotations file.");
    Verbose::progressBegin(1, "Annotation::loadAnnotationSQLite(...) ", (arProbeSets.getCount() / 100000), 1, (arProbeSets.getCount() / 100000));
    int i = 0;
    while (rset.fetch())
    {
        if ((i % 100000) == 0) {Verbose::progressStep(1);}
        i++;
        strChromosome = rset.getString(5);
        if (strChromosome.length() > 0) {cChromosome = strChromosome.getAt(0);}
        str = rset.getString(0);
        str.removeSurroundingQuotes();
        objSearch.setProbeSetName(str);
        iSearchIndex = arProbeSets.binarySearch(objSearch, 0);
                if(iSearchIndex == -1)
                {
                        int i=0;
                         i++;
                }



        if (iSearchIndex == -1) {continue;}

        CNProbeSet* p = arProbeSets.getAt(iSearchIndex);
        p->setProbeSetName(str);
        if (bCyto2)
        {
            p->setProcessFlag(rset.getInteger(6));
            std::string strProbeSetType = rset.getString(8);
            if(strProbeSetType == "SNP")
            {
                p->setTrulySNP(true);
            }
            if(strProbeSetType == "CN")
            { 
                p->setTrulyCN(true);
            }
        }
        else
        {
			if (bCytoScanHD)
			{
				p->setProcessFlag(rset.getInteger(7));
				std::string strProbeSetType = rset.getString(9);
                if (strProbeSetType.empty()) { p->setProcessFlag(0); }
				if(strProbeSetType == "SNP")
				{
					p->setTrulySNP(true);
				}
				if(strProbeSetType == "CN")
			        {	
					p->setTrulyCN(true);
				}
			}
			else
			{
				if (str.startsWith("SNP"))
				{
					p->setProcessFlag(3);
					p->setTrulySNP(true);
				}
				else
				{
					p->setProcessFlag(1);
					p->setTrulyCN(true);
				}
			}
        }

        if ((bCyto2) || (bCytoScanHD))
        {
           if ((p->getChromosome() != iXChromosome) && (cChromosome == 'X'))
           {
		        p->setProcessFlag(0);
           }

           if ((p->getChromosome() != iYChromosome) && (cChromosome == 'Y'))
           {
               p->setProcessFlag(0);
           }
        }

        p->setChromosome(rset.getInteger(1));
        if (cChromosome == 'X') {p->setChromosome(iXChromosome);}
        else if (cChromosome == 'Y') {p->setChromosome(iYChromosome);}
        else if ((cChromosome == 'M') || (cChromosome == 'm')) {p->setChromosome(0);}
        p->setPosition(rset.getInteger(2));
        p->setGCContent((float)rset.getDouble(3));
		if (bCytoScanHD) {p->setReplicateCount(rset.getInteger(8));}
		else {p->setReplicateCount(rset.getInteger(7));}
        p->setPseudoAutosomalRegion(rset.getInteger(4) > 0);
        if (!bCyto2)
        {
            str = rset.getString(6);
            if (bCytoScanHD && str.empty()) { p->setProcessFlag(0); }
            str.removeSurroundingQuotes();
            p->setStyAdapterCode(parseAdapterCode(str, "StyI", vStyAdapters));
            p->setNspAdapterCode(parseAdapterCode(str, "NspI", vNspAdapters));
        }
        else
        {
            p->setStyAdapterCode(-1);
            p->setNspAdapterCode(-1);
        }

        if (bCytoScanHD)
        {
            if (p->getPosition() == 0) { p->setProcessFlag(0); }
            if (p->getGCContent() != p->getGCContent()) { p->setProcessFlag(0); }
        }
    }
    rset.close();
    Verbose::progressEnd(1, "Done.");
    Verbose::out(1, "MAJOR PROGRESS UPDATE: End Loading Annotations file.");

    close();
    }
    catch (SQLiteException& e)
    {
        vStyAdapters.deleteAll();
        vNspAdapters.deleteAll();
        close();
        Verbose::out(1, "*"); Verbose::out(1, strSQL2);  Verbose::out(1, "*");
        throw(Except(e.getMessage()));
    }
    catch (...)
    {
        vStyAdapters.deleteAll();
        vNspAdapters.deleteAll();
        close();
        throw(Except("Unknown exception during Annotation processing."));
    }

    vStyAdapters.deleteAll();
    vNspAdapters.deleteAll();
}

void Annotation::loadAnnotationHDF5(BaseEngine& engine, AffxArray<CNProbeSet>& arProbeSets)
{
    AffxString strFileName = engine.getOpt("annotation-file");
    CNProbeSet objSearch;
    int iSearchIndex = -1;
    AffxString str;
    AffxByteArray ba;
    m_vParams.clear();

    affymetrix_calvin_parameter::ParameterNameValueType param;
    param.SetName(StringUtils::ConvertMBSToWCS("affymetrix-algorithm-param-genome-version-ucsc"));
    param.SetValueAscii("hg18");
    m_vParams.push_back(param);
    param.SetName(StringUtils::ConvertMBSToWCS("affymetrix-algorithm-param-dbsnp-version"));
    param.SetValueAscii("128");
    m_vParams.push_back(param);
    param.SetName(StringUtils::ConvertMBSToWCS("affymetrix-algorithm-param-netaffx-build"));
    param.SetValueAscii("26.1");
    m_vParams.push_back(param);

    affx::File5_File file5;
    affx::File5_Group* group5 = NULL;
    affx::File5_Tsv* tsv5 = NULL;

    file5.open(strFileName, affx::FILE5_OPEN_RO);
    group5 = file5.openGroup("Annotation", affx::FILE5_OPEN);
    tsv5 = group5->openTsv("Annotation", affx::FILE5_OPEN);
    unsigned int uiRowCount = 0;
    while (tsv5->nextLine() == affx::FILE5_OK)
    {
        uiRowCount++;
    }
    tsv5->close();
    delete tsv5;
    int i = 0;
    double d = 0;
    if (uiRowCount > 0)
    {
        tsv5 = group5->openTsv("Annotation", affx::FILE5_OPEN);
        while (tsv5->nextLine() == affx::FILE5_OK)
        {
            tsv5->get(0, 0, &str);
            objSearch.setProbeSetName(str);
            iSearchIndex = arProbeSets.binarySearch(objSearch, 0);
            if (iSearchIndex == -1) {continue;}

            CNProbeSet* p = arProbeSets.getAt(iSearchIndex);
            p->setProbeSetName(str);
            p->setProcessFlag(((str.startsWith("SNP")) || (str.startsWith("rs")) || (str.startsWith("S-")) ? 3 : 1));
            tsv5->get(0, 1, &i); p->setChromosome(i);
            tsv5->get(0, 2, &i); p->setPosition(i);
            tsv5->get(0, 3, &d); p->setGCContent((float)d);
            if (p->getGCContent() > 1) {p->setGCContent(p->getGCContent() / 49.0);}
            tsv5->get(0, 4, &i);
            p->setPseudoAutosomalRegion(i == 1);
            p->setStyAdapterCode(-1);
            p->setNspAdapterCode(-1);
        }
        tsv5->close();
        delete tsv5;
    }

    group5->close();
    delete group5;
    file5.close();

}

char Annotation::parseAdapterCode(AffxString& str, const AffxString& strEnzyme, AffxArray<AffxString>& vAdapters)
{
    AffxString strAdapter;
    char cAdapter = -1;
    int iSearchIndex = -1;
    int iFindIndex = str.indexOf(strEnzyme);
    if (iFindIndex != -1)
    {
        strAdapter = str.substring(iFindIndex + 8, iFindIndex + 21);
        strAdapter = strAdapter.toLowerCase();
        for (unsigned int ui = 0; (ui < strAdapter.length()); ui++)
        {
            if (!strAdapter.isBaseAt(ui))
            {
                if (strAdapter.getAt(ui) == '_') {strAdapter.setAt(ui, 'n');}
                else {Err::errAbort("Adapter in Annotation files cannot be processed.");}
            }
        }
        iSearchIndex = -1;
        for (int i = 0; (i < vAdapters.getCount()); i++)
        {
            if (*vAdapters.getAt(i) == strAdapter)
            {
                iSearchIndex = i;
                break;
            }
        }
        if (iSearchIndex == -1)
        {
            try {strAdapter = strAdapter.reverseComplement();} catch(...) {Err::errAbort("Adapter in Annotation files cannot be processed.");}
            iSearchIndex = -1;
            for (int i = 0; (i < vAdapters.getCount()); i++)
            {
                if (*vAdapters.getAt(i) == strAdapter)
                {
                    iSearchIndex = i;
                    break;
                }
            }
        }
        if (iSearchIndex == -1)
        {
            vAdapters.add(new AffxString(strAdapter));
            cAdapter = (char)(vAdapters.getCount() - 1);
        }
        else
        {
            cAdapter = (char)iSearchIndex;
        }
    }
    return cAdapter;
}

/**
 * @brief Load the probe set names from the restrict list file
 * @param const AffxString& - The signal summary file name
 * @param AffxArray<AffxString>& - The vector to load the string values into
 */
void Annotation::loadProbeSetNamesFromRestrictList(BaseEngine& engine, AffxArray<AffxString>& ar)
{
  AffxString strProbeListFileName = engine.getOpt("probeset-ids");
  if (strProbeListFileName == "") {
    return;
  }
  affx::TsvFile tsv;
  tsv.m_optAutoTrim = true;
  tsv.open(strProbeListFileName);
  std::string strProbeSetName;
  tsv.bind(0,0,&strProbeSetName, affx::TSV_BIND_REQUIRED);
  while (tsv.nextLevel(0) == affx::TSV_OK)  {
    ar.add(new AffxString(strProbeSetName));
  }
  tsv.clear();
  ar.quickSort(0);
  return;
}

int Annotation::loadAnnotationForAdapterTypeNormTran(const std::string& strFileName, ChipLayout& layout, AffxMultiDimensionalArray<char>& mxAdapterProbeTypes, AffxMultiDimensionalArray<char>& mxAdapterNspTypes, AffxMultiDimensionalArray<char>& mxAdapterStyTypes)
{
    int iMaxAdapterCode = 0;
    AffxArray<AffxString> vStyAdapters;
    AffxArray<AffxString> vNspAdapters;

    AffxString strSQL;
    AffxString strSQL1;
    AffxString strSQL2;
    
    try    {
    open(strFileName);
    SQLiteRecordset rset(*this);

    AffxString strProbeSetName;
    AffxString tmpStr;
	bool bCytoScanHD = false;
    strSQL1 = "select snp_interference, probeset_ID, enzyme_fragment from Annotations where in_final_list = 1";
    strSQL2 = "select snp_interference, probeset_ID, enzyme_fragment, probeset_type from Annotations where process_flag > 0";
	try {rset.open(strSQL1); strSQL = strSQL1;} catch(...) {rset.open(strSQL2); strSQL = strSQL2; bCytoScanHD = true;}
    while (rset.fetch())
    {
        bool bSnpInterference = (rset.getInteger(0) == 1);
        strProbeSetName = rset.getString(1);
        strProbeSetName.removeSurroundingQuotes();
        tmpStr = rset.getString(2);
        tmpStr.removeSurroundingQuotes();
        bool bSnp = ((strProbeSetName.find("SNP") != std::string::npos) ? true : false);
		if (bCytoScanHD) {bSnp = ((rset.getString(3) == "SNP") ? true : false);}
        if (layout.containsProbeSet(strProbeSetName))
        {
            ProbeListPacked pList = layout.getProbeListByName(strProbeSetName);
            if (!pList.isNull())
            {
                ProbeSet* ps = ProbeListFactory::asProbeSet(pList);
                for (unsigned int iAtomIndex = 0; (iAtomIndex < ps->atoms.size()); iAtomIndex++)
                {
                    for (unsigned int iProbeIndex = 0; (iProbeIndex < ps->atoms[iAtomIndex]->probes.size()); iProbeIndex++)
                    {
                        if (Probe::isPm(*ps->atoms[iAtomIndex]->probes[iProbeIndex]))
                        {
                            int iIndex = ps->atoms[iAtomIndex]->probes[iProbeIndex]->id;
                            mxAdapterProbeTypes.set(iIndex, (char)(bSnp ? 1 : 2));
                            if (!bSnpInterference)
                            {
                                mxAdapterNspTypes.set(iIndex, parseAdapterCode(tmpStr, "NspI", vNspAdapters));
                                mxAdapterStyTypes.set(iIndex, parseAdapterCode(tmpStr, "StyI", vStyAdapters));
                                iMaxAdapterCode = Max(iMaxAdapterCode, (int)mxAdapterNspTypes.get(iIndex));
                                iMaxAdapterCode = Max(iMaxAdapterCode, (int)mxAdapterStyTypes.get(iIndex));
                            }
                        }
                    }
                }
                delete ps;
            }
        }
    }
    rset.close();

    close();
    }
    catch (SQLiteException& e)
    {
        vStyAdapters.deleteAll();
        vNspAdapters.deleteAll();
        close();
        Verbose::out(1, "*"); Verbose::out(1, strSQL);  Verbose::out(1, "*");
        throw(Except(e.getMessage()));
    }
    catch (...)
    {
        vStyAdapters.deleteAll();
        vNspAdapters.deleteAll();
        close();
        throw(Except("Unknown exception during Annotation processing."));
    }
    vStyAdapters.deleteAll();
    vNspAdapters.deleteAll();

    return iMaxAdapterCode;
}

int Annotation::parseFragmentLength(AffxString& str)
{
    int fragmentLen = -1;
    AffxString strTmp = str;
    int iStartIndex = str.indexOf("//");
    for (int i = 0; i < 2; ++i)
    {
        if (iStartIndex == -1) return -1;
        iStartIndex = str.nextIndexOf("//", iStartIndex);
    }
    if (iStartIndex != -1)
    {
        int iEndIndex = str.nextIndexOf("//", iStartIndex);
        if (iEndIndex == -1) return -1;
        AffxString tmp = str.substring(iStartIndex+2, iEndIndex);
        int startPos = Convert::toInt(tmp.trim());

        iStartIndex = iEndIndex;
        iEndIndex = str.nextIndexOf("//", iStartIndex);
        if (iEndIndex == -1)
            tmp = str.substring(iStartIndex+2);
        else
            tmp = str.substring(iStartIndex+2, iEndIndex);
        int endPos = Convert::toInt(tmp.trim());
        fragmentLen = endPos - startPos;
    }
    return fragmentLen;
}
