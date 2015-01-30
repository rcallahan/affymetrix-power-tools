////////////////////////////////////////////////////////////////
//
// Copyright (C) 2009 Affymetrix, Inc.
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

#include "copynumber/CNTrisomyEngine.h"
//
#include "calvin_files/fusion/src/FusionCHPMultiDataData.h"
#include "calvin_files/utils/src/StringUtils.h"
#include "chipstream/EngineUtil.h"
#include "file/TsvFile/TsvFile.h"
#include "util/PgOptions.h"
//
#include <cstdio>
#include <cstring>
#include <ctime>
#include <iostream>
#include <string>
#include <vector>
//

//
using namespace std;
using namespace affx;
using namespace affymetrix_calvin_data;
using namespace affymetrix_fusion_io;
using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_utilities;

typedef struct _ChromosomeSegmentInformation
{
    int count;
    int start;
    int stop;
    float confidence;
} ChromosomeSegmentInformation;

CNTrisomyEngine::Reg CNTrisomyEngine::reg;

CNTrisomyEngine * CNTrisomyEngine::FromBase(BaseEngine *engine)
{
    if (engine != NULL && engine->getEngineName() == CNTrisomyEngine::EngineName())
        return (CNTrisomyEngine *)engine;
    return NULL;
}

/*
 * Construct the class by defining the options and state.
 */
CNTrisomyEngine::CNTrisomyEngine()
{
    defineOptions();
}

/*
 * Destruct the class
 */
CNTrisomyEngine::~CNTrisomyEngine()
{
}

/*
 * Defines each of the parameters for this engine.
 */
void CNTrisomyEngine::defineOptions()
{
    defineOption("", "in-file", PgOpt::STRING_OPT,
        "Text file specifying the CN results.",
        "");
    defineOption("", "out-file", PgOpt::STRING_OPT,
        "Text file specifying the summary data.",
        "");
    defineOption("", "chr-numbers", PgOpt::STRING_OPT,
        "Chromosomes to report.",
        "");
    defineOption("", "chr-display", PgOpt::STRING_OPT,
        "Chromosome display names.",
        "");
    defineOption("", "syndrome-tests", PgOpt::STRING_OPT,
        "The tests for the syndromes.",
        "");
    defineOption("", "syndrome-names", PgOpt::STRING_OPT,
        "The syndrome names.",
        "");
    defineOption("", "confidence-threshold", PgOpt::DOUBLE_OPT,
        "The threshold for confidence values",
        "0.0");
    defineOption("", "gender-tag", PgOpt::STRING_OPT,
        "The tag in the CHP file specifying the gender call.",
        "Y-gender-call");
    defineOption("", "qc-names", PgOpt::STRING_OPT,
        "The name in the CHP file header specifying the qc metrics.",
        "MAPD,snp-qc");
    defineOption("", "qc-ops", PgOpt::STRING_OPT,
        "The operators for comparing the qc metrics.",
        "le,ge");
    defineOption("", "qc-thrs", PgOpt::STRING_OPT,
        "The thresholds for the qc metrics.",
        "0.27,1.1");
    defineOption("", "test", PgOpt::BOOL_OPT,
        "Flag to run internal unit tests.",
        "false");
}

/*
 * Check the options to make sure that each required one is properly set.
 * If this engine is being run in test mode, then there are no options (other than "test")
 * that are required. For non-test mode, check that the required parameters, or those
 * that do not have a default value are set. If not then abort.
 */
void CNTrisomyEngine::checkOptionsImp()
{
    if (this->getOptBool("test") == false)
    {
        string str;
        str = getOpt("in-file");
        if (str.empty() == true)
            Err::errAbort("Must specify an input file.");
        str = getOpt("out-file");
        if (str.empty() == true)
            Err::errAbort("Must specify an output file.");
        str = getOpt("chr-display");
        if (str.empty() == true)
            Err::errAbort("Must specify chromosomes to analyze.");
        str = getOpt("chr-numbers");
        if (str.empty() == true)
            Err::errAbort("Must specify the display values for each chromosome.");
        str = getOpt("syndrome-names");
        if (str.empty() == true)
            Err::errAbort("Must specify the syndrome names.");
        str = getOpt("syndrome-tests");
        if (str.empty() == true)
            Err::errAbort("Must specify the tests for each syndrome.");
    }
}

/*
 * Convert a chromosome number represented by a string to an integer.
 * This assumes that X and Y are represented by a number, not character.
 */
static u_int8_t StrToChr(const string &chr)
{
    return (u_int8_t) atoi(chr.c_str());
}

/*
 * All up all of the markers for each of the CN states. The map passed in is a map of
 * CN state to a list of segment data (marker count and confidence).
 */
static int TotalMarkers(const map<u_int8_t, list<ChromosomeSegmentInformation> > &chrStates)
{
    int count = 0;
    for (map<u_int8_t, list<ChromosomeSegmentInformation> >::const_iterator it=chrStates.begin(); it!=chrStates.end(); it++)
    {
        for (list<ChromosomeSegmentInformation>::const_iterator segIt=it->second.begin(); segIt!=it->second.end(); segIt++)
            count += segIt->count;
    }
    return count;
}

/*
 * Write the engine options to the output file.
 */
static void writeOptions(std::ofstream &out, Options &options)
{
    vector<string> names;
    vector<PgOpt::PgOptType> types;
    options.getOptionNames(names);
    options.getOptionTypes(types);
    int n = (int) names.size();
    for (int i=0; i<n; i++)
    {
        switch (types[i])
        {
        case PgOpt::BOOL_OPT:
            out << "#%" << names[i] << "=" << (options.getOptBool(names[i])==true ? "true" : "false") << endl;
            break;

        case PgOpt::DOUBLE_OPT:
            out << "#%" << names[i] << "=" << options.getOptDouble(names[i]) << endl;
            break;

        case PgOpt::INT_OPT:
            out << "#%" << names[i] << "=" << options.getOptInt(names[i]) << endl;
            break;

        case PgOpt::STRING_OPT:
            if (options.getOpt(names[i]).empty() == false)
                out << "#%" << names[i] << "=" << options.getOpt(names[i]) << endl;

        case PgOpt::INVALID_OPT:
            break;

        default:
            break;
        }
    }
}

/*
 * Write the QC metrics and the pass/fail flag of the QC to the output file.
 */
static void writeQCMetrics(std::ofstream &out, const string &gender, const vector<string> &qcNames, const vector<float> &qcValues, bool qcTest)
{
    out << "#%gender=" << gender << endl;
    for (int i=0; i<(int)qcNames.size(); i++)
        out << "#%" << qcNames[i] << "=" << qcValues[i] << endl;
    out << "#%qc-test=" << (qcTest == true ? "Pass" : "Fail") << endl;
}

/*
 * Write the summary of each reported chromosome to the output file.
 * The summary consists of, for each chromosome, the number of markers for a given CN state,
 * the fraction of markers of a CN state to the total number of markers, the confidence
 * of that segment. The output also includes the total marker count for the chromosome.
 */
static void writeChromosomeSummary(std::ofstream &out, const vector<string> &chrNumbers, const vector<string> &chrDisplay,
                                   map<u_int8_t, map<u_int8_t, list<ChromosomeSegmentInformation> > > &chrStates)
{
    out.setf(ios::fixed, ios::floatfield);
    out.precision(2);
    int n = (int) chrNumbers.size();
    for (int i=0; i<n; i++)
    {
        const string &chr = chrDisplay[i];
        u_int8_t chrNum = StrToChr(chrNumbers[i]);
        map<u_int8_t, list<ChromosomeSegmentInformation> > &states = chrStates[chrNum];
        int count = TotalMarkers(states);
        for (map<u_int8_t, list<ChromosomeSegmentInformation> >::iterator it=states.begin(); it!=states.end(); it++)
        {
            int cnCount = 0;
            float cnConf = 0.0f;
            for (list<ChromosomeSegmentInformation>::iterator segIt=it->second.begin(); segIt!=it->second.end(); segIt++)
            {
                cnCount += segIt->count;
                cnConf += segIt->confidence * segIt->count;
            }
            out << "#%chr" << chr << "-cn" << (int) it->first << "-count=" << cnCount << endl;
            out << "#%chr" << chr << "-cn" << (int) it->first << "-fraction=" << (float) cnCount / count << endl;
            out << "#%chr" << chr << "-cn" << (int) it->first << "-confidence=" << (float) cnConf / cnCount << endl;
        }
        out << "#%chr" << chr << "-count=" << count << endl;
    }
}

/*
 * Returns true if the chromosome is a sex chromosome.
 */
static bool IsSexChromosome(const string &chr)
{
    return (chr == "X" || chr == "Y");
}

/*
 * Performs a test on a value compared to a threshold.
 * The parameters include the operator to act on the value and its threshold.
 * The operators include:
     ge for greater than or equal to.
     gt for greater than.
     le for less than or equal to.
     lt for less than.
     eq equal to.
     ne not equal to.
*/
static bool TestValue(const string &op, float value, float thr)
{
    if (op == "ge" && value < thr)
    {
        return false;
    }
    else if (op == "gt" && value <= thr)
    {
        return false;
    }
    else if (op == "le" && value > thr)
    {
        return false;
    }
    else if (op == "lt" && value >= thr)
    {
        return false;
    }
    else if (op == "eq" && value != thr)
    {
        return false;
    }
    else if (op == "ne" && value == thr)
    {
        return false;
    }
    return true;
}

/*
 * Split the string into a vector based on the input separator.
 */
static vector<string> split(const string &inputString, const string &sep = ",")
{
    vector<string> tokens;
    size_t substrBegin = 0;
    for (;;)
    {
        size_t substrEnd = inputString.find (sep, substrBegin);
        if (substrEnd == string::npos)
        {
            // No more separators - save what's left, quit.
            string subString = inputString.substr (substrBegin);
            // Avoid returning a null string from a terminating separator or an empty inputString.
            if (! subString.empty())
                tokens.push_back (subString);
            break;
        }
        // Avoid null strings from an initial separator(s).
        if (substrEnd != substrBegin)
            tokens.push_back (inputString.substr (substrBegin, substrEnd - substrBegin) );
        // Continue following the separator
        substrBegin = substrEnd + 1;
    }
    return tokens;
}

/*
 * Extract the QC metrics from the CYCHP file and compare the value to
 * the input thresholds.
 */
static bool TestQcMetrics(FusionCHPMultiDataData *mchp, const vector<string> &qcNames, const vector<string> &qcThrs, const vector<string> &qcOps, vector<float> &qcValues)
{
    bool qcTest = true;
    ParameterNameValueTypeList summary = mchp->GetSummaryParams();
    int n = (int) qcNames.size();
    qcValues.resize(n);
    for (int i=0; i<n; i++)
    {
        float thr = (float)atof(qcThrs[i].c_str());
        const wstring qcName = StringUtils::ConvertMBSToWCS(qcNames[i]);
        for (ParameterNameValueTypeList::iterator it=summary.begin(); it!=summary.end(); it++)
        {
            if (it->GetName() == qcName)
            {
                qcValues[i] = it->GetValueFloat();
                if (TestValue(qcOps[i], qcValues[i], thr) == false)
                    qcTest = false;
                break;
            }
        }
    }
    return qcTest;
}

/*
 * Get the gender from the CYCHP file given a parameter name.
 */
static affx::Gender GetGender(FusionCHPMultiDataData *mchp, const wstring &genderTag)
{
    string female = getGenderString(affx::Female);
    string male = getGenderString(affx::Male);
    ParameterNameValueTypeList summary = mchp->GetSummaryParams();
    for (ParameterNameValueTypeList::iterator it=summary.begin(); it!=summary.end(); it++)
    {
        if (it->GetName() == genderTag)
        {
            if (it->GetValueAscii() == female)
                return affx::Female;
            else if (it->GetValueAscii() == male)
                return affx::Male;
            break;
        }
    }
    return affx::UnknownGender;
}

/*
 * Gather the segments (CN states and confidence values) for each input chromosome.
 * The CN state and confidence are stored in a structure. The structure is stored in a
 * map that is accessed by chromosome and CN state.
 */
static void GatherChromosomeResults(FusionCHPMultiDataData *mchp, float confThr, const vector<string> &chrNumbers, const vector<string> &chrDisplay,
                                    map<u_int8_t, map<u_int8_t, list<ChromosomeSegmentInformation> > > &chrStates)
{
    map<u_int8_t, string> chrs;
    for (int i=0; i<(int) chrNumbers.size(); i++)
        chrs[StrToChr(chrNumbers[i])] = chrDisplay[i];

    int n = mchp->GetEntryCount(SegmentCNMultiDataType);
    ChromosomeSegmentData entry;
    for (int i=0; i<n; i++)
    {
        mchp->GetChromosomeSegmentEntry(SegmentCNMultiDataType, i,  entry);
        if (chrs.find(entry.chr) == chrs.end())
            continue;
        float conf = entry.metrics[1].GetValueFloat();
        if (conf < confThr)
            continue;
        int cnState = (int) (entry.metrics[0].GetValueFloat() + 0.1);
        if (chrStates.find(entry.chr) == chrStates.end())
            chrStates[entry.chr] = map<u_int8_t, list<ChromosomeSegmentInformation> >();
        ChromosomeSegmentInformation info;
        info.count = entry.markerCount;
        info.start = entry.startPosition;
        info.stop = entry.stopPosition;
        info.confidence = conf;
        chrStates[entry.chr][cnState].push_back(info);
    }
}

/*
 * Check if the syndrome is positive given the chromosome segments summary, coverage threshold, gender, and syndrome test.
 * The syntax for the syndrome tests are:
 *   Syndrom tests separated by commas.
 *     Within a syndrome test you have either chromosome or gender tests.
 *     Chromosome tests are indicated by “chr|” and gender tests are indicated by “gender|”.
 *     The definition of the test is after the “|” character.
 *       For gender the value after the “|” character represents the gender from the CYCHP header that must be satisfied.
 *       There is no option for multiple values, that is it will only test for one of the valid genders: male, female, unknown.
 *         For example “gender|male” means the gender must be male. Pretty straight forward.
 *       For chromosome tests the syntax is <chr#>.<eq,lt,le,ge,lt>.<cn state>|<start>|<stop>|<coverage threshold>
 *     You can have multiple tests, separated by a “&” character.
 * For example "chr|13.eq.3|0|100000|0.7&gender|male" means chromosome 13 must have CN state equal to 3 for at least 70% and the gender must be male.
 */
static bool SyndromeTest(const vector<string> &chrNumbers, map<u_int8_t, map<u_int8_t, list<ChromosomeSegmentInformation> > > &chrStates, affx::Gender gender, const string &syndromeTest, float &coverage)
{
    coverage = 0.0f;
    vector<string> tests = split(syndromeTest, "&");
    string genderString = getGenderString(gender);
    for (vector<string>::iterator testsIt=tests.begin(); testsIt!=tests.end(); testsIt++)
    {
        vector<string> test = split(*testsIt, "|");
        if (test[0] == "gender")
        {
            if (test[1] != genderString)
                return false;
        }
        else if (test[0] == "chr")
        {
            vector<string> chrTest = split(test[1], ".");
            int start = atoi(test[2].c_str());
            int stop = atoi(test[3].c_str());
            float covThr = (float)atof(test[4].c_str());
            u_int8_t chr = StrToChr(chrTest[0]);
            string op = chrTest[1];
            u_int8_t cn = StrToChr(chrTest[2]);
            if (chrStates.find(chr) == chrStates.end())
                return false;
            int total = 0;
            int count = 0;
            for (map<u_int8_t, list<ChromosomeSegmentInformation> >::iterator cnIt=chrStates[chr].begin(); cnIt!=chrStates[chr].end(); cnIt++)
            {
                for (list<ChromosomeSegmentInformation>::iterator segIt=cnIt->second.begin(); segIt!=cnIt->second.end(); segIt++)
                {
                    int overlap = min(stop, segIt->stop) - max(start, segIt->start);
                    if (overlap > 0)
                    {
                        total += overlap;
                        if ((op == "eq" && cnIt->first == cn) ||
                            (op == "lt" && cnIt->first < cn) ||
                            (op == "le" && cnIt->first <= cn) ||
                            (op == "gt" && cnIt->first > cn) ||
                            (op == "ge" && cnIt->first >= cn))
                        {
                            count += overlap;
                        }
                    }
                }
            }
            coverage = (float) count / total;
            if (coverage < covThr)
                return false;
        }
    }
    return true;
}

/*
 * Perform unit tests of the static functions.
 */
static void UnitTests()
{
    if (StrToChr("1") != 1)
        Err::errAbort("failed StrToChr1");
    if (StrToChr("24") != 24)
        Err::errAbort("failed StrToChr2");

    ChromosomeSegmentInformation seg={0,0,0,0.0};

    map<u_int8_t, list<ChromosomeSegmentInformation> > states;
    seg.count = 2;
    states[0].push_back(seg);
    seg.count = 3;
    states[1].push_back(seg);
    if (TotalMarkers(states) != 5)
        Err::errAbort("failed TotalMarkers1");
    states.clear();
    if (TotalMarkers(states) != 0)
        Err::errAbort("failed TotalMarkers2");

    if (IsSexChromosome("1") == true)
        Err::errAbort("failed IsSexChromosome1");
    if (IsSexChromosome("2") == true)
        Err::errAbort("failed IsSexChromosome2");
    if (IsSexChromosome("X") == false)
        Err::errAbort("failed IsSexChromosome3");
    if (IsSexChromosome("Y") == false)
        Err::errAbort("failed IsSexChromosome4");

    if (TestValue("eq", 1.0, 1.0) == false)
        Err::errAbort("failed TestValue1");
    if (TestValue("eq", 1.0, 2.0) == true)
        Err::errAbort("failed TestValue2");
    if (TestValue("ne", 1.0, 1.0) == true)
        Err::errAbort("failed TestValue3");
    if (TestValue("ne", 1.0, 2.0) == false)
        Err::errAbort("failed TestValue4");
    if (TestValue("lt", 1.0, 2.0) == false)
        Err::errAbort("failed TestValue5");
    if (TestValue("lt", 3.0, 1.0) == true)
        Err::errAbort("failed TestValue6");
    if (TestValue("le", 1.0, 1.0) == false)
        Err::errAbort("failed TestValue7");
    if (TestValue("le", 1.0, 2.0) == false)
        Err::errAbort("failed TestValue8");
    if (TestValue("le", 2.0, 1.0) == true)
        Err::errAbort("failed TestValue9");

    if (TestValue("gt", 1.0, 2.0) == true)
        Err::errAbort("failed TestValue10");
    if (TestValue("gt", 3.0, 1.0) == false)
        Err::errAbort("failed TestValue11");
    if (TestValue("ge", 1.0, 1.0) == false)
        Err::errAbort("failed TestValue12");
    if (TestValue("ge", 1.0, 2.0) == true)
        Err::errAbort("failed TestValue13");
    if (TestValue("ge", 2.0, 1.0) == false)
        Err::errAbort("failed TestValue14");

    float coverage;
    vector<string> chrNumbers;
    ChromosomeSegmentInformation info = {0, 0, 0, 0.0};

    map<u_int8_t, map<u_int8_t, list<ChromosomeSegmentInformation> > > chrStates;
    affx::Gender gender = affx::Male;
    chrNumbers.push_back("1");
    if (SyndromeTest(chrNumbers, chrStates, gender, "chr|24.eq.1|0|10000000000000|0.70&chr|25.eq.0|0|10000000000000|0.70&gender|female", coverage) == true)
        Err::errAbort("failed SyndromeTest1");
    info.count = 100;
    info.start = 0;
    info.stop = 1000;
    chrStates[1][3].push_back(info);
    if (SyndromeTest(chrNumbers, chrStates, gender, "chr|1.eq.3|0|10000000000000|0.70&gender|male", coverage) == false)
        Err::errAbort("failed SyndromeTest2");
    info.start = 1000;
    chrStates[1][2].push_back(info);
    if (SyndromeTest(chrNumbers, chrStates, gender, "chr|1.eq.3|0|10000000000000|0.70&gender|male", coverage) == true)
        Err::errAbort("failed SyndromeTest3");
    info.count = 100000;
    chrStates[1][4].push_back(info);
    if (SyndromeTest(chrNumbers, chrStates, gender, "chr|1.ge.3|0|10000000000000|0.70&gender|male", coverage) == false)
        Err::errAbort("failed SyndromeTest4");
    if (SyndromeTest(chrNumbers, chrStates, gender, "chr|1.gt.3|0|10000000000000|0.70&gender|male", coverage) == false)
        Err::errAbort("failed SyndromeTest5");
    if (SyndromeTest(chrNumbers, chrStates, gender, "chr|1.le.4|0|10000000000000|0.70&gender|male", coverage) == false)
        Err::errAbort("failed SyndromeTest6");
    if (SyndromeTest(chrNumbers, chrStates, gender, "chr|1.le.2|0|10000000000000|0.70&gender|male", coverage) == true)
        Err::errAbort("failed SyndromeTest7");
}

/*
 * Run the analsis.
 * If the analysis parameters specify unit testing then only run the unit tests.
 * The non-unit test analysis is performed by:
 *   Read the input CYCHP file.
 *     Extract the QC metrics and perform the QC test.
 *     Extract the segment information in a map accessed by chromosome and CN state.
 *     Get the calculated gender.
 *   Create the output file and output the header.
 *   Test each syndrome and output its status.
 */
void CNTrisomyEngine::runImp()
{
    checkOptions();

    // Run the internal tests
    if (getOptBool("test") == true)
    {
        UnitTests();
        return;
    }

    // Read the CHP file.
    string inFile = getOpt("in-file");
    FusionCHPData *chp = FusionCHPDataReg::Read(inFile);
    if (chp == NULL)
        Err::errAbort("Unable to read the CYCHP file.");
    FusionCHPMultiDataData *mchp = FusionCHPMultiDataData::FromBase(chp);

    // Get the QC metrics and apply the thresholds to determine pass/fail
    vector<float> qcValues;
    vector<string> qcNames = split(getOpt("qc-names"));
    vector<string> qcThrs = split(getOpt("qc-thrs"));
    vector<string> qcOps = split(getOpt("qc-ops"));
    bool qcTest = TestQcMetrics(mchp, qcNames, qcThrs, qcOps, qcValues);

    // Gather the CN states.
    map<u_int8_t, map<u_int8_t, list<ChromosomeSegmentInformation> > > chrStates;
    double confThr = getOptDouble("confidence-threshold");
    vector<string> chrNumbers = split(getOpt("chr-numbers"));
    vector<string> chrDisplay = split(getOpt("chr-display"));
    GatherChromosomeResults(mchp, confThr, chrNumbers, chrDisplay, chrStates);

    // Get the gender
    wstring genderTag = StringUtils::ConvertMBSToWCS(getOpt("gender-tag"));
    affx::Gender gender = GetGender(mchp, genderTag);

    // Close the CHP file.
    delete mchp;
    mchp = NULL;

    // Create the output file.
    string outFile = getOpt("out-file");
    ofstream out(outFile.c_str(), ios::out);

    // Output the header.
    writeOptions(out, *this);
    writeQCMetrics(out, getGenderString(gender), qcNames, qcValues, qcTest);
    writeChromosomeSummary(out, chrNumbers, chrDisplay, chrStates);

    // Test each syndrome.
    vector<string> syndromeNames = split(getOpt("syndrome-names"));
    vector<string> syndromeTests = split(getOpt("syndrome-tests"));
    int n = (int) syndromeNames.size();
    float coverage = 0.0f;
    out << "Syndrome" << "\t" << "Status" << "\t" << "Coverage" << endl;
    for (int i=0; i<n; i++)
    {
        bool status = SyndromeTest(chrNumbers, chrStates, gender, syndromeTests[i], coverage);
        out << syndromeNames[i] << "\t";
        if (status == true)
            out << "Positive" << "\t" << coverage << endl;
        else
            out << "Negative" << "\t" << endl;
    }
    out.close();
}
