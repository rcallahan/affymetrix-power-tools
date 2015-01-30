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

#ifdef _MSC_VER
#pragma warning(disable: 4786) // identifier was truncated in the debug information
#endif

//
#include "mas5-stat/src/ExpressionAlgorithmImplementation.h"
//
#include "mas5-stat/src/IntensityFileType.h"
#include "mas5-stat/src/mathLib.h"
#include "mas5-stat/src/pTable.h"
//
#include "calvin_files/utils/src/StringUtils.h"
#include "file/CELFileWriter.h"
#include "file/MSKFileData.h"
#include "rawq/src/RawQ.h"
#include "util/Fs.h"
#include "util/Util.h"
//
#include <cstring>
#include <string.h>


//
const char *ExpressionAbsCallString[]  = {"P", "M", "A", "No Call"};
const char *ExpressionCompCallString[] = {"?", "I", "D", "MI", "MD", "No Change", "No Call"};

using namespace affymetrix_fusion_io;
using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_utilities;
using namespace affxmsk;

//////////////////////////////////////////////////////////////////////

#ifndef FLT_MAX
#define FLT_MAX 3.402823466e+38F 
#endif

#ifndef min
#define min(a,b) (a < b ? a : b)
#endif 

#ifndef max
#define max(a,b) (a > b ? a : b)
#endif 

//////////////////////////////////////////////////////////////////////

double logtwo(double value)
{
	return(log(value)/log(2.0));
}

//////////////////////////////////////////////////////////////////////

double antiLog(double value)
{
	return(pow(2.0, value));
}

//////////////////////////////////////////////////////////////////////

double getDegreeOfFreedom(int nProbeSet)
{
	double df = 0.7 * (nProbeSet - 1);
	if (df < 1.0f)
	{
		df = 1;
	}

	return df;
}

//////////////////////////////////////////////////////////////////////

float vGetT(const int nAtoms, const vector<float> & tTable, float level) 
{
	float t = 0.0f;
	if (nAtoms == 0) {
		return t;
	}
	else if (nAtoms <= (int) tTable.size()) {
		t = tTable[nAtoms - 1];
	}
	else {
		double degreeF = 0.7 * (nAtoms - 1);
		t = (float)tCDFinversed(level, degreeF);
	}
	return t;
}

//////////////////////////////////////////////////////////////////////

int SortFloatsAscending(const void *elem1, const void *elem2)
{
	float v1 = *(float *) elem1;
	float v2 = *(float *) elem2;

	if (v1 > v2)
		return 1;
	else if (v2 > v1)
		return -1;
	else
		return 0;
}

//////////////////////////////////////////////////////////////////////

float computeSquaredDistance(float x1, float y1, float x2, float y2)
{
	float diffx = x1 - x2;
	float diffy = y1 - y2;
	return (float) ( diffx * diffx + diffy * diffy );
}

//////////////////////////////////////////////////////////////////////

float ComputeWeightAtXY(float x, float y, float centerX, float centerY, float smoothFactor)
{
	return 1.0f / (computeSquaredDistance(x , y, centerX, centerY) + smoothFactor);
}

//////////////////////////////////////////////////////////////////////

float trimmedInterpolation(float bwGM, float bLow, float bHigh,	float left, float right)
{
	float x = left;
	if (bwGM >= bHigh)
	{
		x = right;
	}
	else if (bwGM > bLow) 
	{
		float wG = (bwGM - bLow) / (bHigh - bLow);
		x = wG * right + (1.0f - wG) * left;
	}
	return x;
}

//////////////////////////////////////////////////////////////////////

float OneStepBiweightAlgorithm(const vector<float> & x, float c, float epsilon)
{
	if (x.size() == 0)
		return 0.0f;

	float medianValue = median(x);	
	float MAD = medianAbsoluteDeviation(x) * c + epsilon;
	int n = (int) x.size();
	float value=0.0f;
	float weightedSumNumer=0.0f;
	float weightedSumDenom=0.0f;

	for (int i=0; i<n; i++)
	{
		float diff = x[i] - medianValue;
		float u = diff / MAD;
		float uSquare = u * u;
		float oneMinusUSqaure = 1.0f - uSquare;
		if (fabs(u) < 1.0f)
		{
			weightedSumNumer += diff * oneMinusUSqaure * oneMinusUSqaure;
			weightedSumDenom += oneMinusUSqaure * oneMinusUSqaure;
		}
	}
	
	if (weightedSumDenom != 0.0f)
	{
		value = medianValue + weightedSumNumer / weightedSumDenom;
	}

	return value;
}

//////////////////////////////////////////////////////////////////////

float UncertaintyOfEstimate(const vector<float> & x, float c, float epsilon)
{
	if (x.size() == 0)
		return 0.0f;

	float medianValue = median(x);
	float MAD = medianAbsoluteDeviation(x) * c + epsilon;
	int n = (int) x.size();
	float value=0.0f;
	float numer=0.0f;
	float denom=0.0f;
	
	for (int i=0; i<n; i++)
	{
		float diff = x[i] - medianValue;
		float u = diff / MAD;
		float uSquare = u * u;
		float oneMinusUSquare = 1.0f - uSquare;
		if (fabs(u) < 1.0f)
		{
			numer += diff * diff * pow(oneMinusUSquare, 4);
			denom += oneMinusUSquare * (1.0f - 5.0f * uSquare);
		}
	}
	numer = sqrt(numer);
	denom = fabs(denom);

	if (denom != 0.0f)
		value = numer / denom;
	
	return value;
}

//////////////////////////////////////////////////////////////////////

double ComputeEstIntenDiff(const vector<float> & v, float stp) 
{
	double so = 0;
	if (v.size() == 0)
		return so;

	if (v.size() == 1) {
		so = v[0];
	}
	else {
		double meanV = mean(v);
		double s = 0;
		int i;
		for (i = 0; i < (int)v.size(); ++i) {
                        s += (v[i] - meanV) * (v[i] - meanV);
		}
		double std = sqrt(s / (v.size() - 1));	
		double high = meanV + stp * std;
		double low = meanV - stp * std;
		s = 0;
		int ctr = 0;
		for (i = 0; i < (int)v.size(); ++i) {
			if (v[i] <= high && v[i] >= low) {
				ctr ++;
				s += v[i];
			}
		}
		so = s / (double) ctr;
	}
	return so;
}

//////////////////////////////////////////////////////////////////////

void formPTable(const int nT, vector<vector<float> > & vvP) {
	vector<int> mask(nT);
	int i;
	int j;
	int n;
	for (n = 0; n < nT; ++n) {
		mask[n] = 1 << n;
	}
	for (n = 1; n <= nT; ++n) {
		int twoToN = 1 << n;
		vector <int> posRanks(twoToN);
		for (i = 0; i < twoToN; ++i) {
			int sum = 0;
			for (j = 0; j < n; ++j) {
				if (i & mask[j])
					sum += j + 1;
			}
			posRanks[i] = sum;
		}
				
		for (i = 0; i < twoToN; ++i) {
			vector<int> signedRanks(n);
			int w = 0;
			int myCode = 0;
			int j;
			for (j = 0; j < n; ++j) {
				if (i & mask[j]) {
					w += j + 1;
					myCode += 1 << j;
				}
			}
			double tail = 0;
			for (j = 0; j < twoToN; ++j) {
				if (posRanks[j] > w) {
					tail ++;
				}
				else if (posRanks[j] == w) {
					tail += 0.5;
				}
			}
			vvP[n - 1][myCode] =  (float)(tail / (double) twoToN);
		}
	}				
}
/////////////////////////////////////////////////////////////////////////////
// Template Functions
/////////////////////////////////////////////////////////////////////////////
template<class T> ExpResults oneSidedSignRank2(const vector<T> & x, const double alpha) 
{
	// int len = (int) x.size(); // unused pk
	// 1. Ignore all zero differences.
	vector <T> newdiff = x;
	int n = 0;
	int i;
	for (i = 0; i < (int)x.size(); ++i) {
		if (x[i] != 0.0) {
			newdiff[n] = x[i];
			n++;
		}
	}
    if (n == 0) // No non-zero differences.  Output 0.5 as the one-sided p-value and detection is absent.
        return ExpResults (0.5, 0);
	newdiff.resize(n);

	// 2.  Assign integer ranks to the differences.
	vector <T> ranks(n);
	for (i=0; i<n; ++i) {
		ranks[i] = (float)(i + 1);
	}

	// 3. Convert differences to absolute values and sort in ascending order.
	vector <struct ExpResults> absdiff(n);
	for (i = 0; i < n; ++i) {
		absdiff[i].p_value = fabs(newdiff[i]);
		absdiff[i].call = i;
	}
	sort(absdiff.begin(), absdiff.end());

	// 4. If there are ties among absolute differences, all differences in a tie
	//    group are assigned to a rank equal to the average of the integer ranks.
	int nTies = 0;	
    // Avoid cross-platform compatibility problems from approximate floating point arithmetic.
	const double tiny = 2e-09;
	for (i = 0; i < n - 1; ++i) {
		if (fabs(absdiff[i].p_value - absdiff[i+1].p_value) < tiny)
		{
			nTies ++;
			break;
		}
	}
	int tieGroup = 0;
	double doubleVarMod = 0; // modification of variance due to ties.
	if (nTies) {
		i = 0;
		while ( i < n - 1) {
			double initElement = absdiff[i].p_value;
			int tieGroupSize = 1;
			for (int j = i + 1; j < n; ++j) {
				if (fabs (absdiff[j].p_value - initElement) < tiny) {
					tieGroupSize ++;
					if (j == n - 1) {
						i = j;
						tieGroup ++;
						for (int m = j - tieGroupSize + 1; m <= j; ++m) {
							ranks[m] = (2*j - tieGroupSize + 3) / 2.0f;
						}
						doubleVarMod += tieGroupSize *
							((double)tieGroupSize * tieGroupSize - 1);
						break;
					}
				}
				else {
					i = j;
					if (tieGroupSize > 1) {
						tieGroup ++;
						for (int m = j - tieGroupSize; m <= j - 1; ++m) {
							ranks[m] = (2 * j - tieGroupSize + 1) / 2.0f;
						}
						doubleVarMod += tieGroupSize *
							((double)tieGroupSize * tieGroupSize - 1);
					}
					break;
				}
			}
		}
	}
	vector <T> invr(n);
	for (i = 0; i < n; ++i) {
		invr[absdiff[i].call] = ranks[i];
	}

	double w = 0;
	for (i = 0; i < n; ++i) {
		if (newdiff[i] > 0)
			w += invr[i];
	}
	struct ExpResults ans;
	if (n > 11) {
		// Use the asymptotic approximation:
		// S' = [S - n(n+1)/4]/sqrt[n(n+1)(2n+1)/24 - c]
		// where S is the sum of all positive signed ranks
		// and c = sum(b(b*b - 1)/48 is the modification of variance due to ties.
		//
		// p-value = 1 - f(S')
		// where f(S') is the cumulative distribution function of
		// standard normal distribution.
		double dw = w - ((double) n) * (n + 1) / 4.0;
		double denom2 = (((double)n)*(n+1)*(2*n+1) - 0.5*doubleVarMod)/24.0;
		if (denom2 <=0) {
			return ExpResults(0, 0);
		}
		double z = dw / sqrt(denom2);
		ans.p_value = 1 - normalCDF(z);
	}
	else {
		if (nTies == 0)
		{
			int iCode = 0;
			for (i=0; i < n; ++i)
			{
				if (newdiff[i] > 0)
				{
					iCode += 1 << ((int) invr[i] - 1);
				}
			}
			ans.p_value = fGetPValue(n-1, iCode);
		}
		else
		{
			int twoToN = 1 << n;
			vector<int> mask(n);
			for (i = 0; i < n; ++i) {
				mask[i] = 1 << i;
			}
			vector <int> posRanks(twoToN);
			for (i = 0; i < twoToN; ++i) {
				double sum = 0;
				for (int j = 0; j < n; ++j) {
					if (i & mask[j])
						sum += ranks[j];
				}
				posRanks[i] = (int) sum;
			}
			double tail = 0;
			for (i = 0; i < twoToN; ++i) {
				if (posRanks[i] > w) {
					tail ++;
				}
				else if (posRanks[i] == w) {
					tail += 0.5;
				}
			}
			ans.p_value = tail / (double) twoToN;
		}
	}
	ans.call = (ans.p_value < alpha) ? 1 : 0;
	return ans;
}

//////////////////////////////////////////////////////////////////////

template<class T> ExpResults newSignRank(const vector<T> & dif, 
	const double alpha1, const double alpha2) {
	struct ExpResults newPH;
	struct ExpResults oldPH = oneSidedSignRank2(dif, alpha1);
	newPH.p_value = oldPH.p_value;
	if (oldPH.call == 1)
		newPH.call = 2;   
	else if (oldPH.p_value < alpha2)
		newPH.call = 1;
	else if (oldPH.p_value > 1 - alpha1)
		newPH.call = -2;
	else if (oldPH.p_value > 1 - alpha2)
		newPH.call = -1;
	else
		newPH.call = 0;
	return newPH;
}

//////////////////////////////////////////////////////////////////////

template<class T> float mean(const vector<T> & v) {
	if (v.empty()) {
		return -1;
	}
	double sum = 0;
	for (int i = 0; i < (int) v.size(); ++i) {
		sum += v[i];
	}
	return (float)( sum / (double) v.size());
}

//////////////////////////////////////////////////////////////////////

template<class T> float median(const vector<T> & v) {
	vector <T> u = v;
	int len = (int) u.size();
	if (len < 1) {
		return -1;
	}
	sort(u.begin(), u.end());
	int half = len / 2;
	if (len % 2 == 1) {
		return (float)u[half];
	}
	else {
		return ((float)u[half - 1] + (float)u[half]) / 2.0;
	}
}

//////////////////////////////////////////////////////////////////////

template<class T> float stddev(const vector<T> & v) {
	if (v.empty() || v.size() == 1) {
		return -1.0;
	}
	float meanValue = mean(v);
	int n = (int) v.size();
	double sum = 0.0;
	for (int i = 0; i < n; ++i)
	{
		sum += (v[i] - meanValue) * (v[i] - meanValue);
	}
	return sqrt((1.0 / (n - 1)) * sum);
}

//////////////////////////////////////////////////////////////////////

template<class T> float medianAbsoluteDeviation(const vector<T> & x)
{
	float meanValue = median(x);
	int size = (int) x.size();
	vector<T> v(size);
	for (int i=0;i<size;i++)
	{
		v[i] = fabs(x[i] - meanValue);
	}
	return median(v);
}

//////////////////////////////////////////////////////////////////////

template <class T> FloatPair trimMeanAndStd(vector<T> & v, const double p1, const double p2) {
	int total = (int) v.size();
	FloatPair fp(0, 0);
	if (total > 0) {
		sort(v.begin(), v.end());
		int n1 = 0;
		int n2 = (int) floor(total * p2);
		double subtotal = n2;
		if (subtotal > 2) {
			double sum = 0;
			int i;
			for (i = n1; i < n2; ++i) {
				sum += v[i];
			}
			double tMean = sum / subtotal;
			sum = 0;
			for (i = n1; i < n2; ++i) {
                                sum += (v[i] - tMean) * (v[i] - tMean);
			}
			fp.value1 = (float) tMean;
			fp.value2 = (float) sqrt(sum / (subtotal - 1));
		}
		else if (subtotal == 1) {
			fp.value1 = v[n1];
		}
	}
	return fp;
}

//////////////////////////////////////////////////////////////////////

template <class T> double trimMean(const vector<T> & vec, const double p1, const double p2) 
{
	vector<T> whole = vec;
	int total = (int) whole.size();
	if (total == 0)
		return 0.0f;

	sort(whole.begin(), whole.end());

	double dG1 = total * p1;
	double dG2 = total * (1.0 - p2);
	int g1 = (int) floor(dG1);
	int g2 = (int) floor(dG2);
	double r1 = dG1 - g1;
	double r2 = dG2 - g2;
	int last = total - g2 - 1;
	if (last <= 0.0f) { // it is theoretically impossible for last < 0, but
		last = 0; // we add the code here to guarantee proper bounds even if there is any numerical unstability
	}
	double sum = (1.0f - r1) * whole[g1] + (1.0f - r2) * whole[last];
	for (int i = g1 + 1; i < last; ++i) {
		sum += whole[i];
	}
	double subtotal = last - g1 -1;
	if (subtotal <= 0.0f) {
		subtotal = 0.0;
	}
	subtotal += 2.0 - r1 - r2;
	return sum / subtotal;
}

//////////////////////////////////////////////////////////////////////

CExpressionAlgorithmImplementation::CExpressionAlgorithmImplementation()
{
	m_RawQ = 0.0f;
	m_NumResults = 0;
	m_bCompDataExists = false;
	m_ZonesInfo.pZones = NULL;
	m_WriteCell = false;
	m_OutputDirectory = NULL;
}

//////////////////////////////////////////////////////////////////////

CExpressionAlgorithmImplementation::~CExpressionAlgorithmImplementation()
{
	Clear();
}

//////////////////////////////////////////////////////////////////////

void CExpressionAlgorithmImplementation::Clear()
{
	maskedProbes.clear();
	m_ProbeSetNames.clear();
	m_NumResults = 0;
	m_bCompDataExists = false;
    m_ControlInfo.clear();
	try
	{
		m_Cell.Close();
		m_Baseline.Close();
		m_Cdf.Close();
	}
	catch (...)
	{
	}

	vector<AbsStatExpressionProbeSetResultType *>::iterator abs;
	for(abs = m_AbsStatResults.begin(); abs != m_AbsStatResults.end(); ++abs)
		delete[] (*abs);
	m_AbsStatResults.erase(m_AbsStatResults.begin(), m_AbsStatResults.end());
	for(abs = m_BaselineAbsStatResults.begin(); abs != m_BaselineAbsStatResults.end(); ++abs)
		delete[] (*abs);
	m_BaselineAbsStatResults.erase(m_BaselineAbsStatResults.begin(), m_BaselineAbsStatResults.end());
	vector<CompStatExpressionProbeSetResultType *>::iterator comp;
	for(comp = m_CompStatResults.begin(); comp != m_CompStatResults.end(); ++comp)
		delete[] (*comp);
	m_CompStatResults.erase(m_CompStatResults.begin(), m_CompStatResults.end());

	if (m_ZonesInfo.pZones)
		delete [] m_ZonesInfo.pZones;
	m_ZonesInfo.pZones = NULL;
}

//////////////////////////////////////////////////////////////////////

string CExpressionAlgorithmImplementation::GetError() const
{
	return m_Error;
}

//////////////////////////////////////////////////////////////////////

void CExpressionAlgorithmImplementation::SetWriteCell(const bool writeCell)
{
	m_WriteCell = writeCell;
}

//////////////////////////////////////////////////////////////////////

void CExpressionAlgorithmImplementation::SetOutputDirectory(const char* outputDir)
{
	m_OutputDirectory = (char*)outputDir;
}

//////////////////////////////////////////////////////////////////////

void CExpressionAlgorithmImplementation::SetLibPath(const char *libPath)
{
	m_LibPath = libPath;
}

//////////////////////////////////////////////////////////////////////

int CExpressionAlgorithmImplementation::GetNumResults() const
{
	return m_NumResults;
}

//////////////////////////////////////////////////////////////////////

bool CExpressionAlgorithmImplementation::DoesCompDataExists() const
{
	return m_bCompDataExists;
}

//////////////////////////////////////////////////////////////////////

string CExpressionAlgorithmImplementation::GetProbeSetName(int index) const
{
	return m_Cdf.GetProbeSetName(index);
}

//////////////////////////////////////////////////////////////////////

AbsStatExpressionProbeSetResultType *CExpressionAlgorithmImplementation::GetAbsStatResult(int index)
{
	return m_AbsStatResults[index];
}

//////////////////////////////////////////////////////////////////////

AbsStatExpressionProbeSetResultType *CExpressionAlgorithmImplementation::GetBaselineAbsStatResult(int index)
{
	return m_BaselineAbsStatResults[index];
}

//////////////////////////////////////////////////////////////////////

CompStatExpressionProbeSetResultType *CExpressionAlgorithmImplementation::GetCompStatResult(int index)
{
	return (m_bCompDataExists == false ? NULL : m_CompStatResults[index]);
}

//////////////////////////////////////////////////////////////////////

CExpStatAlgSettings &CExpressionAlgorithmImplementation::GetParameters()
{
	return m_Params;
}

//////////////////////////////////////////////////////////////////////

AllZonesInfoType &CExpressionAlgorithmImplementation::GetBackgroundZoneInfo()
{
	return m_ZonesInfo;
}

//////////////////////////////////////////////////////////////////////

bool CExpressionAlgorithmImplementation::IsMasked(FusionCELData *cel, int x, int y)
{
	if (cel->IsMasked(x, y) == true)
		return true;

	if (maskedProbes.find(cel->XYToIndex(x, y)) != maskedProbes.end())
		return true;

	return false;
}

//////////////////////////////////////////////////////////////////////

bool CExpressionAlgorithmImplementation::ReadNormMSKFile()
{
	if (m_Params.NormMethod == CExpStatAlgSettings::NORM_TO_SELECTED_PROBE_SETS)
	{
		// Read the MSK file.
		m_Params.NormGenes.clear();
		CMSKFileData msk;
		msk.SetFileName(m_Params.NormMaskFile.c_str());
		if (msk.Read() == false)
		{
			m_Error = "Unable to read the MSK file.";
			Clear();
			return false;
		}

		// Get the probe setes.
		ProbeSetListConstIt begin;
		ProbeSetListConstIt end;
		ProbeSetListConstIt it;
		msk.GetProbeSetIterators(begin, end);
		for (it=begin; it!=end; ++it)
		{
			m_Params.NormGenes.push_back(*it);
		}
	}
	return true;
}

//////////////////////////////////////////////////////////////////////

bool CExpressionAlgorithmImplementation::ReadScaleMSKFile()
{
	if (m_Params.SFMethod == CExpStatAlgSettings::SCALE_TO_SELECTED_PROBE_SETS)
	{
		// Read the MSK file.
		m_Params.ScaleGenes.clear();
		CMSKFileData msk;
		msk.SetFileName(m_Params.ScaleMaskFile.c_str());
		if (msk.Read() == false)
		{
			m_Error = "Unable to read the MSK file.";
			Clear();
			return false;
		}

		// Get the probe setes.
		ProbeSetListConstIt begin;
		ProbeSetListConstIt end;
		ProbeSetListConstIt it;
		msk.GetProbeSetIterators(begin, end);
		for (it=begin; it!=end; ++it)
		{
			m_Params.ScaleGenes.push_back(*it);
		}
	}
	return true;
}

//////////////////////////////////////////////////////////////////////

bool CExpressionAlgorithmImplementation::ReadProbeMSKFile()
{
	// Check if a MSK file was given.
	maskedProbes.clear();
	if (m_Params.ProbeMaskFile.length() == 0)
		return true;

	// Read the mask file.
	CMSKFileData msk;
	msk.SetFileName(m_Params.ProbeMaskFile.c_str());
	if (msk.Read() == false)
	{
		m_Error = "Unable to read the MSK file.";
		Clear();
		return false;
	}

	// Get each item from the probe pair mask section.
	FusionCDFProbeSetInformation unit;
	FusionCDFProbeGroupInformation blk;
	FusionCDFProbeInformation cel;
	int index;
	int xyindex;
	ProbeSetIndiciesListConstIt begin;
	ProbeSetIndiciesListConstIt end;
	ProbeSetIndiciesListConstIt it;
	msk.GetProbeSetIndiciesIterators(begin, end);
	for (it=begin; it!=end; ++it)
	{
		// Find the corresponding probe set and mask the pair (both PM and MM probe)
		index = m_ProbeSetNames[it->probeSetName];
		m_Cdf.GetProbeSetInformation(index, unit);
		unit.GetGroupInformation(0, blk);
		for (list<int>::const_iterator pairIt=it->indicies.begin(); pairIt!=it->indicies.end(); ++pairIt)
		{
			index = 2*(*pairIt);
			if (index+1 >= blk.GetNumCells())
				continue;

			// Mask the PM probe
			blk.GetCell(index, cel);
			xyindex = m_Cell.XYToIndex(cel.GetX(), cel.GetY());
			maskedProbes.insert(make_pair(xyindex, true));

			// Mask the MM probe.
			blk.GetCell(index+1, cel);
			xyindex = m_Cell.XYToIndex(cel.GetX(), cel.GetY());
			maskedProbes.insert(make_pair(xyindex, true));
		}
	}
	return true;
}

//////////////////////////////////////////////////////////////////////

bool CExpressionAlgorithmImplementation::ReadCelFile(const char *celFile)
{
	// Read the CEL file
	m_Cell.SetFileName(celFile);
	if (m_Cell.Exists() == false)
	{
		m_Error = "The input CEL file does not exist.";
		Clear();
		return false;
	}
	if (m_Cell.Read() == false)
	{
		m_Error = "Unable to read the input CEL file.";
		Clear();
		return false;
	}

    // Modify the saturated intensity if the CEL file was created from the HP scanner.
    if (FromHP(m_Cell) == true)
        m_Params.SaturatedIntensity = m_Params.HPSaturatedIntensity;

	return true;
}

//////////////////////////////////////////////////////////////////////

bool CExpressionAlgorithmImplementation::ReadBaselineCelFile(const char *baseline)
{
	// Read the baseline CEL file
	if (baseline != 0 && baseline[0] != 0)
	{
		m_bCompDataExists = true;
		m_Baseline.SetFileName(baseline);
		if (m_Baseline.Exists() == false)
		{
			m_Error = "The baseline CEL file does not exist.";
			Clear();
			return false;
		}
		if (m_Baseline.Read() == false)
		{
			m_Error = "Unable to read the baseline CEL file.";
			Clear();
			return false;
		}
	}
	return true;
}

//////////////////////////////////////////////////////////////////////

bool CExpressionAlgorithmImplementation::ReadCdfFile(const char *cdfFile)
{
	// Read the CDF file.
	m_Cdf.SetFileName(cdfFile);
	if (m_Cdf.Exists() == false)
	{
		m_Error = "The associated library file (CDF) does not exist.";
		Clear();
		return false;
	}
	if (m_Cdf.Read() == false)
	{
		m_Error = "Unable to read the library file.";
		Clear();
		return false;
	}

	// Check CDF file for all expression units.
	int iunit;
	m_NumResults = m_Cdf.GetHeader().GetNumProbeSets();
	for (iunit=0; iunit<m_NumResults; iunit++)
	{
		if (m_Cdf.GetProbeSetType(iunit) != affxcdf::ExpressionProbeSetType)
		{
			m_Error = "Unable to run algorithm on non-Expression arrays.";
			Clear();
			return false;
		}
	}
	m_NumFeatures = m_Cdf.GetHeader().GetCols() * m_Cdf.GetHeader().GetRows();
	return true;
}

//////////////////////////////////////////////////////////////////////

bool CExpressionAlgorithmImplementation::ReadCdfFile(FusionCELData& cel)
{
	// Read the CDF file.
        string cdfFile = Fs::join(m_LibPath, StringUtils::ConvertWCSToMBS(cel.GetChipType()) + ".CDF");
	m_Cdf.SetFileName(cdfFile.c_str());
	if (m_Cdf.Exists() == false)
	{
		m_Error = "The associated library file (CDF) does not exist.";
		Clear();
		return false;
	}
	if (m_Cdf.Read() == false)
	{
		m_Error = "Unable to read the library file.";
		Clear();
		return false;
	}

	// Check CDF file for all expression units.
	int iunit;
	m_NumResults = m_Cdf.GetHeader().GetNumProbeSets();
	for (iunit=0; iunit<m_NumResults; iunit++)
	{
		if (m_Cdf.GetProbeSetType(iunit) != affxcdf::ExpressionProbeSetType)
		{
			m_Error = "Unable to run algorithm on non-Expression arrays.";
			Clear();
			return false;
		}
	}
	m_NumFeatures = m_Cdf.GetHeader().GetCols() * m_Cdf.GetHeader().GetRows();
	return true;
}

//////////////////////////////////////////////////////////////////////

void CExpressionAlgorithmImplementation::AllocateMemoryForResults()
{
	// Allocate memory for the results
	m_AbsStatResults.resize(m_NumResults, NULL);
	if (m_bCompDataExists == true)
		m_CompStatResults.resize(m_NumResults, NULL);


	// Copy the probe set names and allocate memory.
	string name;
	for (int iunit=0; iunit<m_NumResults; iunit++)
	{
		name = m_Cdf.GetProbeSetName(iunit);
		m_ProbeSetNames[name] = iunit;
		m_AbsStatResults[iunit] = new AbsStatExpressionProbeSetResultType;
		if (m_bCompDataExists == true)
			m_CompStatResults[iunit] = new CompStatExpressionProbeSetResultType;
	}
}

//////////////////////////////////////////////////////////////////////

bool CExpressionAlgorithmImplementation::CheckIntensityData()
{
	// Check if there is any data left.
	if (AllMaskedOut(&m_Cell))
	{
		m_Error = "Unable to compute expression results. The data has all been masked.";
		Clear();
		return false;
	}
	if (IsBlankCellFile(&m_Cell))
	{
		m_Error = "Unable to compute expression results. The data is all zeros.";
		Clear();
		return false;
	}

	// Now check the baseline data.
	if (m_bCompDataExists == true)
	{
		if (AllMaskedOut(&m_Baseline))
		{
			m_Error = "Unable to compute expression results. The baseline data has all been masked.";
			Clear();
			return false;
		}
		if (IsBlankCellFile(&m_Baseline))
		{
			m_Error = "Unable to compute expression results. The baseline data is all zeros.";
			Clear();
			return false;
		}
	}
	return true;
}

//////////////////////////////////////////////////////////////////////

bool CExpressionAlgorithmImplementation::RunStat(const char *celFile, const char *baseline, const char *cdfFile)
{
	// Clear old results.
	Clear();
	m_Error.erase();

	// Read the CEL file
	if (ReadCelFile(celFile) == false)
		return false;

	// Read the baseline CEL file
	if (ReadBaselineCelFile(baseline) == false)
		return false;

	// Read the CDF file.
	if (ReadCdfFile(cdfFile) == false)
		return false;

	// Allocate memory for the results
	AllocateMemoryForResults();

	// Read the MSK file
	if (ReadProbeMSKFile() == false)
		return false;

	// Check the intensity data
	if (CheckIntensityData() == false)
		return false;

	// Read the MSK files for normalization and scaling probe sets.
	if (ReadNormMSKFile() == false || ReadScaleMSKFile() == false)
		return false;

	// Call the algorithms
	ComputeStat();

	return true;
}

//////////////////////////////////////////////////////////////////////

void CExpressionAlgorithmImplementation::ComputeStat()
{
	// Make the expression call.	
	for (int iunit=0; iunit<m_NumResults; iunit++)
		ComputeExpressionStat(&m_Cell, iunit,  m_AbsStatResults[iunit]);

	// Compute the scaled and adjusted intensities for perfect match and mis match cells
	// based on new expression algorithm by Earl Hubbell, Mark Durst.
	vector<vector<bool> > UseAtomE(m_NumResults);
	vector<vector<float> > PMe(m_NumResults);		// Scaled & Adjusted Intensity, SA
	vector<vector<float> > MMe(m_NumResults);		// Scaled & Adjusted Intensity, SA
	vector<vector<FloatPair> > BGe(m_NumResults);		// Smoothing Background. (PM in value1, MM in value2)
	vector<vector<FloatPair> > NoiseE(m_NumResults);	// Noise. (PM in value1, MM in value2)
	vector<float> featureIntene;		// Scaled & Adjusted Intensity by feature

	// Pass in member zone info variable for experiment
	ComputeScaledAdjustedIntensity(&m_Cell, PMe, MMe, UseAtomE, BGe, NoiseE, m_ZonesInfo, featureIntene);

	// Compute the measurement and confidence.
	vector<float> avgMeasurementE(m_NumResults);	// Average measurement.
	vector<vector<float> > PVe(m_NumResults);
	ComputeMeasurement(PMe, MMe, UseAtomE, avgMeasurementE, PVe);

	// Set the measurement to the chip file.
	SetMeasurement(BGe, avgMeasurementE, m_AbsStatResults);

	// Determine the scale factor based on the average measurement.
	m_Params.ScaleFactor = DetermineScaleFactor(m_AbsStatResults);
	if (m_Params.ScaleFactor != 1.0f)
		ModifyIntensities(m_AbsStatResults, m_Params.ScaleFactor);


	// Compute the absolute statistics for the baseline CEL file.
	if (m_bCompDataExists == true)
	{
		// Make the expression call.	
		//vector<AbsStatExpressionProbeSetResultType *> baselineAbsStatResults;
		m_BaselineAbsStatResults.resize(m_NumResults, NULL);
		for (int iunit=0; iunit<m_NumResults; iunit++)
		{
			m_BaselineAbsStatResults[iunit] = new AbsStatExpressionProbeSetResultType;
			ComputeExpressionStat(&m_Baseline, iunit,  m_BaselineAbsStatResults[iunit]);
		}

		// Compute the scaled and adjusted intensities for perfect match and mis match cells
		// based on new expression algorithm by Earl Hubbell, Mark Durst.
		vector<vector<bool> > UseAtomB(m_NumResults);
		vector<vector<float> > PMb(m_NumResults);		// Scaled & Adjusted Intensity, SA
		vector<vector<float> > MMb(m_NumResults);		// Scaled & Adjusted Intensity, SA
		vector<vector<FloatPair> > BGb(m_NumResults);		// Smoothing Background. (PM in value1, MM in value2)
		vector<vector<FloatPair> > NoiseB(m_NumResults);	// Noise. (PM in value1, MM in value2)
		vector<float> featureIntenb;		// Scaled & Adjusted Intensity by feature
		
		// Pass in local zone info variable for baseline
		AllZonesInfoType zonesinfo;
		ComputeScaledAdjustedIntensity(&m_Baseline, PMb, MMb, UseAtomB, BGb, NoiseB, zonesinfo, featureIntenb);
		if (zonesinfo.pZones) delete[] zonesinfo.pZones;

		// Compute the measurement and confidence.
		vector<float> avgMeasurementB(m_NumResults);	// Average measurement.
		vector<vector<float> > PVb(m_NumResults);
		ComputeMeasurement(PMb, MMb, UseAtomB, avgMeasurementB, PVb);

		// Set the measurement to the chip file.
		SetMeasurement(BGb, avgMeasurementB, m_BaselineAbsStatResults);

		// Determine the scale factor based on the average measurement.
		m_Params.BaseScaleFactor = DetermineScaleFactor(m_BaselineAbsStatResults);
		if (m_Params.BaseScaleFactor != 1.0f)
			ModifyIntensities(m_BaselineAbsStatResults, m_Params.BaseScaleFactor);

		// Compute the comparison call.
		CallComparison(m_Params.ScaleFactor, m_Params.BaseScaleFactor, m_BaselineAbsStatResults, BGe, BGb);


		// Determine the normalization factor
		DetermineNormFactor(m_AbsStatResults, m_BaselineAbsStatResults);

		// Normalize the intensities and signal values.
		if (m_Params.NormFactor != 1.0f)
			ModifyIntensities(m_AbsStatResults, m_Params.NormFactor);

		// Scale and Normalize the PV before computing fold change.
		for (int i=0; i<m_NumResults; i++)
		{
			int nProbePair = (int) PVe[i].size();
			for (int j=0; j<nProbePair; j++)
			{
				PVe[i][j] += (float)logtwo(m_Params.ScaleFactor) + (float)logtwo(m_Params.NormFactor);  // PV is log scale.
				PVb[i][j] += (float)logtwo(m_Params.BaseScaleFactor) + (float)logtwo(1.0);
			}
		}


		// Compute the fold change.
		ComputeFoldChange(PVb, PVe, UseAtomB, UseAtomE);
	}

	// Compute some summary statistics
	ComputeAvgMaxMinStdBgNoise(BGe, NoiseE);

	CRawQ rawQ;
	rawQ.SetDefaults();
	rawQ.SetVerticalZones(m_Params.NumberVertZones);
	rawQ.SetHorizontalZones(m_Params.NumberHorZones);
	rawQ.SetPercentBG(m_Params.NumberBGCells);
	m_RawQ = rawQ.ComputeRawQ(m_Cell, m_Cdf);
	if (m_bCompDataExists == true)
		m_BaselineRawQ = rawQ.ComputeRawQ(m_Baseline, m_Cdf);

	ReportCornerControls(affxcdf::CheckerboardPositiveQCProbeSetType);
	ReportCornerControls(affxcdf::CheckerboardNegativeQCProbeSetType);
	ReportCornerControls(affxcdf::CentralCrossPositiveQCProbeSetType);
	ReportCornerControls(affxcdf::CentralCrossNegativeQCProbeSetType);

	// Optionally generate a cel file containing background corrected
	// experiment intensities
	if (m_WriteCell)
	{
		affxcel::CCELFileWriter writeCell;
		writeCell.SetFileName (m_Cell.GetFileName().c_str());
		if (! writeCell.Read())
		{
			cerr << "Unable to read cel file." << endl;
			return;
		}
		// ensure not mapped read-only
		writeCell.EnsureNotMmapped();
		for (int celIx = 0; celIx < m_NumFeatures; ++celIx)
			writeCell.SetIntensity (celIx, featureIntene[celIx]);
		if (strlen (m_OutputDirectory) == 0)
		{
			cerr << "An output directory is required to write a cel file." << endl;
			return;
		}
		string outCell = Fs::join(m_OutputDirectory,Fs::basename(m_Cell.GetFileName()));
		writeCell.SetFileName (outCell.c_str());
		if (! writeCell.WriteTextCel())
		{
			cerr << "Error writing cel file." << endl;
			cerr << writeCell.GetError() << endl;
		}
	}
}

//////////////////////////////////////////////////////////////////////

bool CExpressionAlgorithmImplementation::DetermineRelCallNormFactor(float sfe, float sfb)
{
	FusionCDFProbeSetInformation unit;
	FusionCDFProbeGroupInformation blk;
	FusionCDFProbeInformation pmCell;
	FusionCDFProbeInformation mmCell;

	m_RelNF = 1.0f;
	m_RelNF2 = 1.0f;

	float absNFDiff = (float) fabs(m_Params.NormFactor - 1.0f);
	if ( (m_Params.NormMethod == CExpStatAlgSettings::DEFINED_NORMALIZATION_FACTOR) && (absNFDiff > 0.001f) )
	{
		float baselineSF = sfb;
		float expSF = sfe;
		m_RelNF = m_Params.NormFactor * expSF / baselineSF;
		m_RelNF2 = m_RelNF;
		return true;
	}

	int numAtomsUsed;
	int	numBaseAtomsUsed;
	int numUnitUsed = 0;
	int numBaseUnitUsed = 0;

	// float factorE = sfe; // unused pk
	// float factorB = sfb; // unused pk

	// Loop over all of the units.
	int UnitsPerChip = m_Cdf.GetHeader().GetNumProbeSets();
	vector<vector<float> > diffE(UnitsPerChip);
	vector<vector<float> > diffB(UnitsPerChip);
	vector<float> estIntenDiffE(UnitsPerChip, 0.0f);
	vector<float> estIntenDiffB(UnitsPerChip, 0.0f);
	vector<float> estIntenPosDiffE(UnitsPerChip, 0.0f);
	vector<float> estIntenPosDiffB(UnitsPerChip, 0.0f);
	vector<vector<float> > PMe(UnitsPerChip);
	vector<vector<float> > PMb(UnitsPerChip);
	vector<float> estPMIntenE(UnitsPerChip, 0.0f);
	vector<float> estPMIntenB(UnitsPerChip, 0.0f);

	float stp = m_Params.STP;

	std::string probeSetName;
	for (int iUnit=0; iUnit<UnitsPerChip; iUnit++)
	{
		probeSetName = m_Cdf.GetProbeSetName(iUnit);
		if (m_Params.NormMethod == CExpStatAlgSettings::NORM_TO_SELECTED_PROBE_SETS &&
			this->UseUnitInNormFactor(probeSetName, m_Params.NormGenes) == false)
			continue;

		m_Cdf.GetProbeSetInformation(iUnit, unit);
		unit.GetGroupInformation(0, blk);
		int numCells = blk.GetNumCells();
		int nAtoms = blk.GetNumLists();
		diffE[iUnit].resize(nAtoms);
		diffB[iUnit].resize(nAtoms);
		PMe[iUnit].resize(nAtoms);
		PMb[iUnit].resize(nAtoms);
		numAtomsUsed = 0;
		numBaseAtomsUsed = 0;

		for (int iCell=0; iCell<numCells; iCell+=2)
		{
			blk.GetCell(iCell, pmCell);
			blk.GetCell(iCell+1, mmCell);

			// Check if atom is used in the exp chip.
			if (IsMasked(&m_Cell, pmCell.GetX(), pmCell.GetY()) == false &&
				IsMasked(&m_Cell, mmCell.GetX(), mmCell.GetY()) == false)
			{
				diffE[iUnit][numAtomsUsed] = (
					m_Cell.GetIntensity(pmCell.GetX(), pmCell.GetY()) -
					m_Cell.GetIntensity(mmCell.GetX(), mmCell.GetY()));
				PMe[iUnit][numAtomsUsed] = m_Cell.GetIntensity(pmCell.GetX(), pmCell.GetY());
				numAtomsUsed++;
			}

			// Check if atom is used in the baseline chip.
			if (IsMasked(&m_Baseline, pmCell.GetX(), pmCell.GetY()) == false &&
				IsMasked(&m_Baseline, mmCell.GetX(), mmCell.GetY()) == false)
			{
				diffB[iUnit][numBaseAtomsUsed] = (
					m_Baseline.GetIntensity(pmCell.GetX(), pmCell.GetY()) -
					m_Baseline.GetIntensity(mmCell.GetX(), mmCell.GetY()));
				PMb[iUnit][numBaseAtomsUsed] = m_Baseline.GetIntensity(pmCell.GetX(), pmCell.GetY());
				numBaseAtomsUsed++;
			}
		}

		diffE[iUnit].resize(numAtomsUsed);
		diffB[iUnit].resize(numBaseAtomsUsed);
		PMe[iUnit].resize(numAtomsUsed);
		PMb[iUnit].resize(numBaseAtomsUsed);
		if (numAtomsUsed > 0)
		{
			estIntenDiffE[numUnitUsed] = (float)ComputeEstIntenDiff(diffE[iUnit], stp);
			estIntenPosDiffE[numUnitUsed] = (estIntenDiffE[numUnitUsed] > 0.0f) ? estIntenDiffE[numUnitUsed] : 0.0f;
			estPMIntenE[numUnitUsed] = (float)ComputeEstIntenDiff(PMe[iUnit], stp);
			numUnitUsed++;
		}

		if (numBaseAtomsUsed > 0)
		{
			estIntenDiffB[numBaseUnitUsed] = (float)ComputeEstIntenDiff(diffB[iUnit], stp);
			estIntenPosDiffB[numBaseUnitUsed] = (estIntenDiffB[numBaseUnitUsed] > 0.0f) ? estIntenDiffB[numBaseUnitUsed] : 0.0f;
			estPMIntenB[numBaseUnitUsed] = (float)ComputeEstIntenDiff(PMb[iUnit], stp);
			numBaseUnitUsed++;
		}
	}

	estIntenDiffE.resize(numUnitUsed);
	estIntenPosDiffE.resize(numUnitUsed);
	estIntenDiffB.resize(numBaseUnitUsed);
	estIntenPosDiffB.resize(numBaseUnitUsed);
	estPMIntenE.resize(numUnitUsed);
	estPMIntenB.resize(numBaseUnitUsed);

	float p1 = m_Params.IntensityLowPercent / 100;
	float p2 = 1.0f - m_Params.IntensityHighPercent / 100;
	
	// Calculate the primary normalization factor for PM-MM
	float TMe = (float)trimMean(estIntenDiffE, p1, p2);
	if (TMe <= 0.0f) 
	{
		TMe = (float)trimMean(estIntenPosDiffE, p1, p2);
		if (TMe <= 0.0f)
			return false;
	}

	float TMb = (float)trimMean(estIntenDiffB, p1, p2);
	if (TMb <= 0.0f)
	{
		TMb = (float)trimMean(estIntenPosDiffB, p1, p2);
		if (TMb <= 0.0f)
			return false;
	}

	this->m_RelNF = TMb / TMe;

	// Calculate the primary normalization factor for PM-B
	TMe = (float)trimMean(estPMIntenE, p1, p2);
	TMb = (float)trimMean(estPMIntenB, p1, p2);

	if (TMe > 0.0f)
		this->m_RelNF2 = TMb / TMe;
	else
		return false;

	return true;
}

//////////////////////////////////////////////////////////////////////

void CExpressionAlgorithmImplementation::CallComparison
(
	float sfe,
	float sfb,
	vector<AbsStatExpressionProbeSetResultType *> &baselineAbsStatResults,
	vector<vector<FloatPair> > & BGe,
	vector<vector<FloatPair> > & BGb
)
{
	int numAtomsUsed;
	int numCommonAtoms;
	FusionCDFProbeSetInformation unit;
	FusionCDFProbeGroupInformation blk;
	FusionCDFProbeInformation pmCell;
	FusionCDFProbeInformation mmCell;

	bool bRelCallNormFactor = DetermineRelCallNormFactor(sfe, sfb);

	// float factorE = sfe; // unused pk
	// float factorB = sfb; // unused pk

	// Loop over all of the units.
	int UnitsPerChip = m_Cdf.GetHeader().GetNumProbeSets();
	for (int iUnit=0; iUnit<UnitsPerChip; iUnit++)
	{
		CompStatExpressionProbeSetResultType *pUnitCompResult = GetCompStatResult(iUnit);
		
		AbsStatExpressionProbeSetResultType *pUnitAbsResult = GetAbsStatResult(iUnit);

		// AbsStatExpressionProbeSetResultType *pBaseUnitResult = baselineAbsStatResults[iUnit]; // unused pk

		// Initialize the number of atoms used.
		numAtomsUsed = 0;
		numCommonAtoms = 0;

		// Loop over all atoms to determine the number of probe pairs
		// for which the difference and ratio are significantly greater
		// than the baseline.
		int nAtoms = pUnitAbsResult->NumPairs;
		vector<float> pmB(nAtoms, 0.0);
		vector<float> pmE(nAtoms, 0.0);
		vector<float> mmB(nAtoms, 0.0);
		vector<float> mmE(nAtoms, 0.0);
		vector<float> PMinusBgB(nAtoms, 0.0);
		vector<float> PMinusBgE(nAtoms, 0.0);
		vector<float> diffB(nAtoms, 0.0);
		vector<float> diffE(nAtoms, 0.0);


		m_Cdf.GetProbeSetInformation(iUnit, unit);
		unit.GetGroupInformation(0, blk);
		int numCells = blk.GetNumCells();
		int ctr = 0;
		for (int iCell=0; iCell<numCells; iCell+=2)
		{
			blk.GetCell(iCell, pmCell);
			blk.GetCell(iCell+1, mmCell);

			if (IsMasked(&m_Cell, pmCell.GetX(), pmCell.GetY()) == false &&
				IsMasked(&m_Cell, mmCell.GetX(), mmCell.GetY()) == false &&
				IsMasked(&m_Baseline, pmCell.GetX(), pmCell.GetY()) == false &&
				IsMasked(&m_Baseline, mmCell.GetX(), mmCell.GetY()) == false)
			{
				++numAtomsUsed;

				float expMatchCellI = m_Cell.GetIntensity(pmCell.GetX(), pmCell.GetY());
				float expMismatchCellI = m_Cell.GetIntensity(mmCell.GetX(), mmCell.GetY());
				float baselineMatchCellI = m_Baseline.GetIntensity(pmCell.GetX(), pmCell.GetY());
				float baselineMismatchCellI = m_Baseline.GetIntensity(mmCell.GetX(), mmCell.GetY());

				// Scale back the cell intensities for comparison.
				if (expMatchCellI < m_Params.SaturatedIntensity &&
					expMismatchCellI < m_Params.SaturatedIntensity &&
					baselineMatchCellI < m_Params.SaturatedIntensity &&
					baselineMismatchCellI < m_Params.SaturatedIntensity)
				{
					pmB[ctr] = baselineMatchCellI;
					float mm = baselineMismatchCellI;
					// No need to scale down the background like what we did in MAS since
					// the background is calculated on fly.
					float backGround = BGb[iUnit][iCell/2].value1;
					PMinusBgB[ctr] = pmB[ctr] - backGround;
					diffB[ctr] = pmB[ctr] - mm;

					pmE[ctr] = expMatchCellI;
					mm = expMismatchCellI;
					backGround = BGe[iUnit][iCell/2].value1;
					PMinusBgE[ctr] = pmE[ctr] - backGround;
					diffE[ctr] = pmE[ctr] - mm;
					ctr ++;
				}
			}
		}

		pUnitCompResult->NumCommonPairs = 0;

		if (ctr<=0 || bRelCallNormFactor == false)
		{
			// All intensities are greater than the saturated intensity value.
			// No Call will be returned.
			pUnitCompResult->Change = EXP_NO_COMP_CALL_TYPE;
			pUnitCompResult->ChangePValue = 0;
		}
		else
		{
			pmB.resize(ctr);
			mmB.resize(ctr);
			PMinusBgB.resize(ctr);
			diffB.resize(ctr);
			pmE.resize(ctr);
			mmE.resize(ctr);
			PMinusBgE.resize(ctr);
			diffE.resize(ctr);

			float gamma1 = m_Params.Gamma1H;
			float gamma2 = m_Params.Gamma2H;

			ComputeIntensityDependentSignificances(gamma1, gamma2, pmE, pmB);

			DetermineComparativeCall(pUnitCompResult, diffE, diffB, PMinusBgE, PMinusBgB, m_Params.Gamma1H, m_Params.Gamma2H);

			pUnitCompResult->NumCommonPairs = ctr;
		}
	}
}

//////////////////////////////////////////////////////////////////////

void CExpressionAlgorithmImplementation::DetermineComparativeCall
(
	CompStatExpressionProbeSetResultType *pCompData, 
	vector<float> &diffE,
	vector<float> &diffB,
	vector<float> &PMinusBgE,
	vector<float> &PMinusBgB,
	float gamma1, float gamma2
)
{
	float oneMinusGamma1 = 1.0f - gamma1;
	float oneMinusGamma2 = 1.0f - gamma2;

	vector<float> multipliers(3);
	multipliers[0] = 1.0f / m_Params.Perturbation;
	multipliers[1] = 1.0f;
	multipliers[2] = m_Params.Perturbation;

	vector<float> pValues(3, 0.0f);
	int size = (int) diffE.size();
	for (int i=0; i<3; i++) 
	{
		float nf  = multipliers[i] * m_RelNF;
		float nf2 = multipliers[i] * m_RelNF2;
		int size2 = 2 * size;
		vector<float> nDiff(size2, 0.0f);
		for (int j=0; j < size; j++)
		{
			nDiff[j] = nf * diffE[j] - diffB[j];
		}
		for (int j=size; j < size2; j++)
		{
			nDiff[j] = m_Params.CMultiplier * 
				(nf2 * PMinusBgE[j - size] - PMinusBgB[j - size]);
		}	

		ExpResults result = newSignRank(nDiff, gamma1, gamma2);
		pValues[i] = (float)result.p_value;
	}
				

	int Decision = EXP_NO_CHANGE_CALL_TYPE;

	if (pValues[0] < gamma1 && pValues[1] < gamma1 && pValues[2] < gamma1) 
	{
		Decision = EXP_INCREASE_CALL_TYPE;
	}
	else if (pValues[0] > oneMinusGamma1 && pValues[1] > oneMinusGamma1 && pValues[2] > oneMinusGamma1) 
	{
		Decision = EXP_DECREASE_CALL_TYPE;
	}
	else if (pValues[0] < gamma2 && pValues[1] < gamma2 && pValues[2] < gamma2) 
	{
		Decision = EXP_INCREASE_MODERATE_CALL_TYPE;
	}
	else if (pValues[0] > oneMinusGamma2 && pValues[1] > oneMinusGamma2 && pValues[2] > oneMinusGamma2) 
	{
		Decision = EXP_DECREASE_MODERATE_CALL_TYPE;
	}

	pCompData->Change = Decision;
	float fCriticalPValue = 0.5f;
	if ( pValues[0] < 0.5f && pValues[1] < 0.5f && pValues[2] < 0.5f )
		fCriticalPValue = (float) max(pValues[0], (float) max(pValues[1], pValues[2]));
	else if ( pValues[0] > 0.5f && pValues[1] > 0.5f && pValues[2] > 0.5f )
		fCriticalPValue = (float) min(pValues[0], (float) min(pValues[1], pValues[2]));
	pCompData->ChangePValue = fCriticalPValue;
}

//////////////////////////////////////////////////////////////////////

void CExpressionAlgorithmImplementation::ComputeIntensityDependentSignificances(float &gamma1, float &gamma2,
													 vector<float> & PMe,
													 vector<float> & PMb)
{
	float p1 = m_Params.IntensityLowPercent / 100;
	float p2 = 1.0f - m_Params.IntensityHighPercent / 100;
	
	float TPb = (float)trimMean(PMb, p1, p2);
	float TPe = (float)trimMean(PMe, p1, p2);

	float bc = sqrt(TPb * TPe);
	float bLow = bc * m_Params.BLCoef;
	float bHigh = bc * m_Params.BHCoef;

	int size = (int) PMe.size();
	vector<float> gMean(size);
	for (int i=0; i<size; i++)
	{
		gMean[i] = sqrt(PMb[i] * PMe[i]);
	}
	float c = m_Params.TuningConstantCGammas;
	float epsilon = m_Params.EpsilonGammas;

	float bwGM = OneStepBiweightAlgorithm(gMean, c, epsilon);

	gamma1 = trimmedInterpolation(bwGM, bLow, bHigh, m_Params.Gamma1H, m_Params.Gamma1L);
	gamma2 = trimmedInterpolation(bwGM, bLow, bHigh, m_Params.Gamma2H, m_Params.Gamma2L);
}

//////////////////////////////////////////////////////////////////////

void CExpressionAlgorithmImplementation::ComputeProbeLogRatio
(
	vector<vector<float> > & PVb,
	vector<vector<float> > & PVe,
	vector<vector<bool> > & UseAtomB,
	vector<vector<bool> > & UseAtomE,
	vector<vector<float> > & PLR
)
{
	int nProbeSet = (int) PVe.size();
	for (int i=0; i<nProbeSet; i++)
	{
		int nProbePair = (int) PVe[i].size();
		PLR[i].resize(nProbePair);
		int nProbePairUsed = 0;
		for (int j=0; j<nProbePair; j++)
		{
			if (UseAtomB[i][j] && UseAtomE[i][j])
			{
				PLR[i][nProbePairUsed] = PVe[i][j] - PVb[i][j];
				nProbePairUsed++;
			}
		}
		PLR[i].resize(nProbePairUsed);
	}
}

//////////////////////////////////////////////////////////////////////

void CExpressionAlgorithmImplementation::ComputeAvgLogRatio
(
	vector<vector<float> > & PLR,
	vector<float> & avgLogRatio
)
{
	int nProbeSet = (int) PLR.size();
	float c = m_Params.TuningConstantCAvgLogRatio;
	float epsilon = m_Params.EpsilonAvgLogRatio;
	for (int i=0; i<nProbeSet; i++)
	{
		avgLogRatio[i] = OneStepBiweightAlgorithm(PLR[i], c, epsilon);
	}
}

//////////////////////////////////////////////////////////////////////

void CExpressionAlgorithmImplementation::BuildTTable(float level)
{
	static const int T_TABLE_SIZE = 50;
	m_tTable.clear();
	m_tTable.resize(T_TABLE_SIZE, 0.0f);
	for (int i = 0; i < T_TABLE_SIZE; ++i) 
	{
		double df = 0.7 * i;
		if (df < 1.0) {
			df = 1;
		}
		m_tTable[i] = (float)tCDFinversed(level, df);
	}
}

//////////////////////////////////////////////////////////////////////

void CExpressionAlgorithmImplementation::ComputeFoldChange
(
	vector<vector<float> > & PVb,
	vector<vector<float> > & PVe,
	vector<vector<bool> > & UseAtomB,
	vector<vector<bool> > & UseAtomE
)
{
	// Compute fold change
	int nProbeSet = (int) PVe.size();
	vector<vector<float> > PLR(nProbeSet);
	vector<float> avgLogRatio(nProbeSet);

	ComputeProbeLogRatio(PVb, PVe, UseAtomB, UseAtomE, PLR);
	ComputeAvgLogRatio(PLR, avgLogRatio);

	vector<float> uncertainty(nProbeSet);
	float avgLogRatioLow;
	float avgLogRatioHigh;
	float c = m_Params.TuningConstantCAvgLogRatio;
	float epsilon = m_Params.EpsilonAvgLogRatio;

	// Store the fold changes and confidence intervals to the chip file.
	int UnitsPerChip = m_Cdf.GetHeader().GetNumProbeSets();
	BuildTTable(m_Params.RelConfInterval);

	for (int iUnit=0; iUnit<UnitsPerChip; iUnit++)
	{
		CompStatExpressionProbeSetResultType *pUnitResult =
			(CompStatExpressionProbeSetResultType *) m_CompStatResults[iUnit];

		pUnitResult->SignalLogRatio = avgLogRatio[iUnit];

		int size = (int) PLR[iUnit].size();
		float tvalue = vGetT(size, m_tTable, m_Params.RelConfInterval);		

		uncertainty[iUnit] = UncertaintyOfEstimate(PLR[iUnit], c, epsilon);
		float confidence = tvalue * uncertainty[iUnit];
		avgLogRatioLow = avgLogRatio[iUnit] - confidence;
		avgLogRatioHigh = avgLogRatio[iUnit] + confidence;

		pUnitResult->SignalLogRatioLow = avgLogRatioLow;
		pUnitResult->SignalLogRatioHigh = avgLogRatioHigh;
	}
}

//////////////////////////////////////////////////////////////////////

void CExpressionAlgorithmImplementation::DetermineNormFactor
(
	vector<AbsStatExpressionProbeSetResultType *> &expStatResults,
	vector<AbsStatExpressionProbeSetResultType *> &baselineStatResults
)
{
	int iUnit;
	float normFactor=0.0f;

	// User defined norm factor.
	if (m_Params.NormMethod == CExpStatAlgSettings::DEFINED_NORMALIZATION_FACTOR)
		return;



	// Loop over all of the units.
	int UnitsPerChip = m_Cdf.GetHeader().GetNumProbeSets();
	vector<float> expList(UnitsPerChip);
	vector<float> baseList(UnitsPerChip);
	int unitCount=0;
	int baseunitCount=0;
	string probeSetName;
	for (iUnit=0; iUnit<UnitsPerChip; iUnit++)
	{
		// Determine if unit should be used.
		probeSetName = m_Cdf.GetProbeSetName(iUnit);
		if (m_Params.NormMethod == CExpStatAlgSettings::NORM_TO_ALL_PROBE_SETS ||
			UseUnitInNormFactor(probeSetName, m_Params.NormGenes))
		{
			if (expStatResults[iUnit]->NumUsedPairs != 0)
			{
				expList[unitCount] = expStatResults[iUnit]->Signal;
				unitCount++;
			}

			if (baselineStatResults[iUnit]->NumUsedPairs != 0)
			{
				baseList[baseunitCount] = baselineStatResults[iUnit]->Signal;
				baseunitCount++;
			}
		}
	}

	expList.resize(unitCount);
	baseList.resize(baseunitCount);

	// Compute the ratio of the intensities.
	float p1 = m_Params.IntensityLowPercent / 100;
	float p2 = 1.0f - m_Params.IntensityHighPercent / 100;
	float avgB = (float)trimMean(baseList, p1, p2);
	float avgE = (float)trimMean(expList, p1, p2);

	if (avgE != 0.0f)
	{
		normFactor = avgB / avgE;
		float diffNF = fabs(normFactor - 1.0f);
		if (diffNF < 0.000001f)
			normFactor = 1.0f;
	}
	else
		normFactor = 1.0f;

	// Store the norm factor
	if (unitCount == 0)
		normFactor = 1.0f;

	m_Params.NormFactor = normFactor;
}

//////////////////////////////////////////////////////////////////////

void CExpressionAlgorithmImplementation::ComputeExpressionStat
(
	FusionCELData *pCell,
	int iUnit,
	AbsStatExpressionProbeSetResultType *pUnit
)
{
	FusionCDFProbeSetInformation unit;
	FusionCDFProbeGroupInformation blk;
	FusionCDFProbeInformation pmCell;
	FusionCDFProbeInformation mmCell;

	m_Cdf.GetProbeSetInformation(iUnit, unit);
	unit.GetGroupInformation(0, blk);
	int numCells = blk.GetNumCells();

	// Loop over the cells in the unit
	vector<double> discMinusTau(blk.GetNumLists());
	int ctr = 0;
	int iCell;
	for (iCell=0; iCell<numCells; iCell+=2)
	{
		blk.GetCell(iCell, pmCell);
		blk.GetCell(iCell+1, mmCell);

		if (IsMasked(pCell, pmCell.GetX(), pmCell.GetY()) == false &&
			IsMasked(pCell, mmCell.GetX(), mmCell.GetY()) == false)
		{
			double pmI = pCell->GetIntensity(pmCell.GetX(), pmCell.GetY());
			double mmI = pCell->GetIntensity(mmCell.GetX(), mmCell.GetY());

			// Exclude saturated probe pair.
			if (mmI < m_Params.SaturatedIntensity)
			{
				double sum = pmI + mmI;
				double tau1 = m_Params.Tau;
				if (sum > 0.0f)
					discMinusTau[ctr] = ((pmI - mmI) / sum) - tau1;
				else
					discMinusTau[ctr] =  -tau1;
				ctr ++;
			}
		}
	}

	pUnit->NumPairs = blk.GetNumLists();

	// Compute the absolute call.
	if (ctr <= 0)
	{
		pUnit->Detection = EXP_NO_ABS_CALL_TYPE;
		pUnit->Signal = -1.0f;
		pUnit->DetectionPValue = 0.0f;
		pUnit->NumUsedPairs = 0;
	}
	else
	{
		discMinusTau.resize(ctr);
		ExpResults result;
		float medianValue = 0.0f;
		if (discMinusTau.size() > 0) {
			medianValue = m_Params.Tau + median(discMinusTau);
			result = newSignRank(discMinusTau, (float) m_Params.Alpha1, (float) m_Params.Alpha2);
		}
		else
		{
			result.call = 2;
			result.p_value = 0.0f;
		}
		pUnit->DetectionPValue = (float)result.p_value;
		ComputeAbsoluteCall(pUnit, &result);
		pUnit->NumUsedPairs = ctr;
	}
}

//////////////////////////////////////////////////////////////////////

bool CExpressionAlgorithmImplementation::AllMaskedOut(FusionCELData *pCell)
{
	int nUnits;
	int iUnit;
	int nCells;
	int iCell;
	FusionCDFProbeSetInformation unit;
	FusionCDFProbeGroupInformation blk;
	FusionCDFProbeInformation cell;

	// Get number of units.
	nUnits = m_Cdf.GetHeader().GetNumProbeSets();

	// Loop over all units.
	for (iUnit=0; iUnit<nUnits; ++iUnit)
	{
		if (m_Cdf.GetProbeSetType(iUnit) != affxcdf::ExpressionProbeSetType)
			continue;

		m_Cdf.GetProbeSetInformation(iUnit, unit);
		unit.GetGroupInformation(0, blk);
		nCells = blk.GetNumCells();
		for (iCell=0; iCell<nCells; iCell++)
		{
			blk.GetCell(iCell, cell);
			if (IsMasked(pCell, cell.GetX(), cell.GetY()) == false)
				return false;
		}
	}
	return true;
}

//////////////////////////////////////////////////////////////////////

bool CExpressionAlgorithmImplementation::IsBlankCellFile(FusionCELData *pCell)
{
	const float MINIMUM_CELL_INTENSITY = 0.00001f;
	bool bResult = true;

	int iUnit;
	int nCells;
	int iCell;
	FusionCDFProbeSetInformation unit;
	FusionCDFProbeGroupInformation blk;
	FusionCDFProbeInformation cell;



	// Loop over all of the units.
	int UnitsPerChip = m_Cdf.GetHeader().GetNumProbeSets();
	for (iUnit=0; iUnit<UnitsPerChip; iUnit++)
	{
		if (m_Cdf.GetProbeSetType(iUnit) != affxcdf::ExpressionProbeSetType)
			continue;

		m_Cdf.GetProbeSetInformation(iUnit, unit);
		unit.GetGroupInformation(0, blk);
		nCells = blk.GetNumCells();
		for (iCell=0; iCell<nCells; iCell++)
		{
			blk.GetCell(iCell, cell);
			if (pCell->GetIntensity(cell.GetX(), cell.GetY()) > MINIMUM_CELL_INTENSITY)
			{
				bResult = false;
				break;
			}
		}
	}

	return bResult;
}

//////////////////////////////////////////////////////////////////////

void CExpressionAlgorithmImplementation::ComputeScaledAdjustedIntensity(
											  FusionCELData *pCell,
											  vector<vector<float> > & PM,
											  vector<vector<float> > & MM,
											  vector<vector<bool> > &UseAtom,
											  vector<vector<FloatPair> > & BG,
											  vector<vector<FloatPair> > & Noise,
											  AllZonesInfoType & ZonesInfo,
											  vector<float> & FeatureIntensity)
{
	int iUnit;
	int CellsRemaining;
	int zonex;
	int zoney;
	int NumberZones;
	// int iInten=0; // unused pk

	int *NumberCellsPerZone=NULL;

	// Determine the number of remaining cells in the vertical direction.
	CellsRemaining = m_Cdf.GetHeader().GetCols() % m_Params.NumberVertZones;
	if(CellsRemaining != 0)
	{
		zonex = (m_Cdf.GetHeader().GetCols() + 
				(int)m_Params.NumberVertZones - CellsRemaining) /
				(int)m_Params.NumberVertZones;
	}
	else
	{
		zonex = m_Cdf.GetHeader().GetCols() / m_Params.NumberVertZones;
	}

	// Determine the number of remaining cells in the horizontal direction.
	CellsRemaining = m_Cdf.GetHeader().GetRows() % m_Params.NumberHorZones;
	if(CellsRemaining != 0)
	{
		zoney = (m_Cdf.GetHeader().GetRows() +
				m_Params.NumberHorZones-CellsRemaining) /
				m_Params.NumberHorZones;
	}
	else
	{
		zoney = m_Cdf.GetHeader().GetRows() / m_Params.NumberHorZones;
	}

	// Ensure that there are a match and mismatch cell in the same zone.
	zoney += zoney % 2; //EXPRESSION_ATOMS_PER_CELL;

	// Determine the total number of zones.
	NumberZones = (int)m_Params.NumberVertZones * m_Params.NumberHorZones;

	// Get number of units.
	int NumUnits = m_Cdf.GetHeader().GetNumProbeSets();

	// Allocate space for all atoms intensities and ID's.
	NumberCellsPerZone = new int[NumberZones];

	// Clear arrays.
	memset(NumberCellsPerZone, 0, sizeof(int)*NumberZones);

	// Loop over all units to determine the zone ID's and intensities.
	vector<vector<float> > ZoneCells(NumberZones);
	FusionCDFProbeSetInformation unit;
	FusionCDFProbeGroupInformation blk;
	FusionCDFProbeInformation cell;
	bool bMasked;
	for (iUnit = 0; iUnit < NumUnits; ++iUnit)
	{
		m_Cdf.GetProbeSetInformation(iUnit, unit);

		// Only process expression units.
		if (unit.GetProbeSetType() == affxcdf::ExpressionProbeSetType)
		{
			unit.GetGroupInformation(0, blk);
			int numCells = blk.GetNumCells();

			// Loop over the atoms in the unit
			for (int iCell = 0; iCell < numCells; iCell++)
			{
				blk.GetCell(iCell, cell);
				bMasked = IsMasked(pCell, cell.GetX(), cell.GetY());
				if (bMasked == false)
				{
					int nZone;
					nZone = DetermineZone(cell.GetX(), cell.GetY(), zonex, zoney);
					ZoneCells[nZone].resize(ZoneCells[nZone].size() + 1);
					ZoneCells[nZone][NumberCellsPerZone[nZone]] =
						pCell->GetIntensity(cell.GetX(), cell.GetY());

					if (nZone >= 0 && nZone < NumberZones)
						NumberCellsPerZone[nZone]++;
				}
			}
		}
	}

	// Allocate zones, set smooth factor and set num zones
	ZonesInfo.pZones = new ZoneInfo[NumberZones];
	ZonesInfo.number_zones = NumberZones;
	ZonesInfo.smooth_factor = m_Params.SmoothFactorBG;

	// compute background for each zone
	for (int iZone = 0; iZone < NumberZones; iZone++)
	{
		// Compute the center coordinates of each zone.
		// (x1,y1) is the upper left corner
		// (x2,y2) is the lower right corner
		float x1 = ((int) (iZone % m_Params.NumberVertZones)) * zonex;
		float y1 = ((int) (iZone / m_Params.NumberVertZones)) * zoney;
		float x2 = x1 + zonex;
		float y2 = y1 + zoney;
		ZonesInfo.pZones[iZone].center.x = (x1 + x2) / 2;
		ZonesInfo.pZones[iZone].center.y = (y1 + y2) / 2;

		int iCell=0;
		int numCell = NumberCellsPerZone[iZone];
		ZonesInfo.pZones[iZone].numCell = numCell;

		vector<float> zoneI(numCell);
		vector<int> rank(numCell);

		for (int i = 0; i < numCell; i++)
		{
			float inten = ZoneCells[iZone][i];
			zoneI[iCell] = ModifyIntensitySlightly(inten);
			iCell++;
		}
		float lowBG = 0.0f;
		float highBG = m_Params.NumberBGCells / 100.0f;
		FloatPair fp = trimMeanAndStd(zoneI, lowBG, highBG);
		ZonesInfo.pZones[iZone].background = fp.value1;
		ZonesInfo.pZones[iZone].noise = fp.value2;
	}  
	// End of computing background intensity and noise for each zone.
	// Carried zones and NumberZones as the required information.

	// Compute b(x,y), n(x,y), SA(x,y) which was stored in PM[i][j] and MM[i][j]
	float smoothF = m_Params.SmoothFactorBG;

	FusionCDFProbeInformation pmcell;
	FusionCDFProbeInformation mmcell;
	for (iUnit=0; iUnit<NumUnits; ++iUnit)
	{
		m_Cdf.GetProbeSetInformation(iUnit, unit);

		// Only process expression units.
		if (unit.GetProbeSetType() == affxcdf::ExpressionProbeSetType)
		{
			unit.GetGroupInformation(0, blk);
			int numCells = blk.GetNumCells();
			int nAtoms = blk.GetNumLists();

			PM[iUnit].resize(nAtoms, 0.0);
			MM[iUnit].resize(nAtoms, 0.0);
			UseAtom[iUnit].resize(nAtoms, true);

			BG[iUnit].resize(nAtoms);
			Noise[iUnit].resize(nAtoms);

			// Loop over the atoms in the unit
			int atomCount = 0;
			for (int iCell=0, iAtom=0; iCell<numCells; iCell+=2, ++iAtom)
			{
				blk.GetCell(iCell, pmcell);
				blk.GetCell(iCell+1, mmcell);
				bMasked = (IsMasked(pCell, pmcell.GetX(), pmcell.GetY()) == true) ||
							  (IsMasked(pCell, mmcell.GetX(), mmcell.GetY()) == true);
				if (bMasked == true)
				{
						UseAtom[iUnit][iAtom] = false;
				}

				// Set the background (with smoothing adjustment) for matchcell.
				Coordinate Cellxy;
				Cellxy.x = pmcell.GetX();
				Cellxy.y = pmcell.GetY();
				float WeightedSumBg = 0.0f;
				float WeightedSumNoise = 0.0f;
				float WeightedSumDenom = 0.0f;
				float background = 0.0f;
				float noise = 0.0f;
				int k;
				for (k = 0; k < NumberZones; k++)
				{
					WeightedSumBg    += ComputeWeightAtXY(Cellxy.x, Cellxy.y, ZonesInfo.pZones[k].center.x, ZonesInfo.pZones[k].center.y, smoothF) * ZonesInfo.pZones[k].background;
					WeightedSumNoise += ComputeWeightAtXY(Cellxy.x, Cellxy.y, ZonesInfo.pZones[k].center.x, ZonesInfo.pZones[k].center.y, smoothF) * ZonesInfo.pZones[k].noise;
					WeightedSumDenom += ComputeWeightAtXY(Cellxy.x, Cellxy.y, ZonesInfo.pZones[k].center.x, ZonesInfo.pZones[k].center.y, smoothF); 
				}
				if (WeightedSumDenom != 0.0f)
				{
					background = WeightedSumBg / WeightedSumDenom;
					noise = WeightedSumNoise / WeightedSumDenom;
				}

				BG[iUnit][iAtom].value1 = background;
				Noise[iUnit][iAtom].value1 = noise;

				//float smoothFactor;
				float scaledAdjustedI;
				float inten;
				float modifiedI;

				inten = pCell->GetIntensity((int)Cellxy.x,(int)Cellxy.y);
				modifiedI = ModifyIntensitySlightly(inten);
				scaledAdjustedI = ComputeAdjustedIntensity(modifiedI, background, noise);
				PM[iUnit][iAtom] = scaledAdjustedI;
				
				//////////////////////////////////////////////////////
				// Compute Mis-match intensities
				//////////////////////////////////////////////////////
				Cellxy.x = mmcell.GetX();
				Cellxy.y = mmcell.GetY();
				WeightedSumBg = 0.0f;
				WeightedSumNoise = 0.0f;
				WeightedSumDenom = 0.0f;
				background = 0.0f;
				noise = 0.0f;
				for (k = 0; k < NumberZones; k++)
				{
					WeightedSumBg    += ComputeWeightAtXY(Cellxy.x, Cellxy.y, ZonesInfo.pZones[k].center.x, ZonesInfo.pZones[k].center.y, smoothF) * ZonesInfo.pZones[k].background;
					WeightedSumNoise += ComputeWeightAtXY(Cellxy.x, Cellxy.y, ZonesInfo.pZones[k].center.x, ZonesInfo.pZones[k].center.y, smoothF) * ZonesInfo.pZones[k].noise;
					WeightedSumDenom += ComputeWeightAtXY(Cellxy.x, Cellxy.y, ZonesInfo.pZones[k].center.x, ZonesInfo.pZones[k].center.y, smoothF); 
				}
				if (WeightedSumDenom != 0.0f)
				{
					background = WeightedSumBg / WeightedSumDenom;
					noise = WeightedSumNoise / WeightedSumDenom;
				}

				BG[iUnit][iAtom].value2 = background;
				Noise[iUnit][iAtom].value2 = noise;

				inten = pCell->GetIntensity((int)Cellxy.x,(int)Cellxy.y);
				modifiedI = ModifyIntensitySlightly(inten);
				scaledAdjustedI = ComputeAdjustedIntensity(modifiedI, background, noise);
				MM[iUnit][iAtom] = scaledAdjustedI;

				atomCount++;
			}
			PM[iUnit].resize(atomCount);
			MM[iUnit].resize(atomCount);
		} // if expression type 
	} // for each unit

	// If a cel file with adjusted intensities is requested, loop over all probes.
	if (m_WriteCell)
	{
		FeatureIntensity.resize(m_NumFeatures);
		const int colCount = m_Cdf.GetHeader().GetCols();
		const int rowCount = m_Cdf.GetHeader().GetRows();
		for (int colIx = 0; colIx < colCount; ++colIx)
			for (int rowIx = 0; rowIx < rowCount; ++rowIx)
			{
				float WeightedSumBg = 0.0f;
				float WeightedSumNoise = 0.0f;
				float WeightedSumDenom = 0.0f;
				float background = 0.0f;
				float noise = 0.0f;
				for (int k = 0; k < NumberZones; k++)
				{
					WeightedSumBg    += ComputeWeightAtXY(colIx, rowIx, ZonesInfo.pZones[k].center.x, ZonesInfo.pZones[k].center.y, smoothF) * ZonesInfo.pZones[k].background;
					WeightedSumNoise += ComputeWeightAtXY(colIx, rowIx, ZonesInfo.pZones[k].center.x, ZonesInfo.pZones[k].center.y, smoothF) * ZonesInfo.pZones[k].noise;
					WeightedSumDenom += ComputeWeightAtXY(colIx, rowIx, ZonesInfo.pZones[k].center.x, ZonesInfo.pZones[k].center.y, smoothF); 
				}
				if (WeightedSumDenom != 0.0f)
				{
					background = WeightedSumBg / WeightedSumDenom;
					noise = WeightedSumNoise / WeightedSumDenom;
				}
	
				float inten = pCell->GetIntensity(colIx,rowIx);
				float modifiedI = ModifyIntensitySlightly(inten);
				float scaledAdjustedI = ComputeAdjustedIntensity(modifiedI, background, noise);
				int xyindex = m_Cell.XYToIndex(colIx,rowIx);
				FeatureIntensity[xyindex] = scaledAdjustedI;
			}
	}

	//Clean up
	delete [] NumberCellsPerZone;	
}

//////////////////////////////////////////////////////////////////////

void CExpressionAlgorithmImplementation::ComputeBackGroundZones(
											  FusionCELData *pCell,
											  AllZonesInfoType & ZonesInfo,
											  vector<float> & FeatureIntensity)
{
	int iUnit;
	int CellsRemaining;
	int zonex;
	int zoney;
	int NumberZones;
	// int iInten=0; // unused pk

	int *NumberCellsPerZone=NULL;

	// Determine the number of remaining cells in the vertical direction.
	CellsRemaining = m_Cdf.GetHeader().GetCols() % m_Params.NumberVertZones;
	if(CellsRemaining != 0)
	{
		zonex = (m_Cdf.GetHeader().GetCols() + 
				(int)m_Params.NumberVertZones - CellsRemaining) /
				(int)m_Params.NumberVertZones;
	}
	else
	{
		zonex = m_Cdf.GetHeader().GetCols() / m_Params.NumberVertZones;
	}

	// Determine the number of remaining cells in the horizontal direction.
	CellsRemaining = m_Cdf.GetHeader().GetRows() % m_Params.NumberHorZones;
	if(CellsRemaining != 0)
	{
		zoney = (m_Cdf.GetHeader().GetRows() +
				m_Params.NumberHorZones-CellsRemaining) /
				m_Params.NumberHorZones;
	}
	else
	{
		zoney = m_Cdf.GetHeader().GetRows() / m_Params.NumberHorZones;
	}

	// Ensure that there are a match and mismatch cell in the same zone.
	zoney += zoney % 2; //EXPRESSION_ATOMS_PER_CELL;

	// Determine the total number of zones.
	NumberZones = (int)m_Params.NumberVertZones * m_Params.NumberHorZones;

	// Get number of units.
	int NumUnits = m_Cdf.GetHeader().GetNumProbeSets();

	// Allocate space for all atoms intensities and ID's.
	NumberCellsPerZone = new int[NumberZones];

	// Clear arrays.
	memset(NumberCellsPerZone, 0, sizeof(int)*NumberZones);

	// Loop over all units to determine the zone ID's and intensities.
	vector<vector<float> > ZoneCells(NumberZones);
	FusionCDFProbeSetInformation unit;
	FusionCDFProbeGroupInformation blk;
	FusionCDFProbeInformation cell;
	bool bMasked;
	for (iUnit = 0; iUnit < NumUnits; ++iUnit)
	{
		m_Cdf.GetProbeSetInformation(iUnit, unit);

		// Only process expression units.
		if (unit.GetProbeSetType() == affxcdf::ExpressionProbeSetType)
		{
			unit.GetGroupInformation(0, blk);
			int numCells = blk.GetNumCells();

			// Loop over the atoms in the unit
			for (int iCell = 0; iCell < numCells; iCell++)
			{
				blk.GetCell(iCell, cell);
				bMasked = IsMasked(pCell, cell.GetX(), cell.GetY());
				if (bMasked == false)
				{
					int nZone;
					nZone = DetermineZone(cell.GetX(), cell.GetY(), zonex, zoney);
					ZoneCells[nZone].resize(ZoneCells[nZone].size() + 1);
					ZoneCells[nZone][NumberCellsPerZone[nZone]] =
						pCell->GetIntensity(cell.GetX(), cell.GetY());

					if (nZone >= 0 && nZone < NumberZones)
						NumberCellsPerZone[nZone]++;
				}
			}
		}
	}

	// Allocate zones, set smooth factor and set num zones
	ZonesInfo.pZones = new ZoneInfo[NumberZones];
	ZonesInfo.number_zones = NumberZones;
	ZonesInfo.smooth_factor = m_Params.SmoothFactorBG;

	// compute background for each zone
	for (int iZone = 0; iZone < NumberZones; iZone++)
	{
		// Compute the center coordinates of each zone.
		// (x1,y1) is the upper left corner
		// (x2,y2) is the lower right corner
		float x1 = ((int) (iZone % m_Params.NumberVertZones)) * zonex;
		float y1 = ((int) (iZone / m_Params.NumberVertZones)) * zoney;
		float x2 = x1 + zonex;
		float y2 = y1 + zoney;
		ZonesInfo.pZones[iZone].center.x = (x1 + x2) / 2;
		ZonesInfo.pZones[iZone].center.y = (y1 + y2) / 2;

		int iCell=0;
		int numCell = NumberCellsPerZone[iZone];
		ZonesInfo.pZones[iZone].numCell = numCell;

		vector<float> zoneI(numCell);
		vector<int> rank(numCell);

		for (int i = 0; i < numCell; i++)
		{
			float inten = ZoneCells[iZone][i];
			zoneI[iCell] = ModifyIntensitySlightly(inten);
			iCell++;
		}
		float lowBG = 0.0f;
		float highBG = m_Params.NumberBGCells / 100.0f;
		FloatPair fp = trimMeanAndStd(zoneI, lowBG, highBG);
		ZonesInfo.pZones[iZone].background = fp.value1;
		ZonesInfo.pZones[iZone].noise = fp.value2;
	}

	//Clean up
	delete [] NumberCellsPerZone;	
}


//////////////////////////////////////////////////////////////////////

void CExpressionAlgorithmImplementation::ComputeMeasurement(vector<vector<float> > & PM,
								  vector<vector<float> > & MM,
								  vector<vector<bool> > & UseAtom,
								  vector<float> & avgMeasurement,
								  vector<vector<float> > & PV)
{
	// Compute the typical difference (for each probe set) of log(PM) over log(MM).
	int nProbeSet = (int) PM.size(); //  Number of Probe Set

	// Compute Contrast Value
	vector<vector<float> > CT(nProbeSet);
	ComputeContrastValue(PM, MM, CT);

	// Compute Probe Value
	ComputeProbeValue(PM, CT, PV);

	// Compute Average Log Intensity for probe set i and its confidence interval.
	float c = m_Params.TuningConstantCAvgLogInten;
	float epsilon = m_Params.EpsilonAvgLogInten;

	vector<vector<float> > PVused(nProbeSet);
	GetUsedSet(PV, UseAtom, PVused);

	vector<float> uncertainty(nProbeSet);
	for (int i=0; i<nProbeSet; i++)
	{
		avgMeasurement[i] = OneStepBiweightAlgorithm(PVused[i], c, epsilon);
		avgMeasurement[i] = antiLog(avgMeasurement[i]);
	}
}

//////////////////////////////////////////////////////////////////////

void CExpressionAlgorithmImplementation::GetUsedSet(vector<vector<float> > & InputSet,
						  vector<vector<bool> > & UseAtom,
						  vector<vector<float> > & UsedSet)
{
	int nProbeSet = (int) InputSet.size();
	for (int i=0; i<nProbeSet; i++)
	{
		int size = (int) InputSet[i].size();
		UsedSet[i].resize(size);
		int used = 0;
		for (int j=0; j<size; j++)
		{
			if (UseAtom[i][j])
			{
				UsedSet[i][used] = InputSet[i][j];
				used++;
			}
		}
		UsedSet[i].resize(used);
	}
}

//////////////////////////////////////////////////////////////////////

void CExpressionAlgorithmImplementation::ComputeContrastValue(vector<vector<float> > & PM,
									vector<vector<float> > & MM,
									vector<vector<float> > & CT)
{
	int nProbeSet = (int) PM.size();
	vector<float> SB(nProbeSet);
	ComputeTypicalDifference(PM, MM, SB);
	float ContrastTau = m_Params.ContrastTau;
	for (int i=0; i<nProbeSet; i++)
	{
		int nProbePair = (int) PM[i].size();
		CT[i].resize(nProbePair);
		for (int j=0; j<nProbePair; j++)
		{
			if (MM[i][j] < PM[i][j])
			{
				CT[i][j] = MM[i][j];
			}
			else if ((MM[i][j] >= PM[i][j]) &&
					 (SB[i] > ContrastTau))
			{
				CT[i][j] = PM[i][j] / (float)antiLog(SB[i]);
			}
			else if ((MM[i][j] >= PM[i][j]) &&
					 (SB[i] <= ContrastTau))
			{
				CT[i][j] = PM[i][j] / (float)antiLog(ContrastTau / (1.0 + (ContrastTau - SB[i]) / m_Params.ScaleTau) );
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////

void CExpressionAlgorithmImplementation::ComputeProbeValue(vector<vector<float> > & PM,
								 vector<vector<float> > & CT,
								 vector<vector<float> > & PV)
{
	int nProbeSet = (int) PM.size();
	float delta = m_Params.Delta;
	float correction = 1.0f + m_Params.BiasCorrect;
	for (int i=0; i<nProbeSet; i++)
	{
		int nProbePair = (int) PM[i].size();
		PV[i].resize(nProbePair);
		for (int j=0; j<nProbePair; j++)
		{
			float v = PM[i][j] - CT[i][j];
			if (v < delta)
				PV[i][j] = correction * (float)logtwo(delta);
			else
				PV[i][j] = correction * (float)logtwo(v);
		}
	}
}

//////////////////////////////////////////////////////////////////////

void CExpressionAlgorithmImplementation::ComputeTypicalDifference(vector<vector<float> > & PM,
										vector<vector<float> > & MM,
										vector<float> & SB)
{
	int nProbeSet = (int) PM.size();
	vector<vector<float> > logPM_minus_logMM(nProbeSet);
	float c = m_Params.TuningConstantCSB;
	float epsilon = m_Params.EpsilonSB;

	for (int i=0; i<nProbeSet; i++)
	{
		int nProbePair = (int) PM[i].size(); // Number of Probe Pair in Probe Set i
		logPM_minus_logMM[i].resize(nProbePair);
		for (int j=0; j<nProbePair; j++)
		{
			logPM_minus_logMM[i][j] = logtwo(PM[i][j]) - logtwo(MM[i][j]);
		}
		SB[i] = OneStepBiweightAlgorithm(logPM_minus_logMM[i], c, epsilon);
	}
}

//////////////////////////////////////////////////////////////////////

void CExpressionAlgorithmImplementation::ModifyIntensities
(
	vector<AbsStatExpressionProbeSetResultType *> statResults,
	float factor
)
{
	// Loop over all of the units.
	int UnitsPerChip = m_Cdf.GetHeader().GetNumProbeSets();
	for (int iUnit=0; iUnit<UnitsPerChip; iUnit++)
	{
		// Get the units.
		AbsStatExpressionProbeSetResultType *pUnitResult = statResults[iUnit];
		pUnitResult->Signal = (pUnitResult->Signal * factor);
	}
}

//////////////////////////////////////////////////////////////////////

float CExpressionAlgorithmImplementation::DetermineScaleFactor(vector<AbsStatExpressionProbeSetResultType *> &statResults)
{
	int iUnit;
	float avg = 0.0f;

	
	// User defined norm factor is already stored in the algorithm parameters.
	if (m_Params.SFMethod == CExpStatAlgSettings::DEFINED_SCALING_FACTOR)
		return m_Params.ScaleFactor;
 
	// Loop over all of the units.
	int UnitsPerChip = m_Cdf.GetHeader().GetNumProbeSets();
	int unitCount=0;
	vector<float> intensityList(UnitsPerChip);
	string probeSetName;
	for (iUnit=0; iUnit<UnitsPerChip; iUnit++)
	{
		// Get the units.
		AbsStatExpressionProbeSetResultType *pUnitResult = statResults[iUnit];

		// find signal only for the used genes.
		probeSetName = m_Cdf.GetProbeSetName(iUnit);
		if (m_Params.SFMethod == CExpStatAlgSettings::SCALE_TO_ALL_PROBE_SETS ||
			UseUnitInNormFactor(probeSetName, m_Params.ScaleGenes))
		{
			// Use Measurement to do scale factor
			if (pUnitResult->NumUsedPairs != 0)
			{
				intensityList[unitCount] = pUnitResult->Signal;
				unitCount++;
			}
		}
	}

	intensityList.resize(unitCount);

	// Compute the trimMean
	float p1 = m_Params.IntensityLowPercent / 100;
	float p2 = 1.0f - m_Params.IntensityHighPercent / 100;
	avg = trimMean(intensityList, p1, p2);

	// Store the scale factor
	float sf=1.0f;
	if (unitCount && avg != 0.0f)
		sf = m_Params.TGT / avg;

	//Check for the validity of SF
	if(sf <= 0)
	{
		sf = 1.0f;
	}
	return sf;
}

//////////////////////////////////////////////////////////////////////

bool CExpressionAlgorithmImplementation::UseUnitInNormFactor(string &probeSetName,
								   list<string> &MaskedGenes)
{
	bool UseUnit=false;
	list<string>::iterator iter;
	for( iter = MaskedGenes.begin(); iter != MaskedGenes.end(); ++iter)
	{
		if (probeSetName == *iter)
		{
			UseUnit = true;
			break;		
		}
	}
	return UseUnit;
}

//////////////////////////////////////////////////////////////////////

void CExpressionAlgorithmImplementation::ComputeAbsoluteCall(AbsStatExpressionProbeSetResultType *unit, ExpResults *res)
{
	// Determine if the gene is present.
	if( res->call == 2)
		unit->Detection = EXP_PRESENT_CALL_TYPE;

	// Determine if marginal
	else if ( res->call == 1)
		unit->Detection = EXP_MARGINAL_CALL_TYPE;

	// Otherwise the gene is absent.
	else
		unit->Detection = EXP_ABSENT_CALL_TYPE;
}

//////////////////////////////////////////////////////////////////////

float CExpressionAlgorithmImplementation::ModifyIntensitySlightly(float intensity)
{
	return max(intensity, m_Params.Epsilon);
}

//////////////////////////////////////////////////////////////////////

float CExpressionAlgorithmImplementation::ComputeAdjustedIntensity(float intensity, float background, float noise)
{
	float factoredNoise = (float) noise * m_Params.NoiseFrac;
	float diff = intensity - background;
	float adjustedI = max(max(diff, factoredNoise), 0.5f);
	// AlexC - 1/22/01
	// Code Comments:
	// if too frequent substitution of the noise value, alert might generated.
	// Refer to the page 4 of the Background and Spatial Variation Adjustement spec., 
	// under eq(16), it said that 
	// "Production software should probably alert the user if the noise value is being 
	//  substituted too frequently (indicating that too much data is below the noise level), 
	//  but an appropriate threshold value is not at present known."
	return adjustedI;
}

//////////////////////////////////////////////////////////////////////

int CExpressionAlgorithmImplementation::DetermineZone(int cellx, int celly, int zonex, int zoney)
{
	float fZx = (float) cellx / (float) zonex;
	float fZy = (float) celly / (float) zoney;

	int Zx = (int) floor(fZx);
	int Zy = (int) floor(fZy);

	int zoneID = Zx + Zy * m_Params.NumberVertZones;
	return zoneID;
}

//////////////////////////////////////////////////////////////////////

void CExpressionAlgorithmImplementation::SetMeasurement
(
	vector<vector<FloatPair> > &BG,
	vector<float> &avgMeasurement,
	vector<AbsStatExpressionProbeSetResultType *> statResults
)
{
	// FusionCDFProbeSetInformation *unit=NULL; // unused pk
	// FusionCDFProbeGroupInformation *blk=NULL; // unused pk
	// FusionCDFProbeInformation *cell=NULL; // unused pk

	int UnitsPerChip = m_Cdf.GetHeader().GetNumProbeSets();
	for (int iUnit=0; iUnit<UnitsPerChip; iUnit++)
	{
		AbsStatExpressionProbeSetResultType *pUnitResult = statResults[iUnit];
		if (pUnitResult->Detection == EXP_NO_ABS_CALL_TYPE)
			pUnitResult->Signal = 0.0f;
		else
			pUnitResult->Signal = avgMeasurement[iUnit];
	}
}

//////////////////////////////////////////////////////////////////////

float CExpressionAlgorithmImplementation::ComputeCellBackground
(
	int x,
	int y,
	ZoneInfo *zones,
	int zoneSize,
	float smoothFactorBg
)
{
	float bkg = 0.0f;
	if (zones != NULL)
	{
		float weightedSumBg = 0.0f;
		float weight = 0.0f;
		float weightedSumDenom = 0.0f;
		float background = 0.0f;

		for (int k = 0; k < zoneSize; k++)
		{
			weight = ComputeWeightAtXY(x, y, zones[k].center.x, zones[k].center.y, smoothFactorBg);
			weightedSumBg += weight * zones[k].background;
			weightedSumDenom += weight; 
		}
		if (weightedSumDenom != 0.0f)
		{
			background = weightedSumBg / weightedSumDenom;
		}
		bkg = background;
	}

	return bkg;
}

//////////////////////////////////////////////////////////////////////

void CExpressionAlgorithmImplementation::ComputeAvgMaxMinStdBgNoise
(
	vector<vector<FloatPair> > & BG, 
	vector<vector<FloatPair> > & Noise
)
{
	float maxBg = 0.0f;
	float minBg = 0.0f;
	float avgBg = 0.0f;
	float stddevBg = 0.0f;
	float maxNoise = 0.0f;
	float minNoise = 0.0f;
	float avgNoise = 0.0f;
	float stddevNoise = 0.0f;

	vector<float> bgList;
	vector<float> noiseList;
	int count=0;

	int UnitsPerChip = m_Cdf.GetHeader().GetNumProbeSets();
	for (int iUnit=0; iUnit<UnitsPerChip; iUnit++)
	{
		int nAtoms = (int) BG[iUnit].size();
		bgList.resize(count + nAtoms * 2);
		noiseList.resize(count + nAtoms * 2);
		for (int iAtom=0; iAtom<nAtoms; iAtom++)
		{
			bgList[count] = BG[iUnit][iAtom].value1;
			noiseList[count] = Noise[iUnit][iAtom].value1;
			count++;

			bgList[count] = BG[iUnit][iAtom].value2;
			noiseList[count] = Noise[iUnit][iAtom].value2;
			count++;
		}
	}

	// Compute the min/max/avg/stdev values
	maxBg = *max_element(bgList.begin(), bgList.end());
	minBg = *min_element(bgList.begin(), bgList.end());
	maxNoise = *max_element(noiseList.begin(), noiseList.end());
	minNoise = *min_element(noiseList.begin(), noiseList.end());
	avgBg = mean(bgList);
	avgNoise = mean(noiseList);
	stddevBg = stddev(bgList);
	stddevNoise = stddev(noiseList);


	// Store the noise stats
	m_NoiseStats.avg = avgNoise;
	m_NoiseStats.stdv = stddevNoise;
	m_NoiseStats.min = minNoise;
	m_NoiseStats.max = maxNoise;


	// Store the background stats
	m_BgStats.avg = avgBg;
	m_BgStats.stdv = stddevBg;
	m_BgStats.min = minBg;
	m_BgStats.max = maxBg;
}

//////////////////////////////////////////////////////////////////////

void CExpressionAlgorithmImplementation::ReportCornerControls(affxcdf::GeneChipQCProbeSetType qcType)
{
	// Get the desired qc probe set.
	FusionCDFQCProbeSetInformation qcUnit;
	m_Cdf.GetQCProbeSetInformation(qcType, qcUnit);
	if (qcUnit.GetNumCells() == 0)
		return;

	// Get the average noise.
	const float NOISE_FRAC = 0.5f;
	float avgNoise = m_NoiseStats.avg * NOISE_FRAC;

	// Compute the average intensity
	int size = qcUnit.GetNumCells();
	float avgIntensity = 0.0f;
	FusionCDFQCProbeInformation qcCell;
	for (int iCell=0; iCell<size; iCell++)
	{
		qcUnit.GetProbeInformation(iCell, qcCell);
		float fQCBg = ComputeCellBackground(qcCell.GetX(), qcCell.GetY(),
			m_ZonesInfo.pZones, m_ZonesInfo.number_zones, m_Params.SmoothFactorBG);
		float fQCIntensity = m_Cell.GetIntensity(qcCell.GetX(), qcCell.GetY());
		float fDelta = fQCIntensity - fQCBg;
		if (fDelta <= avgNoise)
			avgIntensity += avgNoise;
		else
			avgIntensity += fDelta;
	}

	// Compute average intensity.
	if (size)
	{
		avgIntensity /= size;
		ControlInformationType info;
		info.avg = avgIntensity;
		info.count = size;
		info.qcType = qcType;
		m_ControlInfo.push_back(info);
	}
}
