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

#if !defined(AFX_EXPRESSIONALGORITHMIMPLEMENTATION_H__F70115DD_BD28_4226_9D03_61ACFD937B05__INCLUDED_)
#define AFX_EXPRESSIONALGORITHMIMPLEMENTATION_H__F70115DD_BD28_4226_9D03_61ACFD937B05__INCLUDED_

/////////////////////////////////////////////////////////////////////////////

#include "mas5-stat/src/ExpStatAlgSettings.h"
#include "mas5-stat/src/ExpTmpl.h"
//
#include "calvin_files/fusion/src/FusionCDFData.h"
#include "calvin_files/fusion/src/FusionCELData.h"
//
#include <cstring>
#include <list>
#include <map>
#include <string>
#include <vector>
//

using namespace std;
using namespace affymetrix_fusion_io;

/////////////////////////////////////////////////////////////////////////////

#define EXP_PRESENT_CALL_TYPE 0
#define EXP_MARGINAL_CALL_TYPE 1
#define EXP_ABSENT_CALL_TYPE 2
#define EXP_NO_ABS_CALL_TYPE 3
// this should be extern
extern const char *ExpressionAbsCallString[];

//////////////////////////////////////////////////////////////////////

typedef struct _AbsStatExpressionProbeSetResultType
{
	float DetectionPValue;
	float Signal; 
	unsigned short NumPairs;
	unsigned short NumUsedPairs;
	unsigned char Detection;
	_AbsStatExpressionProbeSetResultType *operator=(_AbsStatExpressionProbeSetResultType &src)
	{
		memcpy(this, &src, sizeof(_AbsStatExpressionProbeSetResultType));
		return this;
	};
} AbsStatExpressionProbeSetResultType;

////////////////////////////////////////////////////////////////////

#define EXP_UNKNOWN_COMP_CALL_TYPE 0
#define EXP_INCREASE_CALL_TYPE 1
#define EXP_DECREASE_CALL_TYPE 2
#define EXP_INCREASE_MODERATE_CALL_TYPE 3
#define EXP_DECREASE_MODERATE_CALL_TYPE 4
#define EXP_NO_CHANGE_CALL_TYPE 5
#define EXP_NO_COMP_CALL_TYPE 6
//
extern const char *ExpressionCompCallString[];

////////////////////////////////////////////////////////////////////

typedef struct _CompStatExpressionProbeSetResultType
{
	float ChangePValue; 
	float SignalLogRatio; 
	float SignalLogRatioLow; 
	float SignalLogRatioHigh; 
	unsigned short NumCommonPairs;
	unsigned char Change;
	_CompStatExpressionProbeSetResultType *operator=(_CompStatExpressionProbeSetResultType &src)
	{
		memcpy(this, &src, sizeof(_CompStatExpressionProbeSetResultType));
		return this;
	};
} CompStatExpressionProbeSetResultType;

////////////////////////////////////////////////////////////////////

typedef struct {
	float x;
	float y;
} Coordinate;

////////////////////////////////////////////////////////////////////

typedef struct {
	Coordinate center;
	int	  numCell;
	float background;
	float noise;
} ZoneInfo;

////////////////////////////////////////////////////////////////////

typedef struct {
	int number_zones;
	float smooth_factor;
	ZoneInfo *pZones;
} AllZonesInfoType;

////////////////////////////////////////////////////////////////////

typedef struct {
	float avg;
	float stdv;
	float min;
	float max;
} AvgStdvMinMaxType;

////////////////////////////////////////////////////////////////////

typedef struct {
	affxcdf::GeneChipQCProbeSetType qcType;
	float avg;
	int count;
} ControlInformationType;
typedef std::list<ControlInformationType> ControlInformationList;

////////////////////////////////////////////////////////////////////

class CExpressionAlgorithmImplementation  
{
public:
	CExpressionAlgorithmImplementation();
	virtual ~CExpressionAlgorithmImplementation();
	void Clear();

	// Error
	string GetError() const;

	// Library path
	void SetLibPath(const char *libPath);

	// Write cel file
	void SetWriteCell(const bool writeCell);

	// Output directory
	void SetOutputDirectory(const char *outputDir);

	// Results
	int GetNumResults() const;
	bool DoesCompDataExists() const;
	string GetProbeSetName(int index) const;
	AbsStatExpressionProbeSetResultType *GetAbsStatResult(int index);
	CompStatExpressionProbeSetResultType *GetCompStatResult(int index);
	AbsStatExpressionProbeSetResultType *GetBaselineAbsStatResult(int index);

	// Parameters
	CExpStatAlgSettings & GetParameters();

	// Background zone information
	AllZonesInfoType & GetBackgroundZoneInfo();

	// The control information.
	ControlInformationList &GetControlInfo() { return m_ControlInfo; }

	// The raw noise value (RawQ)
	float GetRawQ() { return m_RawQ; }

	// The raw noise value (RawQ) for the baseline file.
	float GetBaselineRawQ() { return m_BaselineRawQ; }

	// The background stats
	AvgStdvMinMaxType &GetBgStats() { return m_BgStats; }

	// The noise stats
	AvgStdvMinMaxType &GetNoiseStats() { return m_NoiseStats; }

	// Algorithm
	bool RunStat(const char *celFile, const char *baseline, const char *cdfFile);

	// The number of rows of features
	int GetRows() { return m_Cdf.GetHeader().GetRows(); }

	// The number of cols of features
	int GetCols() { return m_Cdf.GetHeader().GetCols(); }

	// The chip type
	string GetChipType() { return m_Cdf.GetChipType(); }

	// The CEL object.
	FusionCELData &GetCellData() { return m_Cell; }

	// The baseline CEL object.
	FusionCELData &GetBaselineCellData() { return m_Baseline; }

	bool ReadCdfFile(FusionCELData& celData);

	/* For Average Intensity Measurement and Intensity Confidences Interval */
	void ComputeBackGroundZones(/*IN*/FusionCELData *pCell,
										/*IN/OUT*/ AllZonesInfoType & ZonesInfo,
										/*OUT*/vector<float> & featureIntensity);

private:
  FusionCELData m_Cell;
	FusionCELData m_Baseline;
  FusionCDFData m_Cdf;
	vector<float> m_tTable;
	map<string, int> m_ProbeSetNames;
	string m_Error;
	string m_LibPath;
	int m_NumResults;
	int m_NumFeatures;
	bool m_bCompDataExists;
	vector<AbsStatExpressionProbeSetResultType *> m_AbsStatResults;
	vector<AbsStatExpressionProbeSetResultType *> m_BaselineAbsStatResults;
	vector<CompStatExpressionProbeSetResultType *> m_CompStatResults;
	CExpStatAlgSettings m_Params;
	float m_RelNF;
	float m_RelNF2;
	float m_RawQ;
	float m_BaselineRawQ;
	map<int, bool> maskedProbes;

	// For storing the average/stdv/min/max background
	AvgStdvMinMaxType m_BgStats;

	// For storing the average/stdv/min/max noise
	AvgStdvMinMaxType m_NoiseStats;

	// For storing computed background zone information
	AllZonesInfoType m_ZonesInfo;

	// The list of control info.
	ControlInformationList m_ControlInfo;

	// Do we want to write a cel file
	bool m_WriteCell;

	// Output directory
	char* m_OutputDirectory;

	// Read the baseline CEL file.
	bool ReadBaselineCelFile(const char *baseline);

	// Allocate memory for the results
	void AllocateMemoryForResults();

	// Check the intensity data if all masked or blank.
	bool CheckIntensityData();

	// Compute the MAS5 stats.
	void ComputeStat();

	// Read the CEL file.
	bool ReadCelFile(const char *celFile);

	// Read the CDF file.
	bool ReadCdfFile(const char *cdfFile);

	/* For Average Intensity Measurement and Intensity Confidences Interval */
	void ComputeScaledAdjustedIntensity(/*IN*/FusionCELData *pCell,
										/*OUT*/vector<vector<float> > & PM,
										/*OUT*/vector<vector<float> > & MM,
										/*OUT*/vector<vector<bool> > & UseAtom,
										/*OUT*/vector<vector<FloatPair> > & BG,
										/*OUT*/vector<vector<FloatPair> > & Noise,
										/*IN/OUT*/ AllZonesInfoType & ZonesInfo,
										/*OUT*/vector<float> & featureIntensity);

	/* For Absolute Calls */
	void ComputeExpressionStat(FusionCELData *pCell, int iUnit, AbsStatExpressionProbeSetResultType *); 
	void ComputeAbsoluteCall(AbsStatExpressionProbeSetResultType *, ExpResults *);

	void ComputeMeasurement(/*IN */vector<vector<float> > & PM,
							/*IN */vector<vector<float> > & MM,
							/*IN */vector<vector<bool> > & UseAtom,
							/*OUT*/vector<float> & avgMeasurement,
							///*OUT*/vector<float> & intensityLow,
							///*OUT*/vector<float> & intensityHigh,
							/*OUT*/vector<vector<float> > & PV);
	float ModifyIntensitySlightly(/*IN */float intensity);
	float ComputeAdjustedIntensity(/*IN */float intensity, /*IN */float background, /*IN */float noise);
	void ComputeRanks(/*OUT*/vector<int> & rank, /*IN*/vector<float> & zoneI, /*IN*/int numCell);
	void SetMeasurement(/*IN*/vector<vector<FloatPair> > &BG,
						/*IN*/vector<float> &avgMeasurement,
						/*IN*/vector<AbsStatExpressionProbeSetResultType *> statResults);

	float DetermineScaleFactor(vector<AbsStatExpressionProbeSetResultType *> &statResults);
	void DetermineNormFactor(vector<AbsStatExpressionProbeSetResultType *> &expStatResults,
						     vector<AbsStatExpressionProbeSetResultType *> &baselineStatResults);
	bool UseUnitInNormFactor(string &probeSetName, list<string> &MaskedGenes);

	// Abs call functions.
	int DetermineZone(int cellx, int celly, int zonex, int zoney);
	void ModifyIntensities(vector<AbsStatExpressionProbeSetResultType *> statResults, float factor);

	bool AllMaskedOut(FusionCELData *pCell);
	bool IsBlankCellFile(FusionCELData *pCell);

	void GetUsedSet(/*IN */vector<vector<float> > & InputSet,
					/*IN */vector<vector<bool> > & UseAtom,
					/*OUT*/vector<vector<float> > & UsedSet);

	void ComputeContrastValue(/*IN */vector<vector<float> > & PM,
							  /*IN */vector<vector<float> > & MM,
							  /*OUT*/vector<vector<float> > & CT);
	void ComputeProbeValue(/*IN */vector<vector<float> > & PM,
		  				   /*IN */vector<vector<float> > & CT,
						   /*OUT*/vector<vector<float> > & PV);
	void ComputeTypicalDifference(/*IN */vector<vector<float> > & PM,
		  						  /*IN */vector<vector<float> > & MM,
								  /*OUT*/vector<float> & SB);

	// Comp call functions
	void ComputeFoldChange(vector<vector<float> > & PVb,
						 vector<vector<float> > & PVe,
						 vector<vector<bool> > & UseAtomB,
						 vector<vector<bool> > & UseAtomE);
	void ComputeProbeLogRatio
	(
		vector<vector<float> > & PVb,
		vector<vector<float> > & PVe,
		vector<vector<bool> > & UseAtomB,
		vector<vector<bool> > & UseAtomE,
		vector<vector<float> > & PLR
	);
	void ComputeAvgLogRatio
	(
		vector<vector<float> > & PLR,
		vector<float> & avgLogRatio
	);
	void BuildTTable(float level);
	void CallComparison
	(
		float sfe,
		float sfb,
		vector<AbsStatExpressionProbeSetResultType *> &baselineAbsStatResults,
		vector<vector<FloatPair> > & BGe,
		vector<vector<FloatPair> > & BGb
	);

	void ComputeIntensityDependentSignificances
	(
		float &gamma1,
		float &gamma2,
		vector<float> & PMe,
		vector<float> & PMb
	);

	void DetermineComparativeCall
	(
		CompStatExpressionProbeSetResultType *pCompData, 
		vector<float> &diffE,
		vector<float> &diffB,
		vector<float> &PMinusBgE,
		vector<float> &PMinusBgB,
		float gamma1, float gamma2
	);

	bool DetermineRelCallNormFactor(float sfe, float sfb);

	/* For calculating summary data */
	float ComputeCellBackground(int x, int y, ZoneInfo *zones, int zoneSize, float smoothFactorBg);

	void ComputeAvgMaxMinStdBgNoise(vector<vector<FloatPair> > & BG, vector<vector<FloatPair> > & Noise);

	void ReportCornerControls(affxcdf::GeneChipQCProbeSetType qcType);

	/* Determine if the probe is masked. */
	bool IsMasked(FusionCELData *cel, int x, int y);

	/* Read the contents of the probe MSK file. */
	bool ReadProbeMSKFile();

	/* Read the contents of the MSK file for normalization probe sets. */
	bool ReadNormMSKFile();

	/* Read the contents of the MSK file for scaling probe sets. */
	bool ReadScaleMSKFile();

};

////////////////////////////////////////////////////////////////////

#endif // !defined(AFX_EXPRESSIONALGORITHMIMPLEMENTATION_H__F70115DD_BD28_4226_9D03_61ACFD937B05__INCLUDED_)
