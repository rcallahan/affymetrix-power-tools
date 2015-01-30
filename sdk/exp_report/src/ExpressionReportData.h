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

#ifndef _ExpressionReportData_HEADER_
#define _ExpressionReportData_HEADER_

/*! \file ExpressionReportData.h This file provides storage for expression reports.
 */

#include "exp_report/src/ChangeStats.h"
#include "exp_report/src/DetectionStats.h"
#include "exp_report/src/ExpressionControlResult.h"
#include "exp_report/src/ProbeSetStats.h"
//
#include <cstring>
#include <fstream>
#include <list>
#include <string>
//

namespace ExpressionReport
{

/*! A type to hold a name value pair. */
typedef struct _NameValuePair
{
	/*! The name of the parameter. */
	std::string name;

	/*! The value of the parameter. */
	std::string value;

	/*! Clears the data. */
	void  Clear() { name = ""; value = ""; }

} NameValuePair;

/*! An STL list of name value pairs. */
typedef std::list<NameValuePair> NameValuePairList;

/*! A type to hold a name float value pair. */
typedef struct _NameFloatValuePair
{
	/*! The name of the parameter. */
	std::string name;

	/*! The value of the parameter. */
	float value;

	/*! Clears the data. */
	void  Clear() { name = ""; value = 0.0f; }

} NameFloatValuePair;

/*! An STL list of name float value pairs. */
typedef std::list<NameFloatValuePair> NameFloatValuePairList;

/*! A type to hold an average, std, min and max values. */
typedef struct _AvgStdvMinMax
{
	/*! The average value. */
	float avg;

	/*! The stdv value. */
	float std;

	/*! The minimum value. */
	float min;

	/*! The maximum value. */
	float max;

	/*! Clears the data. */
	void Clear() { avg=0.0f; std=0.0f; min=0.0f; max=0.0f; }

} AvgStdvMinMax;

/*! A type to hold a name, average and count. */
typedef struct _NameAvgCount
{
	/*! The name. */
	std::string name;

	/*! The average value. */
	float avg;

	/*! The count. */
	int count;

	/*! Clears the data. */
	void  Clear() { name = ""; avg=0.0f; count=0; }

} NameAvgCount;

/*! An stl list of name, avg, count objects. */
typedef std::list<NameAvgCount> NameAvgCountList;

/*! Stores the contents of a Expression report. */
class ExpressionReportData
{
public:
	/*! Constructor */
	ExpressionReportData() { Clear(); }

	/*! Destructor */
	~ExpressionReportData() { Clear(); }

protected:

	/*! The array type */
	std::string arrayType;

	/*! The algorithm name */
	std::string algName;

	/*! The name of the CHP file that was used to compute the report. */
	std::string chpFileName;

	/*! A list of algorithm parameters. */
	NameValuePairList algParams;

	/*! The date */
	std::string date;

	/*! The probe pair threshold. */
	int probePairThreshold;

	/*! Flag to indicate if anti-sense probe sets should be reported. */
	bool antiSenseControls;

	/*! The background values. */
	AvgStdvMinMax backgroundStats;

	/*! The noise values. */
	AvgStdvMinMax noiseStats;

	/*! The list of control stats .*/
	NameAvgCountList controlStats;

	/*! The probe set results. */
	ProbeSetStats setStats;

	/*! The increase call stats. */
	ChangeStats increaseStats;

	/*! The decrease call stats. */
	ChangeStats decreaseStats;

	/*! The no change call stats. */
	ChangeStats noChangeStats;

	/*! The spike control results. */
	ExpressionControlResultList spikeStats;

	/*! The housekeeping control results. */
	ExpressionControlResultList housekeepingStats;

    /*! A list of probe set values (name/quantification). */
    NameFloatValuePairList probeSetValues;

public:

	/*! The probe array type. */
	std::string &ArrayType() { return arrayType; }

	/*! The algorithm name */
	std::string &AlgName() { return algName; }

	/*! The name of the CHP file that was used to compute the report. */
	std::string &CHPFileName() { return chpFileName; }

	/*! A list of algorithm parameters. */
	NameValuePairList &AlgParams() { return algParams; }

	/*! The date of the report */
	std::string &Date() { return date; }

	/*! The background values. */
	AvgStdvMinMax &BackgroundStats() { return backgroundStats; }

	/*! The noise values. */
	AvgStdvMinMax &NoiseStats() { return noiseStats; }

	/*! The list of control stats .*/
	NameAvgCountList &ControlStats() { return controlStats; }

	/*! The report results. */
	ProbeSetStats &ProbeSetResults() { return setStats; }

	/*! The increase call stats. */
	ChangeStats &IncreaseStats() { return increaseStats; }

	/*! The decrease call stats. */
	ChangeStats &DecreaseStats() { return decreaseStats; }

	/*! The no change call stats. */
	ChangeStats &NoChangeStats() { return noChangeStats; }

	/*! The spike control results. */
	ExpressionControlResultList &SpikeStats() { return spikeStats; }

	/*! The housekeeping control results. */
	ExpressionControlResultList &HousekeepingStats() { return housekeepingStats; }

	/*! The probe pair threshold. */
	int &ProbePairThreshold() { return probePairThreshold; }

	/*! Flag to indicate if anti-sense probe sets should be reported. */
	bool &AntiSenseControls() { return antiSenseControls; }

    /*! A list of probe set values (name/quantification). */
    NameFloatValuePairList &ProbeSetValues() { return probeSetValues; }

	/*! Clears memory associated with the class */
	void Clear()
	{
		arrayType = "";
		algName = "";
		chpFileName = "";
		algParams.clear();
		date = "";
		probePairThreshold = 0;
		antiSenseControls = true;
		backgroundStats.Clear();
		noiseStats.Clear();
		controlStats.clear();
		setStats.Clear();
		increaseStats.Clear();
		decreaseStats.Clear();
		noChangeStats.Clear();
		spikeStats.clear();
		housekeepingStats.clear();
        probeSetValues.clear();
	}
};

}

#endif
