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
#include "rawq/src/RawQ.h"
//
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <list>
#include <string.h>
#include <string>
#include <vector>
//

using namespace affymetrix_fusion_io;
using namespace std;

#define EXPRESSION_ATOMS_PER_CELL 2
#define DEFAULT_VERT_ZONES 4
#define DEFAULT_HORZ_ZONES 4
#define DEFAULT_PERCENT_BG 2

#define LARGE_FLOAT_NUMBER 999999999999999999.0f

/////////////////////////////////////////////////////////////////////////////

CRawQ::CRawQ()
{
	SetDefaults();
}

/////////////////////////////////////////////////////////////////////////////

CRawQ::~CRawQ()
{
}

/////////////////////////////////////////////////////////////////////////////

void CRawQ::SetDefaults()
{
	m_VertZones = DEFAULT_VERT_ZONES;
	m_HorzZones = DEFAULT_HORZ_ZONES;
	m_PercentBGCells = DEFAULT_PERCENT_BG;
}

/////////////////////////////////////////////////////////////////////////////

int DetermineZone(int x, int y, int zonex, int zoney, int vZones)
{
	float fZx = (float) x / (float) zonex;
	float fZy = (float) y / (float) zoney;
	int Zx = (int)floor(fZx);
	int Zy = (int)floor(fZy);
	int zoneID = Zx + Zy * vZones;
	return zoneID;
}

/////////////////////////////////////////////////////////////////////////////

float CRawQ::ComputeRawQ(FusionCELData &cell, FusionCDFData &cdf)
{
	typedef struct {
		float intensity;
		float stdev;
		int   pixel;
	} CellStatisticsType;


	// Store the number of background cells.
	FusionCDFFileHeader &cdfHeader = cdf.GetHeader();
	float numer = cdfHeader.GetCols() * cdfHeader.GetRows() * m_PercentBGCells;
	float denom = m_VertZones * m_HorzZones * 100.0f;
	int bgCells = (int) (numer / denom);


	// Determine the number of remaining cells in the vertical direction.
	int CellsRemaining = cdfHeader.GetCols() % m_VertZones;
	int zonex;
	int zoney;
	if(CellsRemaining != 0)
		zonex = (cdfHeader.GetCols() + 
				m_VertZones - CellsRemaining) /
				m_VertZones;
	else
		zonex = cdfHeader.GetCols() / m_VertZones;

	// Determine the number of remaining cells in the horizontal direction.
	CellsRemaining = cdfHeader.GetRows() % m_HorzZones;
	if(CellsRemaining != 0)
		zoney = (cdfHeader.GetRows() +
				m_HorzZones-CellsRemaining) /
				m_HorzZones;
	else
		zoney = cdfHeader.GetRows() / m_HorzZones;

	// Ensure that there are a match and mismatch cell in the same zone.
	zoney += zoney % EXPRESSION_ATOMS_PER_CELL;


	// Determine the total number of zones.
	int NumberZones = m_VertZones * m_HorzZones;

	// Allocate memory for storing background data.
	float *zonebg = new float[NumberZones];
	float *zonenoise = new float[NumberZones];
	CellStatisticsType *bgN = new CellStatisticsType[NumberZones*bgCells];

	// Get number of units.
	int NumUnits = cdfHeader.GetNumProbeSets();

	// Determine the total number of atoms in the chip.
	FusionCDFProbeSetInformation unit;
	int totalCells=0;
	for (int iUnit=0; iUnit<NumUnits; ++iUnit)
	{
		if (cdf.GetProbeSetType(iUnit) == affxcdf::ExpressionProbeSetType)
		{
			cdf.GetProbeSetInformation(iUnit, unit);
			totalCells += unit.GetNumCells();
		}
	}

	// Allocate space for all atoms intensities and ID's.
	int *entryIndex = new int[totalCells];
	int *intenszid = new int[totalCells];

	// Clear arrays.
	memset(intenszid, 0, sizeof(int)*totalCells);
	memset(entryIndex, 0, sizeof(int)*totalCells);

	// Loop over all units to determine the zone ID's and intensities.
	int iInten=0;
	for (int iUnit=0; iUnit<NumUnits; ++iUnit)
	{
		// Only process expression units.
		if (cdf.GetProbeSetType(iUnit) != affxcdf::ExpressionProbeSetType)
			continue;

		// Get the PM and MM intensity objects and their zone ids.
		cdf.GetProbeSetInformation(iUnit, unit);
		FusionCDFProbeGroupInformation group;
		int nGroups = unit.GetNumGroups();
		for (int iGroup=0; iGroup<nGroups; iGroup++)
		{
			unit.GetGroupInformation(iGroup, group);
			int nCells = group.GetNumCells();
			FusionCDFProbeInformation probe;
			for (int iCell=0; iCell<nCells; iCell++)
			{
				group.GetCell(iCell, probe);
				entryIndex[iInten] = cell.XYToIndex(probe.GetX(), probe.GetY());
				intenszid[iInten] = DetermineZone(probe.GetX(), probe.GetY(), zonex, zoney, m_VertZones);
				++iInten;
			}
		}
	}

	// compute background for each zone
	for (int iZone=0; iZone<NumberZones; iZone++)
	{
		// Initialize the background.
		for(int bgcnt = 0; bgcnt < bgCells; bgcnt++)
			bgN[bgcnt+(iZone*bgCells)].intensity = LARGE_FLOAT_NUMBER;

		// find the lowest N intensities in each zone
		for (int iInten = 0; iInten < totalCells; iInten++)
		{
			// Only process those intensities in the current zone.
			if(intenszid[iInten] == iZone)
			{
				int index_cnt;
				int index_max;
				for (index_cnt=1, index_max=0;
					 index_cnt < bgCells; index_cnt++)
				{
					if(bgN[index_cnt+(iZone*bgCells)].intensity > bgN[index_max+(iZone*bgCells)].intensity)
						index_max = index_cnt;
				}


				// Store the low intensity.
				float intensity = min(bgN[index_max+(iZone*bgCells)].intensity, cell.GetIntensity(entryIndex[iInten]));
				if (intensity != bgN[index_max+(iZone*bgCells)].intensity)
				{
					bgN[index_max+(iZone*bgCells)].intensity = cell.GetIntensity(entryIndex[iInten]);
					bgN[index_max+(iZone*bgCells)].pixel = cell.GetPixels(entryIndex[iInten]);
					bgN[index_max+(iZone*bgCells)].stdev = cell.GetStdv(entryIndex[iInten]);
				}
			}
		}

		// compute the average
		float bgSum = 0.0f;
		for(int bgcnt = 0; bgcnt < bgCells; bgcnt++)
			bgSum += bgN[bgcnt+(iZone*bgCells)].intensity;
		zonebg[iZone] = bgSum / bgCells;

		// Compute the noise.
		bgSum = 0.0f;
		for(int bgcnt = 0; bgcnt < bgCells; bgcnt++)
			bgSum += bgN[bgcnt+(iZone*bgCells)].stdev / (float)(sqrt((float)bgN[bgcnt+(iZone*bgCells)].pixel));
		zonenoise[iZone] = bgSum / bgCells;
	}

	// Compute the average noise.
	float avgNoise=0.0f;
	int totalZonesUsed = 0;
	for (int iZone=0; iZone<NumberZones; iZone++)
	{
		if (zonenoise[iZone] >= 0)
		{
			avgNoise+=zonenoise[iZone];
			++totalZonesUsed;
		}
	}
	float rawQ = avgNoise/totalZonesUsed;

	//Clean up
	delete [] entryIndex;
	delete [] intenszid;
	delete [] bgN;
	delete [] zonebg;
	delete [] zonenoise;

	// Return the RawQ value
 	return rawQ;
}

/////////////////////////////////////////////////////////////////////////////
