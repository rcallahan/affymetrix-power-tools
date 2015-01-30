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
 * @file   PmAdjusterTest.cpp
 * @author csugne
 * @date   Mon Nov  7 10:04:34 PST 2005
 * 
 * @brief  Testing the PmAdjuster functions.
 * 
 */
#ifndef PMADJUSTERTEST_H
#define PMADJUSTERTEST_H

//
#include "chipstream/ChipLayout.h"
#include "chipstream/GcAdjust.h"
#include "chipstream/MmAdjust.h"
#include "chipstream/PmAdjuster.h"
#include "chipstream/PmOnlyAdjust.h"
#include "chipstream/SparseMart.h"
#include "chipstream/TableReader.h"
#include "util/TableFile.h"
//
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
//
#include <string>
#include <vector>

using namespace std;
/**
 * @class PmAdjusterTest
 * @brief cppunit class for testing conversion functions.
 */
class PmAdjusterTest : public CppUnit::TestFixture {
public:
  PmAdjusterTest() {}

  CPPUNIT_TEST_SUITE( PmAdjusterTest );
  CPPUNIT_TEST( pmOnlyAdjustTest );
  CPPUNIT_TEST( mmAdjustTest );
  CPPUNIT_TEST( gcAdjustTest );
  CPPUNIT_TEST_SUITE_END();

  /// test pm only adjustment
  void pmOnlyAdjustTest();
  /// test mismatch adjustment
  void mmAdjustTest();
  /// test gc median adjustment
  void gcAdjustTest();
  vector<ChipStream *> m_ITrans;
};

// Registers the fixture into the registry
CPPUNIT_TEST_SUITE_REGISTRATION( PmAdjusterTest );

void PmAdjusterTest::gcAdjustTest() {
  float goldAdjustment[] = {17.0, 23.0, 34.5}; // median calculated in R.
  /* Read pgf file. */
  GcAdjust gc;
  ChipLayout layout;
  layout.setNeedGc(true);
  std::set<const char *, Util::ltstr> probeSetsToLoad;
  vector<bool> probeSubset;
  layout.openPgf("input/mini.pgf", 10, 10, probeSetsToLoad, probeSubset,"");
  vector<Probe *> gcControlProbes;
  /* Read in data. */
  vector<bool> probes(100);
  for(int i = 0; i < probes.size(); i++) 
    probes[i] = true;
  SparseMart iMart(probes);

  /* Insert a number of probes... */
  ProbeListPacked pl = layout.getProbeListAtIndex(0);
  ProbeSet *ps = ProbeListFactory::asProbeSet(pl);
  for(int atomIx = 0; atomIx < ps->atoms.size(); atomIx++) {
    Atom &a = *(ps->atoms[atomIx]);
    for(int probeIx = 0; probeIx < a.probes.size(); probeIx++) {
      Probe *p = a.probes[probeIx];
      gcControlProbes.push_back(p);
    }
  }
  gc.setLayout(layout, gcControlProbes);

  /* Read in as IntensityMart. */
  vector<string> fileNames;
  TableReader reader(100);  
  reader.registerIntensityMart(&iMart);
  fileNames.push_back("input/norm-data.txt");
  reader.setFiles(fileNames);
  reader.readFiles();
  /* We're going to have a lot of empty bins, avoid the error messages. */
  Verbose::setLevel(0);
  for(int i = 0; i < iMart.getCelDataSetCount(); i++) {
    assert(i < sizeof(goldAdjustment)/sizeof(goldAdjustment[0]));
    float adjustment = 0;
    float orig = iMart.getProbeIntensity(0,0);
    gc.pmAdjustment(0,i,iMart,m_ITrans, orig, adjustment);
    CPPUNIT_ASSERT( adjustment == goldAdjustment[i] );
  }
  delete ps;
  Verbose::setLevel(1);
}

void PmAdjusterTest::mmAdjustTest() {
  int psIx = 0;
  std::set<const char *, Util::ltstr> probeSetsToLoad;
  vector<bool> probeSubset;
  /* Read pgf file. */
  ChipLayout layout;
  layout.setNeedMismatch(true);
  layout.openPgf("input/mini.pgf", 10, 10, probeSetsToLoad, probeSubset, "");
  MmAdjust mm; 
  mm.setLayout(layout);
  /* Read in data. */
  vector<bool> probes(100);
  for(int i = 0; i < probes.size(); i++) 
    probes[i] = true;
  SparseMart iMart(probes);

  /* Read in data for answer */
  TableFile toNorm('\t','#',false,false);
  toNorm.open("input/norm-data.txt");

  /* Read in as IntensityMart. */
  vector<string> fileNames;
  TableReader reader(100);  
  reader.registerIntensityMart(&iMart);
  fileNames.push_back("input/norm-data.txt");
  reader.setFiles(fileNames);
  reader.readFiles();

  /* Check each probe set, atom, probe... */
  for(psIx = 0; psIx < layout.getProbeSetCount(); psIx++) {
    ProbeListPacked pl = layout.getProbeListAtIndex(0);
    ProbeSet *ps = ProbeListFactory::asProbeSet(pl);
    //    const ProbeSet &ps = layout.getProbeSetIndex(psIx);
    for(int atomIx = 0; atomIx < ps->atoms.size(); atomIx++) {
      Atom &a = *(ps->atoms[atomIx]);
      for(int probeIx = 0; probeIx < a.probes.size(); probeIx++) {
        Probe &p = *(a.probes[probeIx]);
        /* Only look at probes that are PM. */
        if(p.type == Probe::PMST || p.type == Probe::PMAT) {
          for(int chipIx = 0; chipIx < iMart.getCelDataSetCount(); chipIx++) {
            float expected = iMart.getProbeIntensity(p.id+1, chipIx);
            float calculated = 0;
            float orig = iMart.getProbeIntensity(p.id, chipIx);
            mm.pmAdjustment(p.id, chipIx, iMart, m_ITrans, orig, calculated);
            CPPUNIT_ASSERT( expected == calculated );
          }
        }
      }
    }
    delete ps;
  }
}

void PmAdjusterTest::pmOnlyAdjustTest() {
  PmOnlyAdjust pmOnly;
  vector<bool> probes;
  SparseMart iMart(probes);
  CPPUNIT_ASSERT( pmOnly.getType() == "pm-only");
  float adjustment = 0;
  float orig = 100;
  pmOnly.pmAdjustment(0, 0, iMart, m_ITrans, orig, adjustment);
  CPPUNIT_ASSERT( adjustment == 0 );
  orig = 10;
  pmOnly.pmAdjustment(0, 0, iMart, m_ITrans, orig, adjustment);
  CPPUNIT_ASSERT( adjustment == 0 );
  orig = 1;
  pmOnly.pmAdjustment(0, 0, iMart, m_ITrans, orig, adjustment);
  CPPUNIT_ASSERT( adjustment == 0 );
}

#endif 
