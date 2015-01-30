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
 * @file   snp-cluster-store-test.cpp
 * @author Chuck Sugnet
 * @date   Mon Mar 29 10:30:46 2010
 * 
 * @brief Test the functionality of SnpClusterStore and related
 * conversion utilities.
 */

#include "chipstream/SnpClusterStore.h"

#include "util/CalvinChpCheck.h"
#include "util/Convert.h"
#include "util/RegressionSuite.h"
#include "util/RegressionTest.h"
#include "util/Util.h"
#include "util/Verbose.h"
//
#include <cstring>
#include <string>
#include <vector>
//

using namespace std;

class SnpClusterStoreTest : public RegressionSuite {


public:

    int numPassed, numFailed;
    
    SnpClusterStoreTest() {
        numPassed = numFailed = 0;
    }

    void testBrlmmp();
    void testBirdseedPrior();
    void testBirdseedPosterior();
    void testBrlmmpA5();
    void checkVal(double a, double b);


};

void SnpClusterStoreTest::checkVal(double a, double b) {
    if (!Convert::doubleCloseEnough(a, b, 3)) {
        numFailed++;
    }
    else {
        numPassed++;
    }
}

void SnpClusterStoreTest::testBrlmmp() {
    SnpClusterStore store("../../../regression-data/data/idata/lib/Axiom_GW_Hu_SNP/Axiom_GW_Hu_SNP.r2.AxiomGT1.models");

    snp_distribution snp = store.getSnpCluster("AX-11086525", 2);
    checkVal(snp.bb.m, -3.2536);
    checkVal(snp.ab.m, -0.139578);
    checkVal(snp.aa.m, 2.70885);

    snp = store.getSnpCluster("AX-11088371", 1);
    checkVal(snp.bb.m, -1.98809);
    checkVal(snp.ab.m, -0.0154814);
    checkVal(snp.aa.m, 1.94164);
    
    if (store.snpClusterExists("AX-11088467", 1)) {
        numFailed++;
    }
    else {
        numPassed++;
    }
}

void SnpClusterStoreTest::testBirdseedPrior() {
    SnpClusterStore store("../../../regression-data/data/idata/lib/GenomeWideSNP_6/GenomeWideSNP_6.birdseed.models"); 

    snp_distribution snp = store.getSnpCluster("SNP_A-2012629", 2);
    // In this test the values for the aa and bb clusters are swapped relative
    // to the models file.  This is because the file has the labels swapped for
    // the aa and bb clusters, which is then fixed when the file is written to
    // db.
    checkVal(snp.bb.m, 0.1321);
    checkVal(snp.ab.m, 0.4194);
    checkVal(snp.aa.m, 0.713);

    snp = store.getSnpCluster("SNP_A-4267913", 1);
    checkVal(snp.bb.m, 0.2764);
    checkVal(snp.ab.m, 0);
    checkVal(snp.aa.m, 1.3194);
    
    if (store.snpClusterExists("SNP_A-2054822", 1)) {
        numFailed++;
    }
    else {
        numPassed++;
    }
}

void SnpClusterStoreTest::testBirdseedPosterior() {
    SnpClusterStore store("../../../regression-data/data/idata/lib/GenomeWideSNP_6/birdseed.snp-models.txt");

    snp_distribution snp = store.getSnpCluster("SNP_A-1984828", 2);
    checkVal(snp.bb.m, 424.846);
    checkVal(snp.ab.m, 709.288);
    checkVal(snp.aa.m, 951.04);

    snp = store.getSnpCluster("SNP_A-4267766", 1);
    checkVal(snp.bb.m, 177.538);
    checkVal(snp.ab.m, 0);
    checkVal(snp.aa.m, 423.182);
    
    if (store.snpClusterExists("SNP_A-4213095", 1)) {
        numFailed++;
    }
    else {
        numPassed++;
    }
}

void SnpClusterStoreTest::testBrlmmpA5() {
    SnpClusterStore store("../../../regression-data/data/idata/lib/GenomeWideSNP_6/brlmm-p.snp-posteriors.a5");

    snp_distribution snp = store.getSnpCluster("SNP_A-1984828", 2);
    checkVal(snp.bb.m, -2.38306);
    checkVal(snp.ab.m, -0.501296);
    checkVal(snp.aa.m, 1.40112);

    snp = store.getSnpCluster("SNP_A-4267766", 1);
    checkVal(snp.bb.m, -1.60326);
    checkVal(snp.ab.m, 0.0369297);
    checkVal(snp.aa.m, 1.69559);
    
    if (store.snpClusterExists("SNP_A-4213095", 1)) {
        numFailed++;
    }
    else {
        numPassed++;
    }
}

/** Everybody's favorite function. */
int main(int argc, char* argv[]) {
    try {
        SnpClusterStoreTest test;
        test.testBrlmmp();
        test.testBirdseedPrior();
        test.testBirdseedPosterior();
        test.testBrlmmpA5();
        Verbose::out(1, "NumPassed: " + ToStr(test.numPassed) + " NumFailed: " + ToStr(test.numFailed));
        return test.numFailed;
    }
    catch (...) {
        Verbose::out(1,"Unexpected Error: uncaught exception.");
        return 1;
    }
    return 1;
}
