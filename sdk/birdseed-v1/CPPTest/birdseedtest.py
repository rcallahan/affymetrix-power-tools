#!/usr/bin/env python
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2006 by the
# Broad Institute/Massachusetts Institute of Technology. All rights are
# reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for its
# use, misuse, or functionality.

from __future__ import division
import filecmp
import os
import shutil
import subprocess
import sys
import unittest

def compare_confidences(strFile1, strFile2, fDelta=0.00001):
    f1 = open(strFile1)
    f2 = open(strFile2)
    bRet = True
    strHeader1 = f1.readline()
    strHeader2 = f2.readline()
    if strHeader1 != strHeader2:
        print >> sys.stderr, "Header line mismatch in", strFile1, strFile2
        bRet = False
    for iLineNum, strLine1 in enumerate(f1):
        strLine2 = f2.readline()
        lstFields1 = [float(strField) for strField in strLine1.split()[1:]]
        lstFields2 = [float(strField) for strField in strLine2.split()[1:]]
        for i in xrange(len(lstFields1)):
            if abs(lstFields1[i] - lstFields2[i]) > fDelta:
                print >> sys.stderr, "Confidence mismatch at line %d in files %s and %s: %f %f; delta %f" % \
                (iLineNum, strFile1, strFile2, lstFields1[i], lstFields2[i], abs(lstFields1[i] - lstFields2[i]))
                bRet = False
    f1.close()
    f2.close()
    return bRet

def compare_clusters(strFile1, strFile2, fDelta=0.1001):
    """There is a little jitter in the clusters when reading them and then writing them again in the
    forced-cluster case, which is why the delta is so large."""
    f1 = open(strFile1)
    f2 = open(strFile2)
    bRet = True
    for iLineNum, strLine1 in enumerate(f1):
        strLine2 = f2.readline()
        strClusterName1 = strLine1.split(None, 1)[0]
        strClusterName2 = strLine2.split(None, 1)[0]
        if strClusterName1 != strClusterName2:
            print >> sys.stderr, "Cluster name mismatch at line %d in files %s and %s: %s %s" % \
                  (iLineNum, strFile1, strFile2, strClusterName1, strClusterName2)
        lstFields1 = [float(strField) for strFields in strLine1.split()[1:] for strField in strFields.split(";")]
        lstFields2 = [float(strField) for strFields in strLine2.split()[1:] for strField in strFields.split(";")]
        for i in xrange(len(lstFields1)):
            if abs(lstFields1[i] - lstFields2[i]) > fDelta:
                print >> sys.stderr, "Cluster mismatch at line %d in files %s and %s: %f %f; delta %f" % \
                (iLineNum, strFile1, strFile2, lstFields1[i], lstFields2[i], abs(lstFields1[i] - lstFields2[i]))
                bRet = False
    f1.close()
    f2.close()
    return bRet


class BirdseedTest(unittest.TestCase):
    def setUp(self):
        thisDir = os.path.dirname(__file__) + "/"
        self.outputDir = thisDir + "output/"
        if not os.path.isdir(self.outputDir):
            os.mkdir(self.outputDir)
        self.inputDir = thisDir + "input/"
        self.expectedDir = thisDir + "expected/"
        self.executable = thisDir + "../Debug/birdseed"

    def tearDown(self):
        shutil.rmtree(self.outputDir)
        pass

    def test_noGender(self):
        self._test_and_compare("nogender", gender=False)

    def test_GenderMonomorphic(self):
        self._test_and_compare("gender_monomorphic")
        
    def test_Gender2Clusters(self):
        self._test_and_compare("gender_2clusters")
        
    def test_Gender2ClustersSpecialSnps(self):
        self._test_and_compare("gender_2clusters", bSpecialSnps=True)
        
    def test_noGenderForceClusters(self):
        self._test_and_compare("nogender", gender=False, bForceClusters=True)

    def test_GenderMonomorphicForceClusters(self):
        self._test_and_compare("gender_monomorphic", bForceClusters=True)
        
    def test_Gender2ClustersForceClusters(self):
        self._test_and_compare("gender_2clusters", bForceClusters=True)

    def xtest_notEnoughSamples(self):
        """This test is disabled because we now allow a single sample."""
        self._test_failure("singlesample", "Not enough samples to clusters.  6 are needed but there are only 1.",
                           gender=False)

    def _test_and_compare(self, strTestName, correctionFactor=None, gender=True, bForceClusters=False, bSpecialSnps=False):
        lstArgs = [self.executable,
                   "--write-clusters", self.outputDir + strTestName + ".clusters"]

        if bSpecialSnps:
            lstArgs += ["--special-snps", self.inputDir + "BI_SNP.special_snps"]
        else:
            lstArgs += ["--chrX-snps", self.inputDir + "BI_SNP.chrx"]

        if bForceClusters:
            lstArgs += ["-clusters", self.expectedDir + strTestName + ".clusters"]
        else:
            lstArgs += ["--priors-text", self.inputDir + "priors.txt"]

        if correctionFactor is not None:
            lstArgs += ["-c", str(correctionFactor)]

        if gender:
            lstArgs += ["--gender-file", self.inputDir + strTestName + ".gender"]

        lstArgs += [self.inputDir + strTestName + ".intensities",
                    self.outputDir + strTestName + ".calls",
                    self.outputDir + strTestName + ".confidences"]

        print " ".join(lstArgs)
        proc = subprocess.Popen(lstArgs)

        self.assertEqual(proc.wait(), 0, "Non-zero exit status from birdseed")

        for strFileType in [".calls"]:
            strFileName = strTestName + strFileType
            self.assert_(filecmp.cmp(self.expectedDir + strFileName, self.outputDir + strFileName),
                         "Mismatch on file " + strFileName)

        strClustersFile = strTestName + ".clusters"
        self.assert_(compare_clusters(self.expectedDir + strClustersFile, self.outputDir + strClustersFile),
                     "Mismatch on file " + strClustersFile)
        
        strConfidencesFile = strTestName + ".confidences"
        self.assert_(compare_confidences(self.expectedDir + strConfidencesFile, self.outputDir + strConfidencesFile),
                     "Mismatch on file " + strConfidencesFile)
                   
    def _test_failure(self, strTestName, strExpectedErrorMessage, correctionFactor=None, gender=True, bForceClusters=False):
        lstArgs = [self.executable,
                   "-v", "1",
                   "--chrX-snps", self.inputDir + "BI_SNP.chrx"]
        if bForceClusters:
            lstArgs += ["-clusters", self.expectedDir + strTestName + ".clusters"]
        else:
            lstArgs += ["--priors-text", self.inputDir + "priors.txt"]
        if correctionFactor is not None:
            lstArgs += ["-c", str(correctionFactor)]
        if gender:
            lstArgs += ["--gender-file", self.inputDir + strTestName + ".gender"]

        lstArgs += [self.inputDir + strTestName + ".intensities",
                    self.outputDir + strTestName + ".calls",
                    self.outputDir + strTestName + ".confidences"]

        print " ".join(lstArgs)
        proc = subprocess.Popen(lstArgs, stderr=subprocess.PIPE)
        strError = proc.stderr.read()
        self.assertNotEqual(proc.wait(), 0, "Expected non-zero exit status from birdseed")
        self.assertEqual(strError.count(strExpectedErrorMessage), 1,
                         "Expected to see '" + strExpectedErrorMessage + "' in '" +\
                         strError + "'")
        
        
                   
        
if __name__ == '__main__':
    unittest.main()
