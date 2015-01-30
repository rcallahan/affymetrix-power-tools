////////////////////////////////////////////////////////////////
//
// Copyright (C) 2007 The Broad Institute and Affymetrix, Inc.
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
#include "birdseed-dev/ClustersReader.h"
//
#include "broadutil/BroadUtil.h"
#include "util/Verbose.h"
//
#include <cstring>
#include <string.h>
#include <string>
//

using namespace std;
using namespace birdseed::dev;

static void setCluster(Clusters *clusters, const char *strCluster, size_t i)
{
    double covar;
    int ret = sscanf(strCluster, "%lf %lf %lf %lf %lf %lf",
                     &clusters->means[i][0], &clusters->means[i][1], &clusters->vars[i][0], &covar, &clusters->vars[i][1], &clusters->weights[i]);
    if (ret != 6) {
        throw BroadException("Error parsing cluster", __FILE__, __LINE__, strCluster);
    }
    clusters->covar = covar / sqrt(clusters->vars[i][0] * clusters->vars[i][1]);
}

/**
 * Because clusters get filled in backward, if a haploid set is read, the two clusters
 * must be slid down.
 */
static void slideClustersDown(Clusters *clusters)
{
    assert(clusters->k == 2);
    for (size_t i = 0; i < 2; ++i) {
        clusters->means[i] = clusters->means[i+1];
        clusters->vars[i] = clusters->vars[i+1];
        clusters->weights[i] = clusters->weights[i+1];
    }
    clusters->means.resize(2);
    clusters->vars.resize(2);
    clusters->weights.resize(2);
}

TextClustersReader::TextClustersReader(const std::string& textPath)
{
  FILE *fp = fopen_check(textPath.c_str(), "r");
    char buf[10000];
    while (fgets(buf, sizeof(buf), fp) != NULL) {
        if (strlen(buf) >= sizeof(buf) - 1) {
          throw BroadException("Line too long in text clusters file.", __FILE__, __LINE__, textPath.c_str());
        }
        const char *name = strtok(buf, ";");
        if (name == NULL) {
          throw BroadException("Error parsing text clusters file.",  __FILE__, __LINE__, textPath.c_str());
        }
        Clusters clusters(1, MAX_NUM_CLUSTERS, 0);
        // Finny produces priors in the reverse order of what Josh expects.
        for (int i = MAX_NUM_CLUSTERS - 1; i >= 0; --i) {
            const char *clusterString = strtok(NULL, ";");
            if (clusterString == NULL) {
                if (i == 0) {
                    clusters.k = 2;
                    slideClustersDown(&clusters);
                    break;
                }
                throw BroadException("Error parsing text clusters file.",  __FILE__, __LINE__, textPath.c_str());
            }
            
            setCluster(&clusters, clusterString, i);
        }
        if (!clustersMap.insert(ClustersMap::value_type(name, clusters)).second) {
            throw BroadException("Cluster seen more than once text clusters file.",  __FILE__, __LINE__, name);
        }
    }
    if (ferror(fp)) {
      throw BroadException("Error reading text clusters file.", __FILE__, __LINE__, textPath.c_str(), errno);
    }
    fclose(fp);
}

const Clusters *TextClustersReader::getPrior(const std::string& clusterName)
{
  ClustersMap::const_iterator it = clustersMap.find(clusterName);
    if (it == clustersMap.end()) {
        Verbose::warn(2, std::string("Could not find cluster in text clusters file: ") + clusterName);
        return NULL;
    }
    return &(*it).second;
}


/******************************************************************/
/**************************[END OF ClustersReader.cpp]*************************/
/******************************************************************/

/* Emacs configuration
 * Local Variables:
 * mode: C++
 * tab-width:4
 * End:
 */
