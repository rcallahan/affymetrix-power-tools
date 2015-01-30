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

/*
 * FILE ClustersReader.h
 */

#ifndef _CLUSTERSREADER_H
#define _CLUSTERSREADER_H

//
#include "birdseed-dev/Clusters.h"
#include "birdseed-dev/PriorsReader.h"
//
#include <map>

namespace birdseed 
{
namespace dev
{
    class TextClustersReader
    {
      private:
        typedef std::map<std::string, Clusters> ClustersMap;
        ClustersMap clustersMap;

      public:
      TextClustersReader(const std::string& textPath);
        // Please don't quibble about the fact that this method is not named correctly.
        // It makes PriorsReaderTemplate work properly.
      const Clusters *getPrior(const std::string& clusterName);
    };


    typedef PriorsReaderTemplate<TextClustersReader, Clusters> ClustersReader;
};
};


#endif /* _CLUSTERSREADER_H */

/******************************************************************/
/**************************[END OF ClustersReader.h]**********************/
/******************************************************************/
/* Emacs configuration
 * Local Variables:
 * mode: C++
 * tab-width:4
 * End:
 */
