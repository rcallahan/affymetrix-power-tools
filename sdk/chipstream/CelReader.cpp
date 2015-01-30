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
#include "chipstream/CelReader.h"
//
#include "chipstream/DiskIntensityMart.h"
#include "chipstream/IdxGroup.h"
//
#include "calvin_files/data/src/GenericData.h"
#include "calvin_files/fusion/src/FusionCELData.h"
#include "calvin_files/utils/src/StringUtils.h"
#include "util/Convert.h"
#include "util/Fs.h"
#include "util/Err.h"
#include "util/Util.h"
#include "util/Verbose.h"

using namespace affymetrix_fusion_io;
using namespace affymetrix_calvin_utilities;
using namespace affymetrix_calvin_exceptions;

/**
 * @brief Do the heavy lifting of reading data from the cel files
 * and passing to all the streams and intensity marts.
 * @return true on success or false on error.
 */
bool CelReader::readFiles() {
    int size = 0;
    int cel_channel_count = 0;
//   double newChip_times = 0, setProbeIntensity_times = 0, CelListener_newChip_times = 0;
//   clock_t clock1, clock2;
    /* sanity checks. */
    Err::check(m_FileNames.size() > 0, "CelReader::readFiles() - Can't specify 0 files to read.");
    Err::check(m_Streams.size() > 0 || m_IntenMarts.size() > 0 || m_CelListeners.size() > 0,
               "CelReader::readFiles() - Must specify a stream or IntensityMart to read into.");

    /* @todo: allow for more than one intensity mart?? */
    if (m_IntenMarts.size() > 1 && m_Streams.size() > 0) {
        Verbose::warn(0, "CelReader::readFiles -- more than one IntensityMart registered.  Using first one to populate ChipStreams.");
    }

    Verbose::progressBegin(1, "Reading and pre-processing " + ToStr(m_FileNames.size()) + " cel files", m_FileNames.size(), 0, m_FileNames.size());

    // %HACK% making extra marts for each chipstream.  first mart pointer will
    // be for raw intensities, and not passed on to a chipstream.
    if (m_IntenMarts.size() > 0 ) {
        int stream_mart_diff = m_Streams.size() - m_IntenMarts.size();
        if (stream_mart_diff >= 0) {
            for (int m = 0; m < stream_mart_diff + 1; m++) {
                IntensityMart* dm = m_IntenMarts[0]->copyMetaDataToEmptyMart();
                m_IntenMarts.push_back(dm);
            }
        }
    }

	

    std::vector<float> data;
    for (int fileIx=0; fileIx < m_FileNames.size(); fileIx++) {

    // process cel file
    FusionCELData cel;
    Verbose::progressStep(1);
    try {
      std::string tmp_unc_name=Fs::convertToUncPath(m_FileNames[fileIx]);
      cel.SetFileName(tmp_unc_name.c_str());
        if(!cel.Read())
            Err::errAbort("\nCan't read cel file: " + cel.GetFileName() +
                "\n>>> Error reported: " + StringUtils::ConvertWCSToMBS(cel.GetError()));

            // If this is a AGCC CEL file, then get channel names.  If
            // this is a GCOS CEL file, then an empty vector is returned.
            std::vector<std::wstring> data_channels = cel.GetChannels();

        for (int chanIx = 0; chanIx < data_channels.size(); chanIx++) {
          Verbose::out(4,"In cel file '"+m_FileNames[fileIx] +"' found channel '"+StringUtils::ConvertWCSToMBS(data_channels[chanIx])+"'");
        }

            // If this is a GCOS CEL file, then make up any old name for a
            // single channel -- subsequent call to SetActiveDataGroup
            // will be a no-op in this case, and GetIntensities will
            // proceed as normal.
            if (data_channels.empty()) {
                std::wstring temp = L"Default Group";
                data_channels.push_back(temp);
            }

            std::vector<unsigned int> channel_group;
            for (int chanIx = 0; chanIx < data_channels.size(); chanIx++) {
                cel.SetActiveDataGroup(data_channels[chanIx]);

                /* Check number of entries in cel file. */
                if (fileIx == 0)
                    size = cel.GetRows() * cel.GetCols();
                else
                    assert(size == cel.GetRows() * cel.GetCols());
                data.resize(size);
                cel.GetIntensities(0,data);

				// Compute saturation 
				int nsat = 0;
				for (int icel=0; icel< size; icel++)
				{
					if (data[icel] >= 3800)
						nsat++;
				}

				double sat = (double)nsat/(double)size;

				if (m_Saturation.size () < m_FileNames.size ()*data_channels.size())
					m_Saturation.resize ( m_FileNames.size ()*data_channels.size() );

				m_Saturation[data_channels.size()*fileIx+chanIx] = sat;

                /* Load data into our intensity marts. */
                //     clock1 = clock();
                for (int index = 0; index < m_IntenMarts.size(); index++) {

                    m_IntenMarts[index]->setProbeIntensity(data_channels.size()*fileIx+chanIx, data);
                }
                //     clock2 = clock();
                //     setProbeIntensity_times += clock2 - clock1;

                channel_group.push_back(cel_channel_count++);
            }
            /* Pass data to cel listeners */
            for (int index = 0; index < m_CelListeners.size(); index++) {
                m_CelListeners[index]->newChip(&cel);
            }
            m_CelChannels.addGroup("channels", channel_group);
            cel.Close();
        }
        catch(const Except &e) {
            Err::errAbort(ToStr("\n") + e.what());
        }
        catch(const std::bad_alloc &e) {
            Err::errAbort("\nRan out of memory when reading cel file: " + cel.GetFileName() + "\n");
        }
        catch(CalvinException &ce) {
            Err::errAbort("\nAGCC error when reading cel file: " + cel.GetFileName() + "\n"
                          "Description: '" +StringUtils::ConvertWCSToMBS(ce.Description()) + "'");
        }
        catch(const std::exception &e) {
            Err::errAbort("\nException caught. Message is: " + ToStr(e.what()));
        }
        catch (...) {
            Err::errAbort("\nUnknown problem when reading cel file: " + cel.GetFileName() + "\n");
        }
    } // end reading cel files

    Verbose::progressEnd(1, "Done.");

    std::string plural = (m_Streams.size() > 1 ? "s" : "");

    if (m_IntenMarts.size() > 0) {
        for (int martIx = 0; martIx < m_IntenMarts.size(); martIx++) {
            m_IntenMarts[martIx]->setChannelMapping(m_CelChannels);
        }

        if (m_Streams.size() > 0) {
            Verbose::out(1, "Processing " + ToStr(m_Streams.size()) + " chipstream" + plural + ".");
            for (int index = 0; index < m_Streams.size(); index++) {
                m_Streams[index]->newDataSet(m_IntenMarts[index+1]);
            }
        }
    }

    if (m_Streams.size() > 0)
        Verbose::out(1, "Finalizing " + ToStr(m_Streams.size()) + " chipstream"+plural+".");
    for (int index = 0; index < m_Streams.size(); index++) {
        if (m_Streams[index] != NULL)
            m_Streams[index]->endDataSet();
    }

    return true;
}
