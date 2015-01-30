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
 * @file   ChipStream.cpp
 * @author Chuck Sugnet
 * @date   Fri Oct 21 18:04:46 2005
 * 
 * @brief Base class for objects that can initialize themselves one at
 * a time and transform intensity data.
 */

//
#include "chipstream/ChipStream.h"
//
#include <cstdio>
#include <vector>
//

// constructor
ChipStream::ChipStream() : m_ParentStream(NULL), m_TransformedIMart(NULL) {
  /*     setupSelfDoc(*this); */
}

/** 
 * @brief Virtual destructor for a virtual class.
 */
ChipStream::~ChipStream() {
  if (m_Streams.size() == 0) {
    delete m_TransformedIMart;
  }
}

/** 
 * Fill in the information for Self documentation.
 * @param doc - Self documenter to be filled in.
 */
void ChipStream::setupSelfDoc(SelfDoc &doc) {
  doc.setDocName("no-trans");
  doc.setDocDescription("Placeholder chipstream that does no transformation");
  doc.setDocOptions(getDefaultDocOptions());
}

/** 
 * @brief Supply a little how/what/why about the algorithms this
 * class performs and what parameters it takes.
 * @return SelfDoc
 */
SelfDoc ChipStream::explainSelf() { 
  SelfDoc doc;
  setupSelfDoc(doc);
  return doc;
}

/** 
 * @brief This static function should be overridden by child classes
 * to return an object of the correct type initialized correctly
 * with the parameters in the string, string map. All objects
 * created this way should be deleted when finished using.
 * 
 * @param param - Map of key/value pairs to initialize the object.
 * 
 * @return Pointer toCreate object, this should be sub casted as necessary.
 */
SelfCreate *ChipStream::newObject(std::map<std::string,std::string> &param) {
  SelfDoc doc = explainSelf();
  ChipStream *cs = new ChipStream();
  return cs;
}

/** 
 * @Brief Add a stream to the list of those that will receive
 * downstream data.
 * @param stream - ChipStream that wants to be fed our modified
 * data.
 */
void ChipStream::registerStream(ChipStream *stream) {
  m_Streams.push_back(stream);
}
  
/** 
 * @brief register our parent stream that will be passing
 * data to this object.
 * @param stream - who are we getting data from?
 */
void ChipStream::registerParent(ChipStream *stream) {
  m_ParentStream = stream;
}

/** 
 * @brief Get a reference to our parent stream.
 * @return Reference to parent stream.
 */
const ChipStream &ChipStream::getParent() {
  return *m_ParentStream;
}

/** 
 * @brief Method for being passed a new cel file worth of data.
 * @param iMart - repository of intensities for all CEL files in analysis
 */
void ChipStream::newDataSet(IntensityMart* iMart) {
  m_TransformedIMart = iMart;
  chipStreamPassNewChip(iMart);
}

void ChipStream::newChip(std::vector<float> &data) {};

void ChipStream::finishedChips() {};

/** 
 * @brief Signal that no more data is coming (i.e. newChip() will
 * not be called anymore.
 */
void ChipStream::endDataSet() {
  chipStreamEndChips();
}

/** 
 * @brief Set the minimum subset of a chip that this stream needs to
 * see for initialization.
 * @param subset - bitmask of probes.
 */
void ChipStream::setMinimumSubset(std::vector<bool> &subset) {
  unsigned int i = 0;
  for(i = 0; i < m_Streams.size(); i++) 
    m_Streams[i]->setMinimumSubset(subset);
}

/**
 * @brief Placeholder function for transform()
 *
 * @param probeIx - Probe index on chip.
 * @param chipIx - Set of chip indexes from same sample.
 * @param intensity - CEL intensity
 */
float ChipStream::transform(int probeIx, int chipIx, float intensity) {
  return intensity;
}

/**
 * @brief Placeholder function for transform()
 *
 * @param probeIx - Probe index on chip.
 * @param chipIx - Set of chip indexes from same sample.
 * @param intensity - CEL intensity
 * @param board - Blackboard with various state data
 */
float ChipStream::transform(int probeIx, int chipIx, float intensity, PsBoard &board) {
  return transform(probeIx, chipIx, intensity);
}

/**
 * @brief Retrieve intensity from resident IntensityMart
 *
 * @param probeIx - Probe index on chip.
 * @param chipIx - Set of chip indexes from same sample.
 * @param channelIx - CEL file channel to retrieve intensity from
 */
float ChipStream::getTransformedIntensity(probeidx_t probeIx, chipidx_t chipIx, unsigned int channelIx) {
  float intensity = 0.0;

  if (m_TransformedIMart != NULL) {
    intensity = m_TransformedIMart->getProbeIntensity(probeIx, chipIx, channelIx);
  }
  else {
    Err::errAbort("ChipStream::getTransformedIntensity -- associated IntensityMart is NULL");
  }
  return intensity;
}



/** 
 * @brief Set channelCount for this chipstream and all nodes downstream
 * 
 * @param channelCount - number of channels for this chipstream
 */
/*   virtual void setChannelCount(int channelCount) { */
/*     m_ChannelCount = channelCount; */
/*     for(int i = 0; i < m_Streams.size(); i++) { */
/*       m_Streams[i]->setChannelCount(channelCount); */
/*     } */
/*   } */

/** 
 * @brief Method for being passing a new cel file worth of data to children.
 * @param iMart - IntensityMart of intensities to be initialized
 */
void ChipStream::chipStreamPassNewChip(IntensityMart* iMart) {
  unsigned int i = 0;
  for(i = 0; i < m_Streams.size(); i++) {
    m_Streams[i]->newDataSet(iMart);
  }
}

/** 
 * @brief Signal children that no more data is coming (i.e. newChip() will
 * not be called anymore.
 */
void ChipStream::chipStreamEndChips() {
  if (m_Streams.size() > 0) {
    for(unsigned int i = 0; i < m_Streams.size(); i++) {
      m_Streams[i]->endDataSet();
    }
    m_TransformedIMart = NULL;
  }
}

