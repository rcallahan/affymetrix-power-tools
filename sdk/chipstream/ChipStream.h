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
 * @file   ChipStream.h
 * @author Chuck Sugnet
 * @date   Fri Oct 21 18:04:46 2005
 * 
 * @brief Base class for objects that can initialize themselves one at
 * a time and transform intensity data.
 */
#ifndef _CHIPSTREAM_H_
#define _CHIPSTREAM_H_

//
#include "chipstream/AptTypes.h"
#include "chipstream/IntensityMart.h"
#include "chipstream/IntensityTransformer.h"
#include "chipstream/SelfCreate.h"
#include "chipstream/SelfDoc.h"
#include "chipstream/PsBoard.h"
//
#include <cstdio>
#include <vector>
//

/**
 *  ChipStream - base class for objects that can initialize
 * themselves one at a time and transform intensity data.
 */
class ChipStream : public IntensityTransformer, public SelfDoc, public SelfCreate {

public:

  // constructor
  ChipStream();

  /** 
   * @brief Virtual destructor for a virtual class.
   */
  virtual ~ChipStream();

  /** 
   * Fill in the information for Self documentation.
   * @param doc - Self documenter to be filled in.
   */
  static void setupSelfDoc(SelfDoc &doc);

  /** 
   * @brief Supply a little how/what/why about the algorithms this
   * class performs and what parameters it takes.
   * @return SelfDoc
   */
  static SelfDoc explainSelf();

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
  static SelfCreate *newObject(std::map<std::string,std::string> &param);

  /** 
   * @Brief Add a stream to the list of those that will receive
   * downstream data.
   * @param stream - ChipStream that wants to be fed our modified
   * data.
   */
  virtual void registerStream(ChipStream *stream);
  
  /** 
   * @brief register our parent stream that will be passing
   * data to this object.
   * @param stream - who are we getting data from?
   */
  virtual void registerParent(ChipStream *stream);

  /** 
   * @brief Get a reference to our parent stream.
   * @return Reference to parent stream.
   */
  virtual const ChipStream &getParent();

  /** 
   * @brief Method for being passed a new cel file worth of data.
   * @param iMart - repository of intensities for all CEL files in analysis
   */
  virtual void newDataSet(IntensityMart* iMart);

  virtual void newChip(std::vector<float> &data);

  virtual void finishedChips();
  

  /** 
   * @brief Signal that no more data is coming (i.e. newChip() will
   * not be called anymore.
   */
  virtual void endDataSet();

  /** 
   * @brief Set the minimum subset of a chip that this stream needs to
   * see for initialization.
   * @param subset - bitmask of probes.
   */
  virtual void setMinimumSubset(std::vector<bool> &subset);

  /**
   * @brief Placeholder function for transform()
   *
   * @param probeIx - Probe index on chip.
   * @param chipIx - Set of chip indexes from same sample.
   * @param intensity - CEL intensity
   */
  virtual float transform(int probeIx, int chipIx, float intensity);

  /**
   * @brief Placeholder function for transform()
   *
   * @param probeIx - Probe index on chip.
   * @param chipIx - Set of chip indexes from same sample.
   * @param intensity - CEL intensity
   * @param board - Blackboard with various state data
   */
  virtual float transform(int probeIx, int chipIx, float intensity, PsBoard &board);

  /**
   * @brief Retrieve intensity from resident DiskIntensityMart
   *
   * @param probeIx - Probe index on chip.
   * @param chipIx - Set of chip indexes from same sample.
   * @param channelIx - CEL file channel to retrieve intensity from
   */
  float getTransformedIntensity(probeidx_t probeIx, chipidx_t chipIx, unsigned int channelIx = 0);

  /** 
   * Do any setup specific to this ChipStream object using the data and state specified
   * in the blackboar.
   * 
   * @param board Blackboard with data and state.
   */
  virtual void setParameters(PsBoard &board) {
    // no-op for now, some ChipStream objects don't require so don't make pure virtual
  }
  
protected:

  /** 
   * @brief Method for being passing a new cel file worth of data to children.
   * @param iMart - IntensityMart of intensities to be initialized
   */
  void chipStreamPassNewChip(IntensityMart* iMart);


  /** 
   * @brief Signal children that no more data is coming (i.e. newChip() will
   * not be called anymore.
   */
  void chipStreamEndChips();

  /// All of the registered chip streams.
  std::vector<ChipStream *> m_Streams;
  /// Parent of this stream.
  ChipStream *m_ParentStream;
  /// if newChip() is called with a IntensityMart, then
  /// transformed data is also stored in a IntensityMart called
  /// with one
  IntensityMart* m_TransformedIMart;
}; 

#endif /* CHIPSTREAM_H */
