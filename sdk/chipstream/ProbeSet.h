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
 * @file   ProbeSet.h
 * @author Chuck Sugnet
 * @date   Mon May 16 14:29:19 2005
 *
 * @brief  Container class for probe sets.
 */

#ifndef _PROBESET_H_
#define _PROBESET_H_

//
#include "chipstream/AptTypes.h"
//
#include "util/Convert.h"
#include "util/Err.h"
#include "util/Util.h"
//
#include <cstring>
#include <iostream>
#include <string>
#include <vector>
//

// This was because of the '-1' being mapped to unsigned
// char. We should use signed chars.
#define NULLPROBEGC 254

/// these two macros should be used to make it clear what is happening
/// when probeids are mapped between the internal an
/// external worlds. should also note that "_id" is evaled
/// more than once.  (so dont use an expression.)

/// "0" is bad, map it to "-1", leave "<0" unchanged, subtract 1 from the rest.
#define TO_ZEROBASED(_id) ((0<=_id)?_id-1:_id)
/// Add 1 to everthing bigger than -1
#define FROM_ZEROBASED(_id) ((0<=_id)?_id+1:_id)

// Predeclarations as they are both pointing to each other
class ProbeSet;
class Atom;

/**
   @brief Represents an individual feature on a chip.

   <b>NOTE</b> The Probe, Atom and ProbeSet classes were determined to be
   wasting valuable RAM and are being transitioned to the ProbeList
   class. Please code accordingly.

   In principle a probe is a piece of dna of a particular sequence. In practice
   probes are treated differently if the are different x,y spots on the array
   even if they contain exactly the same sequence.

   For the CEL file SDK the order of the probes indexes are:
   \verbatim
   (0,0) (1,0) (2,0) (3,0) ... (numCols,0)
   (0,1) (1,1) (2,1) (3,1) ... (numCols,1)
   (0,2) (1,2) (2,2) (3,2) ... (numCols,2)
   ...
   (0,numRows)       ... (numCols,numRows)
   \endverbatim

   Where in cel file terminology cols are the x dimension and
   rows are the y dimension. So using the standard XYToIndex()
   function in the cel file API you can transform into a single
   vector index i using:
   - i = ((y*numColumns) + x)
   - x = i % numCols
   - y = i / numCols

   This index, i, is used as the probe id in the probe data structure as it is
   unique and provides an easy index into a vector that represent all of the
   data in the array. While it is 0 based there are files that are 1 based
   (notably the PGF file) and sometimes it is necessary to increment or
   decrement the id to accomedate the application or file format.

   So the index (and thus probe ids) numbering when num cols = 10
   \verbatim
   0  1  2  3  4  5  6  7  8  9
   10 11 12 13 14 15 16 17 18 19
   20 21 22 23 24 25 26 27 28 29
   \endverbatim

   Note that these indexes are only for the CEL files as the underlying DAT
   files can be rotated depending on the scanner and application.

 */
class Probe {

public:
  /// Type of probe.
  enum ProbeType {
    PMST      = 0,
    MMST      = 1,
    PMAT      = 2,
    MMAT      = 3,
    GENERICST = 4,
    GENERICAT = 5,
    JUMBOAT   = 6,
    JUMBOST   = 7,
    THERMOAT  = 8,
    THERMOST  = 9,
    TRIGRIDAT = 10,
    TRIGRIDST = 11,
    BLANK     = 12
  };

  /**
   * @brief Default constructor.
   */
  Probe();
  void clear();

  /**
   * @brief Construct probe from a list of strings.
   *
   * @param vec - Collection of strings.
   * @param offSet - Where in vector strings of intrest start.
   */
  Probe(const std::vector<std::string> &vec, unsigned int offSet=0);

  /**
   * @brief Construct probe from a list of char *
   *
   * @param vec - char * with data to make probe.
   * @param offSet - Offset of data in vector.
   */
  Probe(const std::vector<const char *> &vec, unsigned int offSet=0);

  /**
   * @brief String reprentation of probe type.
   * @param t Type of probe.
   *
   * @return Character description of probe type.
   */
  static const char *stringForType(Probe::ProbeType t);

  /**
   * @brief What type of probe does a string indicate?
   * @param cs - Description of probe type.
   *
   * @return Type of probe.
   */
  static ProbeType typeForString(const std::string& cs);

  static bool isPm(const Probe &p) {
    if(p.type == PMST || p.type == PMAT)
      return true;
    return false;
  }

  static bool isMm(const Probe &p) {
    if(p.type == MMST || p.type == MMAT)
      return true;
    return false;
  }

  /**
   * Allocate new memory for a probe and fill in values coppied from
   * supplied probe. Must delete this memory when finished.
   * @param p - probe to be cloned.
   * @return - new probe object allocated on heap, delete when done.
   */
  static Probe* cloneProbe(const Probe &p) {
    /// @todo this should be "probe->copy()"
    Probe *n = new Probe();
    n->id = p.id;
    n->type = p.type;
    n->gcCount = p.gcCount;
    n->m_apid=p.m_apid;

    //    n->parentAtom = NULL;
    return n;
  }

  /**
   * @brief Minimum number of fields that a probe requires.
   * @return Minimum number of fields
   */
  inline unsigned int numFields() {
    const static unsigned int numFields = 4;
    return numFields;
  }

  int getApid() const {
      return m_apid;
  }

  bool isNull();

  /* Data is public for simplicity and speed. */

  /// Id of probe, in practice usually realated to position on array.
  probeid_t id;
  /// Type of probe.
  unsigned char type;
  /// How many G/C nucleotides does the probe have? NULLPROBEGC is null value
  char gcCount;
  /// Atom that this probe is a part of. Memory not owned here
  // Atom* m_parentAtom;

  /// The AnalysisProbeId.
  /// This is set when the a ProbeSet extracted from a ProbeListFactory.
  /// -1 otherwise.
  int m_apid;
};

/**
 * @brief Represents a collection of probes that make up a unit for an
 * intensity measurement of a particular sequence. Typically a PM/MM
 * pair.
 *
 * <b>NOTE</b> The Probe, Atom and ProbeSet classes were determined to be
 * wasting valuable RAM and are being transitioned to the ProbeList
 *  class. Please code accordingly.
 */
class Atom {

public:
   /* Defines the type of probe replication within an atom   */
  enum ReplicationType {
    UnknownRepType,	// 0 Unspecified replication type
    DifferentRepType,	// 1 All probes have different sequences
    MixedRepType,	// 2 Some probes have identical sequences
    IdenticalRepType	// 3 All probes have identical sequences
  };

  /* Constructor. */
  Atom() : id(0), allele_code(0), context_code(0), channel_code(0), replication_type(UnknownRepType) {}

  /*Constructor. */
  Atom(atomid_t identity, int allele, int context, int channel, ReplicationType rep_type) {
    id = identity;
    allele_code = allele;
    context_code = context;
    channel_code = (char)channel;
    replication_type = rep_type;
  }

  /**
   * @brief Destructor.
   */
  ~Atom() {
    clear();
  }

  void clear() {
    for(size_t i = 0; i < probes.size(); i++) {
      delete probes[i];
      probes[i]=NULL;
    }
    probes.clear();
  }

  static Atom *deepCopy(Atom &a) {
    Atom *copy = new Atom();
    copy->id = a.id;
    copy->allele_code=a.allele_code;
    copy->context_code=a.context_code;
    copy->channel_code=a.channel_code;
    copy->replication_type=a.replication_type;

    //    copy->parentPSet = NULL;
    copy->probes.resize(a.probes.size(), NULL);
    std::vector<Probe *>::iterator pIx;
    std::vector<Probe *>::iterator cIx;
    for(pIx = a.probes.begin(), cIx = copy->probes.begin();
        pIx != a.probes.end() && cIx != copy->probes.end();
        ++pIx, ++cIx) {
      *cIx = Probe::cloneProbe(**pIx);
    }
    return copy;
  }

  /* Data is public for simplicity and speed. */
  /// Id of atom.
  atomid_t id;
  /// Code of allele (Ray generated?)
  int allele_code;
  /// Which SNP context is this atom in?  (in an arbratry order of contexts)
  int context_code;
  /// which channel is it on?
  char channel_code;
  /// what is the rep type (probes in atom are replicated, unique, or
  /// mixture of both?)
  ReplicationType replication_type;

  /// Vector of probes that make up this atom. Order is important
  std::vector<Probe *> probes;
  /// Probe set that this atom belongs to.
  //  ProbeSet* parentPSet;

  int getAlleleCode() {
    return allele_code;
  }
  int getContextCode() {
    return context_code;
  }
  int getChannelCode() {
    // %HACK: added this conditional to handle -1 channel code.  This
    // should go away when APT_PROBELISTPACKED_1 branch is checked in.
    // It's only here to help pass regression in the meantime.
    if (channel_code < 0) {
      channel_code = 0;
    }
    return channel_code;
  }
  int getReplicationType() {
    return replication_type;
  }
};

/**
 * @brief Collection of atoms that measure target that should be
 * absent or present at the same time.
 *
 * <b>NOTE</b> The Probe, Atom and ProbeSet classes were determined to be
 * wasting valuable RAM and are being transitioned to the ProbeList
 * class. Please code accordingly.
 */
class ProbeSet {

public:

  /** What sort of probeset is this? */
  /// @todo the code assumes that these enums match those in CDF file
  ///       we should use one or the other.
  enum ProbeSetType {
    Unknown,         // 0 Don't know this type of probeset.
    Expression,      // 1 mRNA expression measurement.
    GenoType,        // 2 Detecting SNPs, aka mapping arrays.
    Resequencing,    // 3 Determining entire sequences.
    Tag,             // 4 Never dealt with these.
    Copynumber,      // 5 How many chromosomes version
    GenoTypeCtrl,    // 6
    ExpressionCtrl,  // 7
    Marker,          // 8
    MultichannelMarker     // 9
  };

  /** Encoding of the pgf file probeset_type column, used by
      the qc report code to identify controls.
  */
  enum ControlType {
    NotControl,     // Not a control.
    NormGeneExon,   // Positive control.
    NormGeneIntron, // Negative control.
    ControlAffx,    // Generic control.
  };

  /**
      Snp chips with 4 blocks are defined to have 4 groups (blocks) of probes
      and the order is defined to be:
      group 1 - probes interrogating forward direction of A allele.
      group 2 - probes interrogating forward direction of B allele.
      group 3 - probes interrogating reverse direction of A allele.
      group 4 - probes interrogating reverse direction of B allele.
  */
  enum SnpOffsets {
    Aforward,
    Bforward,
    Areverse,
    Breverse,
  };

  /// Annotations a specific block.
  enum BlockAnnotation {
    NoAnnotation,
    A_Allele,
    B_Allele,
    A_AlleleForward,
    B_AlleleForward,
    A_AlleleReverse,
    B_AlleleReverse
  };

  /**
   * @brief Constructor.
   */
  ProbeSet() {
    name = NULL;
    psType = Unknown;
    numGroups = 0;
  }

  /**
   * @brief Destructor.
   */
  ~ProbeSet() {
    clear();
  }

  void clear() {
    delete [] name;
    name=NULL;
    for(size_t i = 0; i < atoms.size(); i++) {
      atoms[i]->clear();
      delete atoms[i];
      atoms[i]=NULL;
    }
    atoms.clear();
  };

  /**
   * @brief Generate probeset_type string from enum Type.
   */
  static std::string typeEnumToString(const ProbeSetType type) {
    std::string str;
    switch(type)
      {
      case Expression:
		  str = "Expression";
		  break;
      case GenoType:
		  str = "GenoType";
		  break;
      case Resequencing:
		  str = "Resequencing";
		  break;
      case Tag:
		  str = "Tag";
		  break;
      case Copynumber:
		  str = "Copynumber";
		  break;
      case GenoTypeCtrl:
                  str = "GenoTypeCtrl";
                  break;
      case ExpressionCtrl:
                  str = "ExpressionCtrl";
                  break;
      case Marker:
                  str = "Marker";
                  break;
      case MultichannelMarker:
                  str = "MultichannelMarker";
                  break;
      default:
		  str = "Unknown";
		  break;
    }
    return str;
  }

  /**
   * @brief Determine control type from the probeset_type column
   * of the .pgf file.
   */
  static enum ControlType controlStringToEnum(const std::string& probeset_type) {

    // Require "normgene" for positive and negative control types.
    size_t idx = probeset_type.find("normgene");
    if(idx != std::string::npos)
    {
      idx = probeset_type.find("exon");
      if(idx != std::string::npos)
        return NormGeneExon;

      idx = probeset_type.find("intron");
      if(idx != std::string::npos)
        return NormGeneIntron;
    }

    // Check for "control" and "affx".
    idx = probeset_type.find("control");
    if(idx != std::string::npos)
    {
      idx = probeset_type.find("affx");
      if(idx != std::string::npos)
        return ControlAffx;
    }

    return NotControl;
  }

  /**
   * @brief Generate probeset_type string from the ControlType.
   */
  static std::string controlEnumToString(const enum ControlType controlType) {
    std::string str;
    switch(controlType) {
      case NotControl:
	// Unknown - return null string.
	break;
      case NormGeneExon:
	// Example.
	str = "normgene->exon";
	break;
      case NormGeneIntron:
	str = "normgene->intron";
	break;
      case ControlAffx:
	str = "control->affx";
	break;
      default:
	break;
    }
    return str;
  }

  bool hasMM() const {
    for(uint32_t atomIx = 0; atomIx < atoms.size(); ++atomIx) {
      Atom *a = atoms[atomIx];
      for(uint32_t probeIx = 0; probeIx < a->probes.size(); ++probeIx) {
        if(Probe::isMm(*(a->probes[probeIx]))) {
          return true;
        }
      }
    }
    return false;
  }

  uint64_t getSizeOf() const {
    uint64_t sum = sizeof(ProbeSet);
    if(name != NULL) {
      sum += sizeof(char) * strlen(name);
      sum += sizeof(char *); // for allocator block
    }
    sum += sizeof(short) * atomsPerGroup.size();
    sum += sizeof(short *);
    sum += sizeof(Atom *) * ((atoms.size() * 2) + 1); // add one for new allocator block
    sum += sizeof(Atom) * atoms.size();
    std::vector<Atom *>::const_iterator aIx;
    for(aIx = atoms.begin(); aIx != atoms.end(); ++aIx) {
      const Atom *a = *aIx;
      sum += sizeof(Probe *) * ((a->probes.size() *2) + 1); // add one for new allocator block
      sum += sizeof(Probe) * a->probes.size();
    }
    return sum;
  }

  unsigned int countProbes() const {
    unsigned int count = 0;
    for(size_t a = 0; a < atoms.size(); a++) {
      count += (unsigned int)atoms[a]->probes.size();
    }
    return count;
  }

  unsigned int countPmProbes() const {
    unsigned int atomPmCount = 0;
    for(size_t atomIx = 0; atomIx < atoms.size(); atomIx++) {
      Atom &atom = *(atoms[atomIx]);
      for(size_t probeIx = 0; probeIx < atom.probes.size(); probeIx++) {
        Probe *p = atom.probes[probeIx];
        if(p->type == Probe::PMST || p->type == Probe::PMAT) {
          atomPmCount++;
        }
      }
    }
    return atomPmCount;
  }

  static ProbeSet *deepCopy(ProbeSet &ps) {
    ProbeSet *copy = new ProbeSet();
    if(ps.name != NULL)
      copy->name = Util::cloneString(ps.name);
    copy->psType = ps.psType;
    copy->numGroups = ps.numGroups;
    copy->atomsPerGroup.resize(ps.atomsPerGroup.size());
    for(uint32_t i = 0; i < ps.atomsPerGroup.size(); i++) {
      copy->atomsPerGroup[i] = ps.atomsPerGroup[i];
    }
    copy->atoms.resize(ps.atoms.size(), NULL);
    std::vector<Atom *>::iterator aIx;
    std::vector<Atom *>::iterator cIx;
    for(aIx = ps.atoms.begin(), cIx = copy->atoms.begin();
	aIx != ps.atoms.end() && cIx != copy->atoms.end();
	++aIx, ++cIx) {
      *cIx = Atom::deepCopy(**aIx);
      //      (*cIx)->parentPSet = copy;
    }
    return copy;
  }

  /// Name of probe set.
  char *name;
  /// Atoms that make up this probe set.
  std::vector<Atom *> atoms;
  /// What type of probeset is this?
  ProbeSetType psType;
  /// How many groups (or blocks) are there? Used extensively in genotyping chips.
  unsigned char numGroups;
  /// Number of atoms in each group, used to get the offsets into the atoms vector.
  std::vector<int> atomsPerGroup;
};
	
/**
 * @brief Group of probe sets that should be processed as a single
 * large probe set.
 */
class ProbeSetGroup {

public :

  /**
   * @brief Constructor for making singleton probe set group.
   * @param ps - Probe set to make singleton from.
   */
  ProbeSetGroup(ProbeSet *ps) {
    assert(ps && "ProbeSetGroup() - ps == NULL");
    ownMem = true;
	displayName = NULL;
    if(ps->name != NULL)
      name = Util::cloneString(ps->name);
    else
      name = NULL;
    probeSets.push_back(ps);
  }

  /**
   * @brief Constructor assures that name is NULL initially.
   */
  ProbeSetGroup() {
    name = NULL;
	displayName = NULL;
    ownMem = true;
  }

  /**
   * @brief Destructor.
   */
  ~ProbeSetGroup() {
    if(ownMem) {
      delete [] name;
      for(std::vector<char *>::size_type i = 0; i < probeSetNames.size(); i++)
        delete [] probeSetNames[i];
      for(std::vector<ProbeSet *>::size_type i = 0; i < probeSets.size(); i++)
        delete probeSets[i];
   }
  }

  inline bool pointersLoaded() const {
    return probeSets.size() > 0;
  }

  inline int probesetCount() const {
      int nnames = int(probeSetNames.size());
      return nnames;
  }

  /**
   * How many perfect match probes are there for a particular probeset
   * group?
   * @return - number of perfect match probes for this group and chip.
   */
  int countPmProbes() {
    unsigned int probeIx = 0, psIx, atomIx = 0, atomPmCount = 0;
    /* Count up the number of probes we need. */
    for(psIx = 0; psIx < probeSets.size(); psIx++) {
      const ProbeSet *ps = probeSets[psIx];
      if(ps == NULL) {
        Verbose::out(1,"Can't find probe set for id: " + ToStr(probeSets[psIx]));
        continue;
      }
      for(atomIx = 0; atomIx < ps->atoms.size(); atomIx++) {
        Atom &atom = *(ps->atoms[atomIx]);
        for(probeIx = 0; probeIx < atom.probes.size(); probeIx++) {
          Probe *p = atom.probes[probeIx];
          if(p->type == Probe::PMST || p->type == Probe::PMAT) {
            atomPmCount++;
          }
        }
      }
    }
    return atomPmCount;
  }

  /// Should we try to free memory in destructor?
  bool ownMem;
  /// Name of probe set group.
  char *name;
  const char *displayName;
  /// Vector of probe set pointers
  std::vector<ProbeSet *> probeSets;
  /// Vector of probe set names
  std::vector<char *> probeSetNames;
};

/// Output a probe.
inline std::ostream & operator<<(std::ostream &o, const Probe &p) {
  o << "\t\t" << p.id << '\t' << Probe::stringForType((enum Probe::ProbeType)p.type) << '\t' << p.gcCount << std::endl;
  return o;
}

/// Output an atom.
inline std::ostream & operator<<(std::ostream &o, const Atom &a) {
  unsigned int i = 0;
  o << '\t' << a.id << std::endl;
  for(i = 0; i < a.probes.size(); i++) {
    o << *(a.probes[i]) << '\t' << 0 << std::endl;
  }
  return o;
}
/// Output a probe set.
inline std::ostream & operator<<(std::ostream &o, const ProbeSet &p) {
  unsigned int i = 0;
//  o << p.name << '\t' << p.type << '\t' << std::endl;
  o << p.name << '\t' << std::endl;
  for(i = 0; i < p.atoms.size(); i++) {
    o << *(p.atoms[i]);
  }
  return o;
}

#endif /* _PROBESET_H_ */
