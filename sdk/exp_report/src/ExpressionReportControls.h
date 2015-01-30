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

#ifndef _ExpressionReportControls_HEADER_
#define _ExpressionReportControls_HEADER_

/*! \file ExpressionReportControls.h Defines classes to store information about 3/m/5 prime controls. */

#include <cstring>
#include <list>
#include <string>
//

namespace ExpressionReport
{

/*! This class is used to store the index to the 3', M and 5' control probe sets. */
class ExpressionControl
{
    /*! The number of sets per control (3', 5', middle) */
	#define NUM_SETS_PER_CONTROL 3

	/*! Enumeration to describe the 3', M and 5' probe sets */
	public: enum ControlValueType { THREE_PRIME_PROBE_SET, MIDDLE_PROBE_SET, FIVE_PRIME_PROBE_SET };

	/*! The name of the control */
	private: std::string name;

	/*! Gets/Sets the name of the control. */
	public: std::string &Name() { return name; }

	/*! The indicies to the CHP file for the control probe sets. */
    private: int psIndex[NUM_SETS_PER_CONTROL];

    /*! Returns the index to the CHP file for the control probe set.
	 * @param valueType The control type
	 * @return The index to the control probe set
	 */
    public: int GetProbeSetIndex(ControlValueType valueType)
    {
        return psIndex[(int)valueType];
    }

    /*! Sets the index to the CHP file for the control probe set.
	 * @param valueType The control type
	 * @param index The index to the control probe set.
	 */
    public: void SetProbeSetIndex(ControlValueType valueType, int index)
    {
        psIndex[(int)valueType] = index;
    }

	/*! Comparison operator
	 * @param n Control name to compare
	 * @return Equality flag
	 */
	public: bool operator == (const std::string &n )
	{
		return (name == n);
	}

	/*! Inequality operator
	 * @param n Control name to compare
	 * @return Inequality flag
	 */
	public: bool operator != (const std::string &n )
	{
		return (name != n);
	}

	/*! Comparison operator
	 * @param c Control to compare
	 * @return Equality flag
	 */
	public: bool operator == (ExpressionControl &c )
	{
		return (name == c.Name());
	}

	/*! Inequality operator
	 * @param c Control to compare
	 * @return Inequality flag
	 */
	public: bool operator != ( ExpressionControl &c )
	{
		return (name != c.Name());
	}

	/*! Returns a flag indicating that the probe set index is valid.
	 * @param valueType The control type
	 * @return Flag indicating valid index
	 */
	public: bool HasValue(ControlValueType valueType)
	{
        return (psIndex[(int)valueType] != -1);
	}

	/*! Initializes the class members. */
	public: ExpressionControl()
	{
		Clear();
	}

	/*! Initializes the class members. */
	public: void Clear()
	{
		name = "";
		for (int i=0; i<NUM_SETS_PER_CONTROL; i++)
		{
            psIndex[i] = -1;
		}
	}
};

/*! An array of controls. */
typedef std::list<ExpressionControl> ExpressionControlList;

/*! This class contains a list of control objects. */
class ExpressionControls
{
	/*! 
	The name of the probe array type associated with the list of controls.
	*/
	private: std::string probeArrayType;

	/*! 
	The name of the probe array type associated with the list of controls.
	*/
	public: std::string &ProbeArrayType() { return probeArrayType; }

	/*! Equality operator
	 * @param controls Control to compare
	 * @param arrayType Name to compare
	 * @return Equality flag
	 */
	public: bool operator == ( const std::string &arrayType )
	{
		return (probeArrayType == arrayType);
	}

	/*! Equality operator
	 * @param c Control to compare
	 * @return Equality flag
	 */
	public: bool operator == ( ExpressionControls &c )
	{
		return (probeArrayType == c.ProbeArrayType());
	}

	/*! Inequality operator
	 * @param c Control to compare
	 * @return Inequality flag
	 */
	public: bool operator != ( ExpressionControls &c )
	{
		return (probeArrayType != c.ProbeArrayType());
	}

    /*! The array of spike controls. */
    private: ExpressionControlList spikeControls;

    /*! The array of spike controls. */
	public: ExpressionControlList &SpikeControls() { return spikeControls; }
    
    /*! The array of housekeeping controls. */
    private: ExpressionControlList housekeepingControls;

    /*! The array of housekeeping controls. */
	public: ExpressionControlList &HousekeepingControls() { return housekeepingControls; }

	/*! Initializes the members.
	 * @param arrayType The name of the probe array type
	 */
	public: ExpressionControls(const std::string &arrayType)
	{
		probeArrayType = arrayType;
    }

	/*! Initializes the members. */
	public: ExpressionControls()
	{
		Clear();
	}

	/*! Initializes the members. */
	public: void Clear()
	{
		probeArrayType = "";
		housekeepingControls.clear();
		spikeControls.clear();
    }
};

/*! The list of probe array type controls. */
typedef std::list<ExpressionControls> ExpressionControlsList;

}

#endif

