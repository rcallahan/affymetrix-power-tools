////////////////////////////////////////////////////////////////
//
// Copyright (C) 2009 Affymetrix, Inc.
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
#ifndef ARRAYELEMENT_H
#define ARRAYELEMENT_H

// ArrayLevel stores an array of ArrayElements.  Each ArrayElement
// corresponds to a summary of a block of residuals.  Specifically
// size -> the total number of probe cells in the block
// count -> how many probe cells in the block do not count as missing data
// mean -> the mean of the data in the block
// back -> an estimate of the background intensity

// This class is implemented as a struct.  The constructors are chosen
// on the basis of convenience of use in the ArrayLevel class


class ArrayElement {
public:
    ArrayElement() : mean(0.0), count(0), size(0), back(0.0) { }

    ArrayElement(int sz) : mean(0.0), count(0), size(sz), back(0.0) { }

    ArrayElement(double m, int n, int sz) :
        mean(m), count(n), size(sz), back(m) { }

    double mean;
    int count;
    int size;
    double back;
};

#endif
