/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "axesRotation.H"
#include "dictionary.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace coordinateRotations
{

    defineTypeName(axes);

    // Standard short name
    addNamedToRunTimeSelectionTable
    (
        coordinateRotation,
        axes,
        dictionary,
        axes
    );

    // Long name - Compat 1806
    addAliasToRunTimeSelectionTable
    (
        coordinateRotation,
        axes,
        dictionary,
        axes,
        axesRotation,
        1806
    );

} // End namespace coordinateRotations
} // End namespace Foam


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::tensor Foam::coordinateRotations::axes::rotation
(
    const vector& axis1,
    const vector& axis2,
    axisOrder order
)
{
    const scalar magAxis1(mag(axis1));
    scalar magAxis2(mag(axis2));

    if (magAxis1 < ROOTVSMALL)
    {
        FatalErrorInFunction
            << "Dominant coordinate axis cannot have zero length"
            << nl << endl
            << abort(FatalError);
    }

    const vector ax1(axis1 / magAxis1);  // normalise
    vector ax2(axis2);

    if (magAxis2 < ROOTVSMALL)
    {
        // axis2 == Zero : Use best-guess for the second axis.
        ax2 = findOrthogonal(axis1);
    }

    ax2.removeCollinear(ax1);

    magAxis2 = mag(ax2);

    if (magAxis2 < SMALL)
    {
        WarningInFunction
            << "axis1, axis2 appear to be co-linear: "
            << axis1 << ", " << axis2 << "  Revert to guessing axis2"
            << nl << endl;

        ax2 = findOrthogonal(axis1);

        ax2.removeCollinear(ax1);

        magAxis2 = mag(ax2);

        if (magAxis2 < SMALL)
        {
            FatalErrorInFunction
                << "Could not find an appropriate second axis"
                << nl << endl
                << abort(FatalError);
        }
    }

    ax2 /= magAxis2;  // normalise


    // The local axes are columns of the rotation matrix

    tensor rotTensor;

    switch (order)
    {
        case E1_E2:
        {
            rotTensor.col<0>(ax1);
            rotTensor.col<1>(ax2);
            rotTensor.col<2>(ax1^ax2);
            break;
        }
        case E2_E3:
        {
            rotTensor.col<0>(ax1^ax2);
            rotTensor.col<1>(ax1);
            rotTensor.col<2>(ax2);
            break;
        }
        case E3_E1:
        case E3_E1_COMPAT:
        {
            rotTensor.col<0>(ax2);
            rotTensor.col<1>(ax1^ax2);
            rotTensor.col<2>(ax1);
            break;
        }
    }

    return rotTensor;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::coordinateRotations::axes::read(const dictionary& dict)
{
    if
    (
        dict.readIfPresent("e1", axis1_)
     && dict.readIfPresent("e2", axis2_)
    )
    {
        order_ = E1_E2;
    }
    else if
    (
        dict.readIfPresent("e2", axis1_)
     && dict.readIfPresent("e3", axis2_)
    )
    {
        order_ = E2_E3;
    }
    else if
    (
        dict.readIfPresent("e3", axis1_)
     && dict.readIfPresent("e1", axis2_)
    )
    {
        order_ = E3_E1;
    }
    else if
    (
        dict.readIfPresent("axis", axis1_)
     && dict.readIfPresent("direction", axis2_)
    )
    {
        order_ = E3_E1_COMPAT;
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "No entries of the type (e1, e2) or (e2, e3) or (e3, e1) found"
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coordinateRotations::axes::axes()
:
    coordinateRotation(),
    axis1_(0,0,1),  // e3 = global Z
    axis2_(1,0,0),  // e1 = global X
    order_(E3_E1)
{}


Foam::coordinateRotations::axes::axes(const axes& crot)
:
    coordinateRotation(),
    axis1_(crot.axis1_),
    axis2_(crot.axis2_),
    order_(crot.order_)
{}


Foam::coordinateRotations::axes::axes
(
    const vector& axis1,
    const vector& axis2,
    axisOrder order
)
:
    coordinateRotation(),
    axis1_(axis1),
    axis2_(axis2),
    order_(order)
{}


Foam::coordinateRotations::axes::axes(const vector& axis)
:
    coordinateRotations::axes(axis, Zero, E3_E1_COMPAT)  // Guess second axis
{}


Foam::coordinateRotations::axes::axes(const dictionary& dict)
:
    coordinateRotations::axes()
{
    read(dict);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::coordinateRotations::axes::clear()
{
    axis1_ = vector(0,0,1);  // primary axis (e3, global Z)
    axis2_ = vector(1,0,0);  // secondary axis (e1, global X)
    order_ = E3_E1;
}


Foam::tensor Foam::coordinateRotations::axes::R() const
{
    return axes::rotation(axis1_, axis2_, order_);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::coordinateRotations::axes::write(Ostream& os) const
{
    switch (order_)
    {
        case E1_E2:
            os << "e1: " << axis1_ << " e2: " << axis2_;
            break;
        case E2_E3:
            os << "e2: " << axis1_ << " e3: " << axis2_;
            break;
        case E3_E1:
            os << "e1: " << axis2_ << " e3: " << axis1_;
            break;
        case E3_E1_COMPAT:
            os << "axis: " << axis1_ << " direction: " << axis2_;
            break;
    }
}


void Foam::coordinateRotations::axes::writeEntry
(
    const word& keyword,
    Ostream& os
) const
{
    // We permit direct embedding of the axes specification without
    // requiring a sub-dictionary.

    const bool subDict = !keyword.empty();

    if (subDict)
    {
        os.beginBlock(keyword);
        os.writeEntry("type", type());
    }

    switch (order_)
    {
        case E1_E2:
        {
            os.writeEntry("e1", axis1_);
            os.writeEntry("e2", axis2_);
            break;
        }
        case E2_E3:
        {
            os.writeEntry("e2", axis1_);
            os.writeEntry("e3", axis2_);
            break;
        }
        case E3_E1:
        {
            os.writeEntry("e1", axis2_);
            os.writeEntry("e3", axis1_);
            break;
        }
        case E3_E1_COMPAT:
        {
            os.writeEntry("axis", axis1_);
            os.writeEntry("direction", axis2_);
            break;
        }
    }

    if (subDict)
    {
        os.endBlock();
    }
}


// ************************************************************************* //
