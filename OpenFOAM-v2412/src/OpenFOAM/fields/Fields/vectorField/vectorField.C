/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "vectorField.H"

// * * * * * * * * * * * * * * * Specializations * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Note: the mag of vector and division is written out to avoid any casting
// between float and double.
//
// This enables specialization for floatVector and doubleVector independent
// of the definition of 'scalar' or 'vector' - useful for mixed-precision
// operation.

template<>
void Field<Vector<float>>::normalise()
{
    typedef float cmptType;

    constexpr float tol = floatScalarROOTVSMALL;

    for (Vector<cmptType>& v : *this)
    {
        // Foam::mag but using cmptType instead of scalar
        cmptType s
        (
            ::sqrt(magSqr(v.x()) + magSqr(v.y()) + magSqr(v.z()))
        );

        if (s < tol)
        {
            v.x() = 0;
            v.y() = 0;
            v.z() = 0;
        }
        else
        {
            v.x() /= s;
            v.y() /= s;
            v.z() /= s;
        }
    }
}


template<>
void Field<Vector<double>>::normalise()
{
    typedef double cmptType;

    constexpr double tol = doubleScalarROOTVSMALL;

    for (Vector<cmptType>& v : *this)
    {
        // Foam::mag but using cmptType instead of scalar
        cmptType s
        (
            ::sqrt(magSqr(v.x()) + magSqr(v.y()) + magSqr(v.z()))
        );

        if (s < tol)
        {
            v.x() = 0;
            v.y() = 0;
            v.z() = 0;
        }
        else
        {
            v.x() /= s;
            v.y() /= s;
            v.z() /= s;
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
