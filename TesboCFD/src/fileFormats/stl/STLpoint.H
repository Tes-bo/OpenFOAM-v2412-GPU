/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

Class
    Foam::STLpoint

Description
    A vertex point or facet normal representation for STL files.

\*---------------------------------------------------------------------------*/

#ifndef Foam_STLpoint_H
#define Foam_STLpoint_H

#include "point.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class STLpoint Declaration
\*---------------------------------------------------------------------------*/

class STLpoint
:
    public Vector<float>
{
public:

    // Constructors

        //- Default construct
        STLpoint() = default;

        //- Construct from single-precision point
        STLpoint(const Vector<float>& p)
        :
            Vector<float>(p)
        {}

        //- Construct from double-precision point
        STLpoint(const Vector<double>& p)
        :
            Vector<float>(float(p.x()), float(p.y()), float(p.z()))
        {}

        //- Construct from single-precision components
        STLpoint(float x, float y, float z)
        :
            Vector<float>(x, y, z)
        {}

        //- Construct from double-precision components
        STLpoint(double x, double y, double z)
        :
            Vector<float>(float(x), float(y), float(z))
        {}


    // Member Operators

        #ifdef WM_DP
        //- Conversion to double-precision point
        operator point() const
        {
            return point(x(), y(), z());
        }
        #endif
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
