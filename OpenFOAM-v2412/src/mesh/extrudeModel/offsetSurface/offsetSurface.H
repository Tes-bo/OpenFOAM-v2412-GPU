/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::extrudeModels::offsetSurface

Description
    Extrudes by interpolating points from one surface to the other. Surfaces
    have to be topologically identical i.e. one has to be an offsetted version
    of the other.

    The offsetSurfaceCoeffs dictionary entries
    \table
        Property      | Description                         | Required | Default
        baseSurface   | Base surface file name              | yes |
        offsetSurface | Offset surface file name            | yes |
        project       | Project onto offsetted surface      | no  | false
    \endtable

\*---------------------------------------------------------------------------*/

#ifndef offsetSurface_H
#define offsetSurface_H

#include "point.H"
#include "extrudeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// Forward declarations
class triSurface;
class triSurfaceSearch;

namespace extrudeModels
{

/*---------------------------------------------------------------------------*\
                Class extrudeModels::offsetSurface Declaration
\*---------------------------------------------------------------------------*/

class offsetSurface
:
    public extrudeModel
{
    // Private data

        //- The base surface
        autoPtr<triSurface> baseSurfPtr_;

        //- The offset surface
        autoPtr<triSurface> offsetSurfPtr_;

        //- The search engine
        autoPtr<triSurfaceSearch> baseSearchPtr_;

        //- The search engine
        autoPtr<triSurfaceSearch> offsetSearchPtr_;

        //- Project onto offsetted surface?
        const bool project_;


public:

    //- Runtime type information
    TypeName("offsetSurface");

    // Constructors

        //- Construct from dictionary
        explicit offsetSurface(const dictionary& dict);


    //- Destructor
    virtual ~offsetSurface();


    // Member Operators

        //- Return point
        point operator()
        (
            const point& surfacePoint,
            const vector& surfaceNormal,
            const label layer
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace extrudeModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
