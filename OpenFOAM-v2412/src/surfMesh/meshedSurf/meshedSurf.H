/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::meshedSurf

Description
    Abstract definition of a meshed surface defined by faces and points.

\*---------------------------------------------------------------------------*/

#ifndef Foam_meshedSurf_H
#define Foam_meshedSurf_H

#include "pointField.H"
#include "faceList.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class meshedSurf Declaration
\*---------------------------------------------------------------------------*/

class meshedSurf
{
public:

    // Forward Declarations
    class emptySurface;


    // Constructors

        //- Default construct
        constexpr meshedSurf() noexcept = default;


    //- Destructor
    virtual ~meshedSurf() = default;


    // Member Functions

        //- The points used for the surface
        virtual const pointField& points() const = 0;

        //- The faces used for the surface
        virtual const faceList& faces() const = 0;

        //- Per-face zone/region information.
        //  Default is zero-sizes (ie, no zone information)
        virtual const labelList& zoneIds() const
        {
            return labelList::null();
        }

        //- Per-face identifier (eg, element Id)
        virtual const labelList& faceIds() const
        {
            return labelList::null();
        }
};


/*---------------------------------------------------------------------------*\
                  Class meshedSurf::emptySurface Declaration
\*---------------------------------------------------------------------------*/

//- A concrete meshedSurf class without faces, points, etc.
class meshedSurf::emptySurface
:
    public meshedSurf
{
public:

    // Constructors

        //- Default construct
        constexpr emptySurface() noexcept = default;


    //- Destructor
    virtual ~emptySurface() = default;


    // Member Functions

        //- Zero-sized point field
        virtual const pointField& points() const
        {
            return pointField::null();
        }

        //- Zero-sized face list
        virtual const faceList& faces() const
        {
            return faceList::null();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
