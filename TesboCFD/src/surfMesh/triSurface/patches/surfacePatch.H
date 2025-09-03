/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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
    Foam::surfacePatch

Description
    A 'patch' on surface as subset of triSurface.

SourceFiles
    surfacePatch.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_surfacePatch_H
#define Foam_surfacePatch_H

#include "geometricSurfacePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class surfZone;

/*---------------------------------------------------------------------------*\
                        Class surfacePatch Declaration
\*---------------------------------------------------------------------------*/

class surfacePatch
:
    public geometricSurfacePatch
{
    // Private Data

        //- Size of this patch in the face list
        label size_;

        //- Start label of this patch in the face list
        label start_;


public:

    // Generated Methods

        //- Copy construct
        surfacePatch(const surfacePatch&) = default;

        //- Copy assignment
        surfacePatch& operator=(const surfacePatch&) = default;


    // Constructors

        //- Default construct, with zero start, size, index=-1
        surfacePatch();

        //- Default construct, with zero start, size, and specified index
        explicit surfacePatch(const label index);

        //- Construct from components
        surfacePatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const word& geometricType = word::null
        );

        //- Construct from dictionary
        surfacePatch
        (
            const word& name,
            const dictionary& dict,
            const label index
        );


    // Member Functions

        //- Return start label of this patch in the face list
        label start() const noexcept
        {
            return start_;
        }

        //- Return start label of this patch in the face list
        label& start() noexcept
        {
            return start_;
        }

        //- Return size of this patch in the face list
        label size() const noexcept
        {
            return size_;
        }

        //- Return size of this patch in the face list
        label& size() noexcept
        {
            return size_;
        }

        //- Write dictionary, includes surrounding braces
        //  \warning Prior to 2020-01 was identical to operator<< output
        void write(Ostream& os) const;


    // Member Operators

        //- Conversion to surfZone representation
        explicit operator surfZone() const;


    // Housekeeping

        //- Removed(2020-01) Construct from Istream
        //  \deprecated(2020-01) - unused, inconsistent
        surfacePatch(Istream& is, const label index) = delete;

        //- Deprecated(2020-01) Construct from components
        //  \deprecated(2020-01) - order inconsistent with other identifiers
        FOAM_DEPRECATED_FOR(2020-01, "construct name/size/start/index/type")
        surfacePatch
        (
            const word& geometricType,
            const word& name,
            const label size,
            const label start,
            const label index
        )
        :
            surfacePatch(name, size, start, index, geometricType)
        {}

        //- Deprecated(2020-01) Ostream output
        //  \deprecated(2020-01) - Ostream output
        FOAM_DEPRECATED_FOR(2020-01, "write() or operator<<")
        void writeDict(Ostream& os) const
        {
            write(os);
        }
};


// Global Operators

//- Compare patches for equality
bool operator==(const surfacePatch& a, const surfacePatch& b);

//- Compare patches for inequality
bool operator!=(const surfacePatch& a, const surfacePatch& b);

//- Write name, geometricType, size, start
Ostream& operator<<(Ostream& os, const surfacePatch& obj);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
