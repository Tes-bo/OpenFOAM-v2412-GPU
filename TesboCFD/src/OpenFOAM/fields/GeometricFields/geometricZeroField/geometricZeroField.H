/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::geometricZeroField

Description
    A class representing the concept of a GeometricField of 0 used to avoid
    unnecessary manipulations for objects which are known to be zero at
    compile-time.

    Used for example as the density argument to a function written for
    compressible to be used for incompressible flow.

\*---------------------------------------------------------------------------*/

#ifndef Foam_geometricZeroField_H
#define Foam_geometricZeroField_H

#include "zeroFieldField.H"
#include "dimensionSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class geometricZeroField Declaration
\*---------------------------------------------------------------------------*/

class geometricZeroField
:
    public zero
{
public:

    // Public Typedefs

        typedef zeroField Internal;
        typedef zeroField Patch;
        typedef zeroFieldField Boundary;
        typedef zero cmptType;


    // Constructors

        //- Default construct
        geometricZeroField() noexcept = default;


    // Member Functions

        const dimensionSet& dimensions() const noexcept
        {
            return dimless;
        }

        zeroField field() const noexcept
        {
            return zeroField{};
        }

        zeroField oldTime() const noexcept
        {
            return zeroField{};
        }

        Internal internalField() const noexcept
        {
            return Internal{};
        }

        zeroField primitiveField() const noexcept
        {
            return zeroField{};
        }

        Boundary boundaryField() const noexcept
        {
            return Boundary{};
        }

        // Same as internalField()
        Internal v() const noexcept
        {
            return Internal{};
        }


    // Member Operators

        scalar operator[](const label) const noexcept
        {
            return 0;
        }

        Internal operator()() const noexcept
        {
            return Internal{};
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Global Operators

#include "geometricZeroFieldI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
