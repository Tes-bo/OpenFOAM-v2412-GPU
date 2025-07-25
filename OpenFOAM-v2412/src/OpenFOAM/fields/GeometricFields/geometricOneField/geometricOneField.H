/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020-2022 OpenCFD Ltd.
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
    Foam::geometricOneField

Description
    A class representing the concept of a GeometricField of 1 used to avoid
    unnecessary manipulations for objects which are known to be one at
    compile-time.

    Used for example as the density argument to a function written for
    compressible to be used for incompressible flow.

\*---------------------------------------------------------------------------*/

#ifndef Foam_geometricOneField_H
#define Foam_geometricOneField_H

#include "oneFieldField.H"
#include "dimensionSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class geometricOneField Declaration
\*---------------------------------------------------------------------------*/

class geometricOneField
:
    public one
{
public:

    // Public Typedefs

        typedef oneField Internal;
        typedef oneField Patch;
        typedef oneFieldField Boundary;
        typedef one cmptType;


    // Constructors

        //- Default construct
        geometricOneField() noexcept = default;


    // Member Functions

        const dimensionSet& dimensions() const noexcept
        {
            return dimless;
        }

        oneField field() const noexcept
        {
            return oneField{};
        }

        oneField oldTime() const noexcept
        {
            return oneField{};
        }

        Internal internalField() const noexcept
        {
            return Internal{};
        }

        oneField primitiveField() const noexcept
        {
            return oneField{};
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

        one operator[](const label) const noexcept
        {
            return one{};
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

#include "geometricOneFieldI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
