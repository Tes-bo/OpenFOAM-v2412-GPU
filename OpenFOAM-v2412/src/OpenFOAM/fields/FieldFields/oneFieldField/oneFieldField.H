/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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
    Foam::oneFieldField

Description
    A class representing the concept of a field of oneFields used to
    avoid unnecessary manipulations for objects which are known to be one at
    compile-time.

    Used for example as the density argument to a function written for
    compressible to be used for incompressible flow.

\*---------------------------------------------------------------------------*/

#ifndef Foam_oneFieldField_H
#define Foam_oneFieldField_H

#include "oneField.H"
#include "Pair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class oneFieldField Declaration
\*---------------------------------------------------------------------------*/

class oneFieldField
:
    public one
{
public:

    // Constructors

        //- Default construct
        oneFieldField() noexcept = default;


    // Member Operators

        //- Return the field
        oneField operator[](const label) const noexcept
        {
            return oneField{};
        }

        //- Return the field element
        one operator[](const labelPair&) const noexcept
        {
            return one{};
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
