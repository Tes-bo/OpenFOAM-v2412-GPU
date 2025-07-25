/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::patchInteractionData

Description
    Helper class for the LocalInteraction patch interaction model

\*---------------------------------------------------------------------------*/

#ifndef patchInteractionData_H
#define patchInteractionData_H

#include "wordRe.H"
#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declaration
class patchInteractionData;
Istream& operator>>(Istream& is, patchInteractionData& pid);

/*---------------------------------------------------------------------------*\
                    Class patchInteractionData Declaration
\*---------------------------------------------------------------------------*/

class patchInteractionData
{
    // Private Data

        //- Interaction type name
        word interactionTypeName_;

        //- Patch name(s)
        wordRe patchName_;

        //- Elasticity coefficient
        scalar e_;

        //- Restitution coefficient
        scalar mu_;


public:

    // Constructor

        //- Default construct
        patchInteractionData();


    // Member Functions

        //- The interaction type name
        const word& interactionTypeName() const noexcept
        {
            return interactionTypeName_;
        }

        //- The patch name(s)
        const wordRe& patchName() const noexcept
        {
            return patchName_;
        }

        //- The elasticity coefficient
        scalar e() const noexcept
        {
            return e_;
        }

        //- The restitution coefficient
        scalar mu() const noexcept
        {
            return mu_;
        }


    // IO Operators

        //- Istream operator
        friend Istream& operator>>(Istream& is, patchInteractionData& pid);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
