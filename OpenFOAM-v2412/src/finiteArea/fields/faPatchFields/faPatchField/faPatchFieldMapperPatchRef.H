/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
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
    Foam::faPatchFieldMapperPatchRef

Description

\*---------------------------------------------------------------------------*/

#ifndef faPatchFieldMapperPatchRef_H
#define faPatchFieldMapperPatchRef_H

#include "faPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class faPatchFieldMapperPatchRef Declaration
\*---------------------------------------------------------------------------*/

class faPatchFieldMapperPatchRef
:
    public faPatchFieldMapper
{
    // Private Data

        const faPatch& sourcePatch_;
        const faPatch& targetPatch_;


    // Private Member Functions

        //- No copy construct
        faPatchFieldMapperPatchRef(const faPatchFieldMapperPatchRef&) = delete;

        //- No copy assignment
        void operator=(const faPatchFieldMapperPatchRef&) = delete;


public:

    // Constructors

        //- Construct from components
        faPatchFieldMapperPatchRef(const faPatch& source, const faPatch& target)
        :
            sourcePatch_(source),
            targetPatch_(target)
        {}


    // Member functions

        const faPatch& sourcePatch() const
        {
            return sourcePatch_;
        }

        const faPatch& targetPatch() const
        {
            return targetPatch_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
