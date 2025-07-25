/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013 OpenFOAM Foundation
    Copyright (C) 2020-2021 OpenCFD Ltd.
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
    Foam::coupleGroupIdentifier

Description
    Encapsulates using "patchGroups" to specify coupled patch

SourceFiles
    coupleGroupIdentifier.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_coupleGroupIdentifier_H
#define Foam_coupleGroupIdentifier_H

#include "word.H"
#include "label.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class dictionary;
class polyMesh;
class polyPatch;

/*---------------------------------------------------------------------------*\
                    Class coupleGroupIdentifier Declaration
\*---------------------------------------------------------------------------*/

class coupleGroupIdentifier
{
    // Private Data

        //- Name of patchGroup
        word name_;


    // Private Member Functions

        //- Find other patch in specified mesh.
        //  \return index of patch or -1
        label findOtherPatchID
        (
            const polyMesh& mesh,
            const polyPatch& thisPatch
        ) const;


public:

    // Generated Methods

        //- Default construct
        coupleGroupIdentifier() = default;


    // Constructors

        //- Construct from patchGroup name
        explicit coupleGroupIdentifier(const word& patchGroupName)
        :
            name_(patchGroupName)
        {}

        //- Construct from dictionary
        explicit coupleGroupIdentifier(const dictionary& dict);


    // Member Functions

        //- Name of patchGroup
        const word& name() const noexcept { return name_; }

        //- The patchGroup has a non-empty name
        bool good() const noexcept { return !name_.empty(); }

        //- The patchGroup has a non-empty name
        bool valid() const noexcept { return good(); }

        //- Find other patch in same region.
        //  \return index of patch or -1.
        label findOtherPatchID(const polyPatch& thisPatch) const;

        //- Find other patch and region.
        //  Returns index of patch and sets otherRegion to name of region.
        //  FatalError if patch not found
        label findOtherPatchID
        (
            const polyPatch& thisPatch,
            word& otherRegion
        ) const;

        //- Write the coupleGroup dictionary entry
        void write(Ostream& os) const;
};


// Global Operators

//- Write the coupleGroup dictionary entry
Ostream& operator<<(Ostream& os, const coupleGroupIdentifier& ident);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
