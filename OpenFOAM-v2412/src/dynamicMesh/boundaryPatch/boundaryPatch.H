/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020-2023 OpenCFD Ltd.
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
    Foam::boundaryPatch

Description
    Like polyPatch but without reference to mesh.
    Used in boundaryMesh to hold data on patches.
    The patchIdentifier::index is set, but not used.

SourceFiles
    boundaryPatch.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_boundaryPatch_H
#define Foam_boundaryPatch_H

#include "patchIdentifier.H"
#include "labelRange.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class boundaryPatch Declaration
\*---------------------------------------------------------------------------*/

class boundaryPatch
:
    public patchIdentifier
{
    // Private Data

        //- The size of the patch
        label size_;

        //- Start label of the patch
        label start_;

public:

    // Generated Methods: copy construct, copy assignment


    // Constructors

        //- Construct from components
        boundaryPatch
        (
            const word& name,
            const label index,
            const label size,
            const label start,
            const word& physicalType = word::null
        );

        //- Construct from dictionary
        boundaryPatch
        (
            const word& name,
            const dictionary& dict,
            const label index
        );

        //- Copy construct, resetting the index
        boundaryPatch(const boundaryPatch& p, const label index);


        //- Clone
        autoPtr<boundaryPatch> clone() const
        {
            return autoPtr<boundaryPatch>::New(*this);
        }


    // Member Functions

        //- The start of the patch
        label start() const noexcept { return start_; }

        //- The start of the patch (modifiable)
        label& start() noexcept { return start_; }

        //- The size of the patch
        label size() const noexcept { return size_; }

        //- The size of the patch (modifiable)
        label& size() noexcept { return size_; }

        //- Return start/size range of this patch
        labelRange range() const noexcept
        {
            return labelRange(start_, size_);
        }


        //- Write dictionary entries (without surrounding braces)
        virtual void write(Ostream& os) const;
};


// Global Operators

//- Write boundaryPatch as dictionary entries (without surrounding braces)
Ostream& operator<<(Ostream&, const boundaryPatch& p);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
