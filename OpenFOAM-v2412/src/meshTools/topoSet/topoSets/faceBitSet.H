/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2024 OpenCFD Ltd.
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
    Foam::faceBitSet

Description
    A special purpose topoSet with the face labels stored as a bitSet.
    It does not correspond to a faceSet either (no associated IOobject).

SourceFiles
    faceBitSet.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_faceBitSet_H
#define Foam_faceBitSet_H

#include "topoBitSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class faceBitSet Declaration
\*---------------------------------------------------------------------------*/

class faceBitSet
:
    public topoBitSet
{
public:

    //- Runtime type information
    TypeNameNoDebug("faceBitSet");


    // Constructors

        //- Construct with nFaces elements,
        //- all elements unset or initial value
        explicit faceBitSet(const polyMesh& mesh, const bool val = false);

        //- Copy construct from bitset, resizing to nFaces elements as required
        faceBitSet(const polyMesh& mesh, const bitSet& bits);

        //- Move construct from bitset, resizing to nFaces elements as required
        faceBitSet(const polyMesh& mesh, bitSet&& bits);


    //- Destructor
    virtual ~faceBitSet() = default;


    // Member Functions

        //- Sync faceBitSet across coupled patches.
        virtual void sync(const polyMesh& mesh);

        //- Return max index+1.
        virtual label maxSize(const polyMesh& mesh) const;

        //- Update any stored data for new labels.
        virtual void updateMesh(const mapPolyMesh& morphMap);

        //- Update any stored data for mesh redistribution.
        virtual void distribute(const mapDistributePolyMesh& map);

        //- Write maxLen items with label and coordinates.
        virtual void writeDebug
        (
            Ostream& os,
            const primitiveMesh& mesh,
            const label maxLen
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
