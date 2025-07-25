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
    Foam::cellBitSet

Description
    A special purpose topoSet with the cell labels stored as a bitSet.
    It does not correspond to a cellSet either (no associated IOobject).

    The cellBitSet::select method is embedded in function objects
    such as Foam::functionObjects::ensightWrite with actions like the
    following:

    \verbatim
    selection
    {
        box
        {
            action  use;
            source  box;
            box     (-0.1 -0.01 -0.1) (0.1 0.30 0.1);
        }
        dome
        {
            action  add;
            shape   sphere;
            origin  (-0.1 -0.01 -0.1);
            radius  0.25;
        }
        centre
        {
            action  subtract;
            source  sphere;
            origin  (-0.1 -0.01 -0.1);
            radius  0.1;
        }
        blob
        {
            action  add;
            source  surface;
            surface triSurfaceMesh;
            name    blob.stl;
        }
    }
    \endverbatim

SourceFiles
    cellBitSet.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_cellBitSet_H
#define Foam_cellBitSet_H

#include "topoBitSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class dictionary;

/*---------------------------------------------------------------------------*\
                           Class cellBitSet Declaration
\*---------------------------------------------------------------------------*/

class cellBitSet
:
    public topoBitSet
{
public:

    //- Runtime type information
    TypeNameNoDebug("cellBitSet");


    // Constructors

        //- Construct with nCells elements,
        //- all elements unset or initial value
        explicit cellBitSet(const polyMesh& mesh, const bool val = false);

        //- Copy construct from bitset, resizing to nCells elements as required
        cellBitSet(const polyMesh& mesh, const bitSet& bits);

        //- Move construct from bitset, resizing to nCells elements as required
        cellBitSet(const polyMesh& mesh, bitSet&& bits);


    //- Destructor
    virtual ~cellBitSet() = default;


    // Static Functions

        //- Return a cell selection according to the dictionary
        //- specification of actions.
        static bitSet select
        (
            const polyMesh& mesh,
            const dictionary& dict,
            const bool verbosity = false
        );


    // Member Functions

        //- Sync cellBitSet across coupled patches.
        virtual void sync(const polyMesh& mesh)
        {}

        //- Return max index+1.
        virtual label maxSize(const polyMesh& mesh) const;

        //- Update any stored data for new labels.
        virtual void updateMesh(const mapPolyMesh& morphMap)
        {}

        //- Update any stored data for mesh redistribution.
        virtual void distribute(const mapDistributePolyMesh& map)
        {}

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
