/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2016-2024 OpenCFD Ltd.
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
    Foam::cellSet

Description
    A collection of cell labels.

SourceFiles
    cellSet.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_cellSet_H
#define Foam_cellSet_H

#include "topoSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class cellSet Declaration
\*---------------------------------------------------------------------------*/

class cellSet
:
    public topoSet
{
    // Private Member Functions

        //- No copy construct
        cellSet(const cellSet&) = delete;


public:

    //- Runtime type information
    TypeNameNoDebug("cellSet");


    // Constructors

        //- Construct from IOobject. No checking.
        explicit cellSet(const IOobject& io);

        //- Construct empty (no-read) with IOobject information
        cellSet(const IOobject& io, const Foam::zero);

        //- Construct from polyMesh (registry) and name.
        //- Checks for valid cell ids.
        cellSet
        (
            const polyMesh& mesh,
            const word& name,
            IOobjectOption::readOption rOpt = IOobjectOption::MUST_READ,
            IOobjectOption::writeOption wOpt = IOobjectOption::NO_WRITE,
            IOobjectOption::registerOption reg = IOobjectOption::LEGACY_REGISTER
        );

        //- Construct from polyMesh (registry), name and registration option
        cellSet
        (
            const polyMesh& mesh,
            const word& name,
            IOobjectOption::registerOption reg,
            IOobjectOption::readOption rOpt = IOobjectOption::MUST_READ,
            IOobjectOption::writeOption wOpt = IOobjectOption::NO_WRITE
        )
        :
            cellSet(mesh, name, rOpt, wOpt, reg)
        {}

        //- Construct empty (no-read) with initial labelHashSet capacity
        cellSet
        (
            const polyMesh& mesh,
            const word& name,
            const label initialCapacity,
            IOobjectOption::writeOption wOpt = IOobjectOption::NO_WRITE
        );

        //- Copy construct (no-read) from existing set
        cellSet
        (
            const polyMesh& mesh,
            const word& name,
            const topoSet& set,
            IOobjectOption::writeOption wOpt = IOobjectOption::NO_WRITE
        );

        //- Copy construct (no-read) from labelHashSet,
        //- with search for IOobject instance.
        cellSet
        (
            const polyMesh& mesh,
            const word& name,
            const labelHashSet& labels,
            IOobjectOption::writeOption wOpt = IOobjectOption::NO_WRITE
        );

        //- Move construct (no-read) from labelHashSet,
        //- with search for IOobject instance.
        cellSet
        (
            const polyMesh& mesh,
            const word& name,
            labelHashSet&& labels,
            IOobjectOption::writeOption wOpt = IOobjectOption::NO_WRITE
        );

        //- Copy construct (no-read) from labels,
        //- with search for IOobject instance.
        cellSet
        (
            const polyMesh& mesh,
            const word& name,
            const labelUList& labels,
            IOobjectOption::writeOption wOpt = IOobjectOption::NO_WRITE
        );


        //- Construct from objectRegistry and name.
        //- Used for tetMesh cellSet only.
        cellSet
        (
            const Time&,
            const word& name,
            IOobjectOption::readOption rOpt = IOobjectOption::MUST_READ,
            IOobjectOption::writeOption wOpt = IOobjectOption::NO_WRITE
        );

        //- Construct empty (no-read) with initial labelHashSet capacity
        //- from objectRegistry.
        //- Used for tetMesh cellSet only.
        cellSet
        (
            const Time&,
            const word& name,
            const label initialCapacity,
            IOobjectOption::writeOption wOpt = IOobjectOption::NO_WRITE
        );

        //- Copy construct (no-read) from labelHashSet.
        //- Used for tetMesh cellSet only.
        cellSet
        (
            const Time&,
            const word& name,
            const labelHashSet& labels,
            IOobjectOption::writeOption wOpt = IOobjectOption::NO_WRITE
        );


    // Factory Methods

        //- Read and return contents. Intermediate IOobject is not registered
        static labelHashSet readContents
        (
            const polyMesh& mesh,
            const word& name
        );


    //- Destructor
    virtual ~cellSet() = default;


    // Member functions

        //- Sync cellSet across coupled patches.
        virtual void sync(const polyMesh& mesh)
        {}

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
            const primitiveMesh&,
            const label maxLen
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
