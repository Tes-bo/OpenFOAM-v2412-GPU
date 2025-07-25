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
    Foam::pointSet

Description
    A set of point labels.

SourceFiles
    pointSet.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_pointSet_H
#define Foam_pointSet_H

#include "topoSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class pointSet Declaration
\*---------------------------------------------------------------------------*/

class pointSet
:
    public topoSet
{
public:

    //- Runtime type information
    TypeNameNoDebug("pointSet");


    // Constructors

        //- Construct from IOobject
        explicit pointSet(const IOobject& io);

        //- Construct empty (no-read) with IOobject information
        pointSet(const IOobject& io, const Foam::zero);

        //- Construct from polyMesh (registry) and name
        pointSet
        (
            const polyMesh& mesh,
            const word& name,
            IOobjectOption::readOption rOpt = IOobjectOption::MUST_READ,
            IOobjectOption::writeOption wOpt = IOobjectOption::NO_WRITE,
            IOobjectOption::registerOption reg = IOobjectOption::LEGACY_REGISTER
        );

        //- Construct from polyMesh (registry), name and registration option
        pointSet
        (
            const polyMesh& mesh,
            const word& name,
            IOobjectOption::registerOption reg,
            IOobjectOption::readOption rOpt = IOobjectOption::MUST_READ,
            IOobjectOption::writeOption wOpt = IOobjectOption::NO_WRITE
        )
        :
            pointSet(mesh, name, rOpt, wOpt, reg)
        {}

        //- Construct empty (no-read) with initial labelHashSet capacity.
        pointSet
        (
            const polyMesh& mesh,
            const word& name,
            const label initialCapacity,
            IOobjectOption::writeOption wOpt = IOobjectOption::NO_WRITE
        );

        //- Copy construct (no-read) from existing set
        pointSet
        (
            const polyMesh& mesh,
            const word& name,
            const topoSet& set,
            IOobjectOption::writeOption wOpt = IOobjectOption::NO_WRITE
        );

        //- Copy construct (no-read) from labelHashSet,
        //- with search for IOobject instance.
        pointSet
        (
            const polyMesh& mesh,
            const word& name,
            const labelHashSet& labels,
            IOobjectOption::writeOption wOpt = IOobjectOption::NO_WRITE
        );

        //- Move construct (no-read) from labelHashSet,
        //- with search for IOobject instance.
        pointSet
        (
            const polyMesh& mesh,
            const word& name,
            labelHashSet&& labels,
            IOobjectOption::writeOption wOpt = IOobjectOption::NO_WRITE
        );

        //- Copy construct (no-read) from labels,
        //- with search for IOobject instance.
        pointSet
        (
            const polyMesh& mesh,
            const word& name,
            const labelUList& labels,
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
    virtual ~pointSet() = default;


    // Member Functions

        //- Sync set across coupled patches. Adds coupled points to set.
        virtual void sync(const polyMesh& mesh);

        //- Return max index+1.
        virtual label maxSize(const polyMesh& mesh) const;

        //- Update any stored data for new labels
        virtual void updateMesh(const mapPolyMesh& morphMap);

        //- Update any stored data for mesh redistribution.
        virtual void distribute(const mapDistributePolyMesh&);

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
