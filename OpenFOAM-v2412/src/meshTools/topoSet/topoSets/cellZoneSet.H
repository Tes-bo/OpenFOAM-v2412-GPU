/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::cellZoneSet

Description
    Like cellSet but -reads data from cellZone -updates cellZone when writing.

SourceFiles
    cellZoneSet.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_cellZoneSet_H
#define Foam_cellZoneSet_H

#include "cellSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class cellZoneSet Declaration
\*---------------------------------------------------------------------------*/

class cellZoneSet
:
    public cellSet
{
    // Private Data

        const polyMesh& mesh_;

        labelList addressing_;

public:

    //- Runtime type information
    TypeNameNoDebug("cellZoneSet");


    // Constructors

        //- Construct from objectRegistry and name
        cellZoneSet
        (
            const polyMesh& mesh,
            const word& name,
            IOobjectOption::readOption rOpt = IOobjectOption::MUST_READ,
            IOobjectOption::writeOption wOpt = IOobjectOption::NO_WRITE
        );

        //- Construct empty (no-read) with initial labelHashSet capacity.
        cellZoneSet
        (
            const polyMesh& mesh,
            const word& name,
            const label initialCapacity,
            IOobjectOption::writeOption wOpt = IOobjectOption::NO_WRITE
        );

        //- Copy construct (no-read) from existing set
        cellZoneSet
        (
            const polyMesh& mesh,
            const word& name,
            const topoSet& set,
            IOobjectOption::writeOption wOpt = IOobjectOption::NO_WRITE
        );



    //- Destructor
    virtual ~cellZoneSet() = default;


    // Member functions

        const labelList& addressing() const noexcept
        {
            return addressing_;
        }

        labelList& addressing() noexcept
        {
            return addressing_;
        }


        //- Sort addressing and make cellSet part consistent with addressing
        void updateSet();

        //- Invert contents.
        //  Insert all members [0,maxLen) which were not in set
        virtual void invert(const label maxLen);

        //- Subset contents. Only elements present in both sets remain.
        virtual void subset(const labelUList& elems);

        //- Subset contents. Only elements present in both sets remain.
        virtual void subset(const topoSet& set);

        //- Add given elements to the set
        virtual void addSet(const labelUList& elems);

        //- Add given elements to the set
        virtual void addSet(const topoSet& set);

        //- Subtract given elements from the set
        virtual void subtractSet(const labelUList& elems);

        //- Subtract given elements from the set
        virtual void subtractSet(const topoSet& set);

        //- Sync cellSet across coupled patches; update cellZone from cellSet
        virtual void sync(const polyMesh& mesh);

        //- Write maxLen items with label and coordinates.
        virtual void writeDebug
        (
            Ostream& os,
            const primitiveMesh&,
            const label maxLen
        ) const;

        //- Write cellZone using stream options
        virtual bool writeObject
        (
            IOstreamOption streamOpt,
            const bool writeOnProc
        ) const;

        //- Update any stored data for new labels
        virtual void updateMesh(const mapPolyMesh& morphMap);

        //- Return max index+1.
        virtual label maxSize(const polyMesh& mesh) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
