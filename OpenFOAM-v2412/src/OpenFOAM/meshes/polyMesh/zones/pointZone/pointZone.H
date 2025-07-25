/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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
    Foam::pointZone

Description
    A subset of mesh points.

    The labels of points in the zone can be obtained from the addressing()
    list.

    For quick check whether a point belongs to the zone use the lookup
    mechanism in pointZoneMesh, where all the zoned points are registered
    with their zone number.

SourceFiles
    pointZone.C
    pointZoneNew.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_pointZone_H
#define Foam_pointZone_H

#include "zone.H"
#include "pointZoneMeshFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class pointZone;
Ostream& operator<<(Ostream& os, const pointZone& zn);


/*---------------------------------------------------------------------------*\
                          Class pointZone Declaration
\*---------------------------------------------------------------------------*/

class pointZone
:
    public zone
{
    // Private Data

        //- Reference to zone list
        const pointZoneMesh& zoneMesh_;

public:

    // Static Data Members

        //- The name associated with the zone-labels dictionary entry
        //- ("pointLabels")
        static const char * const labelsName;


    //- Runtime type information
    TypeName("pointZone");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            pointZone,
            dictionary,
            (
                const word& name,
                const dictionary& dict,
                const label index,
                const pointZoneMesh& zm
            ),
            (name, dict, index, zm)
        );


    // Constructors

        //- No copy construct
        pointZone(const pointZone&) = delete;

        //- Construct an empty zone - name="", index=0
        explicit pointZone(const pointZoneMesh& zm);

        //- Construct an empty zone with specified name and index
        pointZone
        (
            const word& name,
            const label index,
            const pointZoneMesh& zm
        );

        //- Construct from components
        pointZone
        (
            const word& name,
            const labelUList& addr,
            const label index,
            const pointZoneMesh& zm
        );

        //- Construct from components, transferring addressing
        pointZone
        (
            const word& name,
            labelList&& addr,
            const label index,
            const pointZoneMesh& zm
        );

        //- Construct from dictionary
        pointZone
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const pointZoneMesh& zm
        );

        //- Construct empty with original zone information (name, index, groups)
        //- and mesh reference. Optionally specify a new index.
        pointZone
        (
            const pointZone& originalZone,
            const Foam::zero,
            const pointZoneMesh& zm,
            const label newIndex = -1
        );

        //- Construct empty with original zone information (name, groups),
        //- resetting the index and zone mesh reference.
        pointZone
        (
            const pointZone& originalZone,
            const Foam::zero,
            const label index,
            const pointZoneMesh& zm
        );

        //- Construct with original zone information (name, groups),
        //- resetting the point addressing, the index and zone mesh reference.
        pointZone
        (
            const pointZone& originalZone,
            const labelUList& addr,
            const label index,
            const pointZoneMesh& zm
        );

        //- Construct with a new index and zone mesh information, the name
        //- of the original zone, (move) resetting the point addressing.
        pointZone
        (
            const pointZone& originalZone,
            labelList&& addr,
            const label index,
            const pointZoneMesh& zm
        );


        //- Construct and return a clone, resetting the zone mesh
        virtual autoPtr<pointZone> clone(const pointZoneMesh& zm) const
        {
            return autoPtr<pointZone>::New(*this, *this, index(), zm);
        }

        //- Construct and return a clone, resetting the point list
        //- and zone mesh
        virtual autoPtr<pointZone> clone
        (
            const pointZoneMesh& zm,
            const label index,
            const labelUList& addr
        ) const
        {
            return autoPtr<pointZone>::New(*this, addr, index, zm);
        }


    // Selectors

        //- Return a pointer to a new point zone
        //  created on freestore from dictionary
        static autoPtr<pointZone> New
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const pointZoneMesh& zm
        );


    //- Destructor
    virtual ~pointZone() = default;


    // Member Functions

        //- Return reference to the zone mesh
        const pointZoneMesh& zoneMesh() const noexcept { return zoneMesh_; }

        //- The addressing (point IDs) used for the zone
        const labelList& addressing() const noexcept
        {
            return static_cast<const labelList&>(*this);
        }

        //- Helper function to re-direct to zone::localID(...)
        label whichPoint(const label globalPointID) const;


        //- Check zone definition. Return true if in error.
        virtual bool checkDefinition(const bool report = false) const;

        //- Check whether zone is synchronised across coupled boundaries.
        //  \return True if any errors.
        virtual bool checkParallelSync(const bool report = false) const;

        //- Correct patch after moving points
        virtual void movePoints(const pointField&)
        {}

        //- Write dictionary
        virtual void writeDict(Ostream& os) const;


    // Assign addressing

        //- Move reset addressing from another zone.
        virtual void resetAddressing(pointZone&& zn);

        //- Copy reset addressing from another zone
        virtual void resetAddressing(const pointZone& zn);

        //- Copy assign addressing
        virtual void resetAddressing(const labelUList& addr);

        //- Move assign addressing
        virtual void resetAddressing(labelList&& addr);

        //- Assign addressing, clearing demand-driven data
        void operator=(const pointZone& zn);

        //- Assign addressing, clearing demand-driven data
        void operator=(const labelUList& addr);

        //- Assign addressing, clearing demand-driven data
        void operator=(labelList&& addr);


    // I-O

        //- Ostream Operator
        friend Ostream& operator<<(Ostream& os, const pointZone& zn);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
