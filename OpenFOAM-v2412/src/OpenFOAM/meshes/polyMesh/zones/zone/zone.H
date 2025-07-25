/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::zone

Description
    Base class for mesh zones.

    A zone is a list of labels (eg, cells, faces, points) with
    a name and associated with an index within another list.

SourceFiles
    zone.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_zone_H
#define Foam_zone_H

#include "zoneIdentifier.H"
#include "labelList.H"
#include "typeInfo.H"
#include "Map.H"
#include "pointFieldFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class zone;
class dictionary;
Ostream& operator<<(Ostream& os, const zone& zn);

/*---------------------------------------------------------------------------*\
                           Class zone Declaration
\*---------------------------------------------------------------------------*/

class zone
:
    public zoneIdentifier,
    public labelList
{
    // Private Data

        //- Demand-driven: map of labels in zone for fast location lookup
        mutable std::unique_ptr<Map<label>> lookupMapPtr_;


public:

    //- Runtime type information
    TypeNameNoDebug("zone");


    // Constructors

        //- Default construct
        zone();

        //- Construct an empty zone with specified name and index
        zone(const word& name, const label index);

        //- Copy construct from components
        zone
        (
            const word& name,
            const labelUList& addr,
            const label index
        );

        //- Move construct from components
        zone
        (
            const word& name,
            labelList&& addr,
            const label index
        );

        //- Construct from dictionary
        zone
        (
            const word& name,
            const dictionary& dict,
            const word& labelsName,
            const label index
        );

        //- Construct given the name of the original zone (name is used)
        //- with new addressing and index (-ve: retain zone index).
        zone
        (
            const zone& originalZone,
            const labelUList& addr,
            const label index
        );

        //- Construct given the name of the original zone (name is used)
        //- and (move) resetting addressing and index (-ve: retain zone index).
        zone
        (
            const zone& originalZone,
            labelList&& addr,
            const label index
        );


    //- Destructor
    virtual ~zone() = default;


    // Member Functions

        //- Demand-driven: the look-up map from global to local id
        const Map<label>& lookupMap() const;

        //- Lookup local address in zone for given global index.
        //  \return the local address, or -1 if the item is not in the zone
        label localID(const label globalID) const;

        //- The addressing used by the zone
        const labelList& addressing() const noexcept
        {
            return static_cast<const labelList&>(*this);
        }

        //- Clear addressing
        //- (remove lookup maps and other auxiliary information)
        virtual void clearAddressing();

        //- Clear primitive addressing
        virtual void clearPrimitives();

        //- Check zone definition. Return true if in error.
        virtual bool checkDefinition(const bool report = false) const = 0;

        //- Check zone definition with max size given. Return true if in error.
        virtual bool checkDefinition
        (
            const label maxSize,
            const bool report = false
        ) const;

        //- Correct patch after moving points
        virtual void movePoints(const pointField& pts)
        {}

        //- Write
        virtual void write(Ostream& os) const;

        //- Write dictionary
        virtual void writeDict(Ostream& os) const = 0;


    // I-O

        //- Ostream Operator
        friend Ostream& operator<<(Ostream& os, const zone& zn);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
