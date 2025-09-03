/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
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
    Foam::surfZone

Description
    A surface zone on a MeshedSurface.

    Similar in concept to a faceZone, but the face list is contiguous.

SourceFiles
    surfZone.C

\*---------------------------------------------------------------------------*/

#ifndef surfZone_H
#define surfZone_H

#include "surfZoneIdentifier.H"
#include "labelRange.H"
#include "autoPtr.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class surfZone Declaration
\*---------------------------------------------------------------------------*/

class surfZone
:
    public surfZoneIdentifier
{
    // Private Data

        //- Size of this group in the face list
        label size_;

        //- Start label of this group in the face list
        label start_;


public:

    // Generated Methods

        //- Copy construct
        surfZone(const surfZone&) = default;

        //- Copy assignment
        surfZone& operator=(const surfZone&) = default;


    // Constructors

        //- Default construct, with zero start, size, index
        surfZone();

        //- Construct with name, size. With zero start, index
        surfZone(const word& name, const label size);

        //- Construct from components
        surfZone
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const word& geometricType = word::null
        );

        //- Construct from dictionary
        surfZone
        (
            const word& name,
            const dictionary& dict,
            const label index
        );

        //- Copy construct, resetting the index
        surfZone(const surfZone&, const label index);

        //- Return clone
        autoPtr<surfZone> clone() const
        {
            NotImplemented;
            return nullptr;
        }

        static autoPtr<surfZone> New(Istream& is)
        {
            const word name(is);
            const dictionary dict(is);

            return autoPtr<surfZone>::New(name, dict, 0);
        }


    // Member Functions

        //- The start label of this zone in the face list
        label start() const
        {
            return start_;
        }

        //- Modifiable start label of this zone in the face list
        label& start()
        {
            return start_;
        }

        //- The size of this zone in the face list
        label size() const
        {
            return size_;
        }

        //- Modifiable size of this zone in the face list
        label& size()
        {
            return size_;
        }

        //- The start/size range of this zone
        labelRange range() const
        {
            return labelRange(start_, size_);
        }

        //- Write dictionary, includes surrounding braces
        //  \warning Prior to 2020-01 was identical to operator<< output
        void write(Ostream& os) const;


    // Housekeeping

        //- Removed(2020-01) Construct from Istream
        //  \deprecated(2020-01) - unused, inconsistent
        surfZone(Istream& is, const label index) = delete;

        //- Deprecated(2020-01) Write dictionary
        //  \deprecated(2020-01) - Write dictionary
        FOAM_DEPRECATED_FOR(2020-01, "write() or operator<<")
        void writeDict(Ostream& os) const
        {
            write(os);
        }
};


// Global Operators

//- Compare surfZones for equality
bool operator==(const surfZone& a, const surfZone& b);

//- Compare surfZones for inequality
bool operator!=(const surfZone& a, const surfZone& b);


//- Read as dictionary
Istream& operator>>(Istream& is, surfZone& obj);

//- Write as dictionary
Ostream& operator<<(Ostream&, const surfZone& obj);


// Global Functions

//- The labelRange of a surfZone
template<>
struct labelRangeOp<surfZone>
{
    labelRange operator()(const surfZone& zone) const
    {
        return zone.range();
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
