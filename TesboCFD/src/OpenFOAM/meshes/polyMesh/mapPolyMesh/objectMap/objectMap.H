/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2024 OpenCFD Ltd.
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
    Foam::objectMap

Description
    An objectMap is a pair of labels defining the mapping of an object from
    another object, e.g. a cell mapped from a point.

SourceFiles
    objectMapI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_objectMap_H
#define Foam_objectMap_H

#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class objectMap;
inline bool operator==(const objectMap& a, const objectMap& b);
inline bool operator!=(const objectMap& a, const objectMap& b);
inline Ostream& operator<<(Ostream&, const objectMap&);
inline Istream& operator>>(Istream&, objectMap&);


/*---------------------------------------------------------------------------*\
                          Class objectMap Declaration
\*---------------------------------------------------------------------------*/

class objectMap
{
    // Private Data

        //- Object index
        label index_;

        //- Master object indices
        labelList objects_;


public:

    // Constructors

        //- Default construct, with index=-1 and no objects
        objectMap() noexcept : index_(-1) {}

        //- Copy construct from components
        inline objectMap(const label index, const labelUList& master);

        //- Move construct from components
        inline objectMap(const label index, labelList&& master);

        //- Construct from Istream
        inline explicit objectMap(Istream& is);


    // Member Functions

        //- True if there are no objects
        bool empty() const noexcept { return objects_.empty(); }

        //- The number of objects
        label size() const noexcept { return objects_.size(); }

        //- The object index
        label& index() noexcept { return index_; }

        //- The object index
        label index() const noexcept { return index_; }

        //- The master object indices
        labelList& masterObjects() noexcept { return objects_; }

        //- The master object indices
        const labelList& masterObjects() const noexcept { return objects_; }


    // Friend Operators

        friend bool operator==(const objectMap& a, const objectMap& b);
        friend bool operator!=(const objectMap& a, const objectMap& b);


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const objectMap&);
        friend Istream& operator>>(Istream&, objectMap&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "objectMapI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
