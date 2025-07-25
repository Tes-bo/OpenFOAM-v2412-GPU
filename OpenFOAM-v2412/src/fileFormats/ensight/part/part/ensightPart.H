/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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
    Foam::ensightPart

Description
    Base class for ensightCells, ensightFaces, ensightOutputSurfaces.

SourceFiles
    ensightPart.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_ensightPart_H
#define Foam_ensightPart_H

#include "ensightGeoFile.H"
#include "labelList.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class ensightPart Declaration
\*---------------------------------------------------------------------------*/

class ensightPart
{
    // Private Data

        //- Part index within a list.
        //  The ensight part number is typically this value +1.
        label index_;

        //- OpenFOAM identifier (patch index, zone index, etc).
        //  An unused identifier is -1
        label identifier_;

        //- Part name (or description)
        string name_;

        //- Linear list of element ids (face/cell)
        //  Sub-sectioning by element type is done by derived classes
        labelList address_;


protected:

        //- Element addressing
        const labelList& addressing() const noexcept
        {
            return address_;
        }

        //- Element addressing
        labelList& addressing() noexcept
        {
            return address_;
        }

        //- Clear element addressing
        void clear()
        {
            address_.clear();
        }

        //- Increase addressing by specified offset value
        //  Eg, change local to global id
        void incrAddressing(const label off);

        //- Decrease addressing by specified offset value
        //  Eg, change global to local id
        void decrAddressing(const label off);


public:

    //- Declare type-name, virtual type (without debug switch)
    TypeNameNoDebug("ensightPart");


    // Constructors

        //- Default construct. Index=0, identifier = -1
        ensightPart();

        //- Default construct, with description/partName
        explicit ensightPart(const string& description);


    //- Destructor
    virtual ~ensightPart() = default;


    // Member Functions

        //- The index in a list (0-based)
        label index() const noexcept
        {
            return index_;
        }

        //- The index in a list (0-based)
        label& index() noexcept
        {
            return index_;
        }

        //- OpenFOAM identifier (patch, zone, etc), -1 when not in use.
        label identifier() const noexcept
        {
            return identifier_;
        }

        //- OpenFOAM identifier (patch, zone, etc), -1 when not in use.
        label& identifier() noexcept
        {
            return identifier_;
        }

        //- Processor-local test for any elements.
        bool empty() const noexcept
        {
            return address_.empty();
        }

        //- Processor-local size of all elements.
        label size() const noexcept
        {
            return address_.size();
        }

        //- The part name or description
        const string& name() const noexcept
        {
            return name_;
        }

        //- Change the part name or description
        void rename(const string& value)
        {
            name_ = value;
        }

        //- Change the part name or description
        void rename(string&& value)
        {
            name_ = std::move(value);
        }


    // Output

        //- Write information about the object as a dictionary,
        //- optionally write all element addresses
        virtual void writeDict(Ostream& os, const bool full=false) const
        {}


    // Member Operators

        //- Processor-local element id from linear-list of addresses.
        label operator[](const label i) const
        {
            return address_[i];
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Housekeeping

//- Deprecated(2020-02) - use ensightOutput or member write() methods
//  \deprecated(2020-02) - use ensightOutput or member write() methods
void operator<<(ensightGeoFile&, const ensightPart&) = delete;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
