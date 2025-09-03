/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2023 OpenCFD Ltd.
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
    Foam::rawIOField

Description
    Like IOField but falls back to raw IFstream if no header found.
    Optionally reads average value. For use in MappedFile container.

SourceFiles
    rawIOField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_rawIOField_H
#define Foam_rawIOField_H

#include "IOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class rawIOField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class rawIOField
:
    public regIOobject,
    public Field<Type>
{
    // Private Data

        //- Has an average value
        bool hasAverage_;

        //- The average for the field (zero if not used)
        Type average_;


    // Private Member Functions

        //- Read contents and average
        void readIOcontents(Istream&, IOobjectOption::readOption readAverage);

        //- Read if IOobject flags set. Return true if read.
        bool readIOcontents(IOobjectOption::readOption readAverage);


public:

    //- The underlying content type
    typedef Field<Type> content_type;

    //- Runtime type information
    TypeName("rawField");


    // Constructors

        //- Default copy construct
        rawIOField(const rawIOField&) = default;

        //- Construct from IOobject and specified request for reading average
        //- (default: NO_READ)
        explicit rawIOField
        (
            const IOobject& io,
            IOobjectOption::readOption readAverage = IOobjectOption::NO_READ
        );

        //- Construct from IOobject and specified request for reading average
        rawIOField(const IOobject& io, const bool readAverage);


    // Factory Methods

        //- Read and return contents. The IOobject will not be registered
        static Field<Type> readContents(const IOobject& io);


    //- Destructor
    virtual ~rawIOField() = default;


    // Member Functions

        //- The writeData method for regIOobject write operation
        virtual bool writeData(Ostream& os) const;

        //- The field content
        const Field<Type>& field() const noexcept
        {
            return static_cast<const Field<Type>&>(*this);
        }

        //- The field content, non-const access
        Field<Type>& field() noexcept
        {
            return static_cast<Field<Type>&>(*this);
        }

        //- The average value (if any)
        const Type& average() const noexcept
        {
            return average_;
        }

        //- Has an average value
        bool hasAverage() const noexcept
        {
            return hasAverage_;
        }

        //- Set an average value
        void setAverage(const Type& val)
        {
            hasAverage_ = true;
            average_ = val;
        }

        //- Clear average value
        void clearAverage()
        {
            hasAverage_ = false;
            average_ = Zero;
        }


    // Member Operators

        //- Copy or move assignment of entries
        using Field<Type>::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "rawIOField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
