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
    Foam::IOField

Description
    A primitive field of type \<T\> with automated input and output.

SourceFiles
    IOField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_IOField_H
#define Foam_IOField_H

#include "Field.H"
#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class IOField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class IOField
:
    public regIOobject,
    public Field<Type>
{
    // Private Member Functions

        //- Read if IOobject flags set and 'on-proc' is true.
        //  Return true if read (only accurate when readOnProc == true).
        bool readIOcontents(bool readOnProc = true);

public:

    //- The underlying content type
    typedef Field<Type> content_type;

    //- Runtime type information
    TypeName("Field");


    // Constructors

        //- Default copy construct
        IOField(const IOField&) = default;

        //- Construct from IOobject
        explicit IOField(const IOobject& io);

        //- Construct from IOobject, with local processor conditional reading
        IOField(const IOobject& io, const bool readOnProc);

        //- Construct from IOobject and zero size (if not read)
        IOField(const IOobject& io, Foam::zero);

        //- Construct from IOobject and field size (if not read)
        IOField(const IOobject& io, const label len);

        //- Construct from IOobject and copy of List/Field content
        IOField(const IOobject& io, const UList<Type>& content);

        //- Construct by transferring the Field content
        IOField(const IOobject& io, Field<Type>&& content);

        //- Construct by copying/moving tmp content
        IOField(const IOobject& io, const tmp<Field<Type>>& tfld);


    // Factory Methods

        //- Read and return contents. The IOobject will not be registered
        static Field<Type> readContents(const IOobject& io);


    //- Destructor
    virtual ~IOField() = default;


    // Member Functions

        //- The writeData method for regIOobject write operation
        virtual bool writeData(Ostream& os) const;


    // Member Operators

        //- Copy assignment of entries
        void operator=(const IOField<Type>& rhs);

        //- Copy or move assignment of entries
        using Field<Type>::operator=;
};


/*---------------------------------------------------------------------------*\
                         Class IOFieldRef Declaration
\*---------------------------------------------------------------------------*/

//- A IOField wrapper for writing external data.
template<class Type>
class IOFieldRef
:
    public regIOobject
{
    // Private Data

        //- Reference to the external content
        refPtr<Field<Type>> contentRef_;


public:

    //- The underlying content type
    typedef Field<Type> content_type;


    //- Type is identical to IOField
    virtual const word& type() const
    {
        return IOField<Type>::typeName;
    }


    // Generated Methods

        //- No default construct
        IOFieldRef() = delete;

        //- No copy construct
        IOFieldRef(const IOFieldRef&) = delete;

        //- No copy assignment
        void operator=(const IOFieldRef&) = delete;


    // Constructors

        //- Construct from IOobject and const data reference
        IOFieldRef(const IOobject& io, const Field<Type>& content);


    //- Destructor
    virtual ~IOFieldRef() = default;


    // Member Functions

        //- Allow cast to const content
        //  Fatal if content is not set
        operator const Field<Type>&() const
        {
            return contentRef_.cref();
        }

        //- The writeData method for regIOobject write operation
        //  Fatal if content is not set
        virtual bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "IOField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
