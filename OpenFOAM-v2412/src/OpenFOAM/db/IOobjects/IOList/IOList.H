/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::IOList

Description
    A List of objects of type \<T\> with automated input and output.

SourceFiles
    IOList.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_IOList_H
#define Foam_IOList_H

#include "List.H"
#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class IOList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class IOList
:
    public regIOobject,
    public List<T>
{
    // Private Member Functions

        //- Read if IOobject flags set. Return true if read.
        bool readIOcontents();

public:

    //- The underlying content type
    typedef List<T> content_type;

    //- Runtime type information
    TypeName("List");


    // Constructors

        //- Default copy construct
        IOList(const IOList&) = default;

        //- Construct from IOobject
        explicit IOList(const IOobject& io);

        //- Construct from IOobject and zero size (if not read)
        IOList(const IOobject& io, Foam::zero);

        //- Construct from IOobject and list size (if not read)
        IOList(const IOobject& io, const label len);

        //- Construct from IOobject and a copy of content
        IOList(const IOobject& io, const UList<T>& content);

        //- Construct by transferring the content
        IOList(const IOobject& io, List<T>&& content);


    // Factory Methods

        //- Read and return contents. The IOobject is never registered
        static List<T> readContents(const IOobject& io);


    //- Destructor
    virtual ~IOList() = default;


    // Member Functions

        //- The writeData method for regIOobject write operation
        virtual bool writeData(Ostream& os) const;


    // Member Operators

        //- Copy assignment of entries
        void operator=(const IOList<T>& rhs);

        //- Copy or move assignment of entries
        using List<T>::operator=;
};


/*---------------------------------------------------------------------------*\
                          Class IOListRef Declaration
\*---------------------------------------------------------------------------*/

//- A IOList wrapper for writing external data.
template<class T>
class IOListRef
:
    public regIOobject
{
    // Private Data

        //- Reference to the external content
        refPtr<List<T>> contentRef_;


public:

    //- The underlying content type
    typedef List<T> content_type;


    //- Type is identical to IOList
    virtual const word& type() const
    {
        return IOList<T>::typeName;
    }


    // Generated Methods

        //- No default construct
        IOListRef() = delete;

        //- No copy construct
        IOListRef(const IOListRef&) = delete;

        //- No copy assignment
        void operator=(const IOListRef&) = delete;


    // Constructors

        //- Construct from IOobject and const data reference
        IOListRef(const IOobject& io, const List<T>& content);


    //- Destructor
    virtual ~IOListRef() = default;


    // Member Functions

        //- Allow cast to const content
        //  Fatal if content is not set
        operator const List<T>&() const
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
    #include "IOList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
