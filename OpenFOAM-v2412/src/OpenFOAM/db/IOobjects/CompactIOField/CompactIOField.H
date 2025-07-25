/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020-2024 OpenCFD Ltd.
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
    Foam::CompactIOField

Description
    A Field of objects of type \<T\> with automated input and output using
    a compact storage. Behaves like IOField except when binary output in
    case it writes a CompactListList.

    Useful for fields of small subfields e.g. in lagrangian

SourceFiles
    CompactIOField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_CompactIOField_H
#define Foam_CompactIOField_H

#include "IOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class T, class BaseType> class CompactIOField;

template<class T, class BaseType> Istream& operator>>
(
    Istream&,
    CompactIOField<T, BaseType>&
);

template<class T, class BaseType> Ostream& operator<<
(
    Ostream&,
    const CompactIOField<T, BaseType>&
);

/*---------------------------------------------------------------------------*\
                        Class CompactIOField Declaration
\*---------------------------------------------------------------------------*/

template<class T, class BaseType>
class CompactIOField
:
    public regIOobject,
    public Field<T>
{
    // Private Member Functions

        //- Read if IOobject flags set and 'on-proc' is true.
        //- Reads according to the header type.
        //  Return true if read (only accurate when readOnProc == true).
        bool readIOcontents(bool readOnProc = true);

public:

    //- Runtime type information
    TypeName("FieldField");


    // Constructors

        //- Default copy construct
        CompactIOField(const CompactIOField&) = default;

        //- Construct from IOobject
        explicit CompactIOField(const IOobject& io);

        //- Construct from IOobject, with local processor conditional reading
        CompactIOField(const IOobject& io, const bool readOnProc);

        //- Construct from IOobject and zero size (if not read)
        CompactIOField(const IOobject& io, Foam::zero);

        //- Construct from IOobject and size (if not read)
        CompactIOField(const IOobject& io, const label len);

        //- Construct from IOobject and a List/Field content
        CompactIOField(const IOobject& io, const UList<T>& content);

        //- Construct by transferring the Field contents
        CompactIOField(const IOobject& io, Field<T>&& content);


    //- Destructor
    virtual ~CompactIOField() = default;


    // Member Functions

        //- Write using stream options
        virtual bool writeObject
        (
            IOstreamOption streamOpt,
            const bool writeOnProc
        ) const;

        virtual bool writeData(Ostream& os) const;


    // Member Operators

        //- Copy assignment of entries
        void operator=(const CompactIOField<T, BaseType>& rhs);

        //- Copy or move assignment of entries
        using Field<T>::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CompactIOField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
