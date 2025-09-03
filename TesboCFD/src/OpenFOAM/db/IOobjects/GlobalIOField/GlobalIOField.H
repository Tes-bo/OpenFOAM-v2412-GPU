/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2017 OpenFOAM Foundation
    Copyright (C) 2016-2024 OpenCFD Ltd.
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
    Foam::GlobalIOField

Description
    IOField with global data (so optionally read from master)

SourceFiles
    GlobalIOField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_GlobalIOField_H
#define Foam_GlobalIOField_H

#include "IOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class GlobalIOField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class GlobalIOField
:
    public regIOobject,
    public Field<Type>
{
public:

    //- The underlying content type
    typedef Field<Type> content_type;

    //- Runtime type information
    TypeName("Field");


    // Constructors

        //- Default copy construct
        GlobalIOField(const GlobalIOField&) = default;

        //- Construct from IOobject
        explicit GlobalIOField(const IOobject& io);

        //- Construct from IOobject and field size (if not read)
        GlobalIOField(const IOobject& io, const label len);

        //- Construct from IOobject and copy of List/Field content
        GlobalIOField(const IOobject& io, const UList<Type>& content);

        //- Construct by transferring the Field content
        GlobalIOField(const IOobject& io, Field<Type>&& content);

        //- Construct by copying/moving tmp content
        GlobalIOField(const IOobject& io, const tmp<Field<Type>>& tfld);


    // Factory Methods

        //- Read and return contents. The IOobject is never registered
        static Field<Type> readContents(const IOobject& io);


    //- Destructor
    virtual ~GlobalIOField() = default;


    // Member Functions

        //- The object is global
        virtual bool global() const
        {
            return true;
        }

        //- Return complete path + object name if the file exists
        //- either in the case/processor or case, otherwise null
        virtual fileName filePath() const
        {
            return globalFilePath(type());
        }

        //- The readData method for regIOobject read operation
        virtual bool readData(Istream& is);

        //- The writeData method for regIOobject write operation
        bool writeData(Ostream& os) const;


    // Member Operators

        //- Copy assignment of entries
        void operator=(const GlobalIOField<Type>& rhs);

        //- Copy or move assignment of entries
        using Field<Type>::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "GlobalIOField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
