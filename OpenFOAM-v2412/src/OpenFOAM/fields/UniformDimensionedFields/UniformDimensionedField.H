/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::UniformDimensionedField

Description
    Dimensioned<Type> registered with the database as a registered IOobject
    which has the functionality of a uniform field and allows values from the
    top-level code to be passed to boundary conditions etc.

    Is a 'global' field, same on all processors

SourceFiles
    UniformDimensionedField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_UniformDimensionedField_H
#define Foam_UniformDimensionedField_H

#include "regIOobject.H"
#include "dimensionedType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class UniformDimensionedField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class UniformDimensionedField
:
    public regIOobject,
    public dimensioned<Type>
{
public:

    //- Runtime type information
    TypeName("UniformDimensionedField");


    // Constructors

        //- Construct from components. Either reads or uses supplied value.
        //  The name of the dimensioned\<Type\> defines the name of the
        //  field (if non-empty)
        UniformDimensionedField
        (
            const IOobject& io,
            const dimensioned<Type>& dt
        );

        //- Construct from components. Either reads or uses supplied value.
        //  The name of the IOobject defines the name of the field.
        //  Parameter ordering as per dimensioned\<Type\>
        UniformDimensionedField
        (
            const IOobject& io,
            const dimensionSet& dims,
            const Type& val
        );

        //- Construct from components. Either reads or uses supplied value.
        //  The name of the IOobject defines the name of the field.
        //  Parameter ordering as per DimensionedField\<Type\>
        UniformDimensionedField
        (
            const IOobject& io,
            const Type& val,
            const dimensionSet& dims
        );

        //- Construct as copy
        UniformDimensionedField(const UniformDimensionedField<Type>&);

        //- Construct from IOobject. Either reads or sets dimensionless zero
        explicit UniformDimensionedField(const IOobject& io);


    //- Destructor
    virtual ~UniformDimensionedField() = default;


    // Member Functions

        //- Use name from dimensioned\<Type\>, not from regIOobject
        virtual const word& name() const
        {
            return dimensioned<Type>::name();
        }

        //- The readData method required for regIOobject read operation
        virtual bool readData(Istream&);

        //- The writeData method required for regIOobject write operation
        bool writeData(Ostream&) const;

        //- Is object global
        virtual bool global() const
        {
            return true;
        }

        //- Return complete path + object name if the file exists
        //- either in the case/processor or case otherwise null
        virtual fileName filePath() const
        {
            return globalFilePath(type());
        }


    // Member Operators

        //- Assign name, dimensions and value of the dimensioned\<Type\>
        void operator=(const UniformDimensionedField<Type>& rhs);

        //- Assign name, dimensions and value to the dimensioned\<Type\>
        void operator=(const dimensioned<Type>& rhs);

        //- Return value
        const Type& operator[](const label) const noexcept
        {
            return dimensioned<Type>::value();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Global file type for UniformDimensionedField
template<class Type>
struct is_globalIOobject<UniformDimensionedField<Type>> : std::true_type {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "UniformDimensionedField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
