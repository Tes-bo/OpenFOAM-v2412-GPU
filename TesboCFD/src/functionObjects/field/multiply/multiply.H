/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::functionObjects::multiply

Group
    grpFieldFunctionObjects

Description
    Multiplies a given list of (at least two or more) fields and outputs the
    result into a new field.

    \verbatim
        fieldResult = field1 * field2 * ... * fieldN
    \endverbatim

    Operands:
    \table
      Operand      | Type                          | Location
      input        | vol\<Type\>Field(s) <!--
               --> |$FOAM_CASE/\<time\>/\<inpField\>s
      output field | vol\<Type\>Field    <!--
               --> | $FOAM_CASE/\<time\>/\<outField\>
    \endtable

    where \c \<Type\>=Scalar/Vector/SphericalTensor/SymmTensor/Tensor.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    multiply1
    {
        // Mandatory entries (unmodifiable)
        type    multiply;
        libs    (fieldFunctionObjects);

        // Mandatory (inherited) entry (runtime modifiable)
        fields  (<field1> <field2> ... <fieldN>);

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Reqd | Dflt
      type         | Type name: multiply                | word | yes | -
      libs         | Library name: fieldFunctionObjects | word | yes | -
      fields       | Names of the operand fields        | wordList | yes | -
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link fieldsExpression.H \endlink

See also
    - Foam::functionObject
    - Foam::functionObjects::fieldsExpression
    - Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    multiply.C
    multiplyTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_multiply_H
#define functionObjects_multiply_H

#include "fieldsExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class multiply Declaration
\*---------------------------------------------------------------------------*/

class multiply
:
    public fieldsExpression
{
    //- Helper struct to define valid operations
    template<class Type1, class Type2>
    struct is_valid_op
    {
        static constexpr bool value =
            (pTraits<Type1>::rank == 0 || pTraits<Type2>::rank == 0)
         || (pTraits<Type1>::rank == 1 && pTraits<Type2>::rank == 1);
    };


    // Private Member Functions

        template<class Type>
        bool initialiseResult(const word& fieldName);

        template<class Type>
        bool multiplyResult(const word& fieldName, bool& processed);

        template<class Type1, class Type2>
        typename std::enable_if
        <
            is_valid_op<Type1, Type2>::value, bool
        >::type
        multiplyFieldType
        (
            GeometricField<Type1, fvPatchField, volMesh>& result,
            const word& fieldName,
            bool& processed
        );

        template<class Type1, class Type2>
        typename std::enable_if
        <
            !is_valid_op<Type1, Type2>::value, bool
        >::type
        multiplyFieldType
        (
            GeometricField<Type1, fvPatchField, volMesh>& result,
            const word& fieldName,
            bool& processed
        );

        //- Multiply the list of fields and return true if successful
        virtual bool calc();


public:

    friend class fieldsExpression;


    //- Runtime type information
    TypeName("multiply");


    // Constructors

        //- Construct from Time and dictionary
        multiply
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        multiply(const multiply&) = delete;

        //- No copy assignment
        void operator=(const multiply&) = delete;


    //- Destructor
    virtual ~multiply() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "multiplyTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
