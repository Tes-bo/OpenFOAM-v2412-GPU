/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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
    Foam::functionObjects::components

Group
    grpFieldFunctionObjects

Description
    Extracts the components of elements of a field and outputs the result
    into new fields, where the fields possess the same sizes and dimensions.

    Depending on the base type of the input field, the number of output fields
    and output suffixes are as follows (the base type cannot be scalar since
    scalars have no components):

    \table
      Base type         | Number of components | Suffixes
      vector            | 3                    | (x y z)
      sphericalTensor   | 3                    | (x y z)
      symmTensor        | 6                    | (xx xy xz yy yz zz)
      tensor            | 9                    | (xx xy xz yx yy yz zx zy zz)
    \endtable

    Operands:
    \table
      Operand           | Type                                        | Location
      input             | {vol,surface}\<Type\>Field (Type != Scalar) <!--
                    --> | $FOAM_CASE/\<time\>/\<inpField\>
      output file       | -                                           | -
      output field      | {vol,surface}ScalarField for each component <!--
                    --> | $FOAM_CASE/\<time\>/\<outField\>\<suffix\>
    \endtable

    where \c \<Type\>=Scalar/Vector/SphericalTensor/SymmTensor/Tensor.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    components1
    {
        // Mandatory entries (unmodifiable)
        type            components;
        libs            (fieldFunctionObjects);

        // Mandatory (inherited) entry (runtime modifiable)
        field           <field>;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Req'd | Dflt
      type         | Type name: components              | word |  yes  | -
      libs         | Library name: fieldFunctionObjects | word |  yes  | -
      field        | Name of the operand field          | word |  yes  | -
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link fieldExpression.H \endlink

    Minimal example by using the \c postProcess utility:
    \verbatim
        postProcess -func "components(<field>)"
    \endverbatim

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - ExtendedCodeGuide::functionObjects::field::components

SourceFiles
    components.C
    componentsTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_components_H
#define Foam_functionObjects_components_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class components Declaration
\*---------------------------------------------------------------------------*/

class components
:
    public fieldExpression
{
    // Private Member Data

        //- List of the component field names
        wordList resultNames_;


    // Private Member Functions

        //- Calculate the components of the field with the specified type
        //- and register the result
        template<class GeoFieldType>
        bool calcComponents(const GeoFieldType& field);

        //- Calculate the components of the field with the specified
        //- element type and register the result
        template<class Type>
        bool calcComponents();

        //- Calculate the components of the field
        //  \return True if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("components");


    // Constructors

        //- Construct from Time and dictionary
        components
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        components(const components&) = delete;

        //- No copy assignment
        void operator=(const components&) = delete;


    //- Destructor
    virtual ~components() = default;


    // Member Functions

        //- Write the component fields
        virtual bool write();

        //- Clear the component fields from the objectRegistry
        virtual bool clear();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
