/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2020-2021 OpenCFD Ltd.
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
    Foam::functionObjects::div

Group
    grpFieldFunctionObjects

Description
    Computes the divergence of an input field.

    Operands:
    \table
      Operand          | Type                 | Location
      input  | {surfaceScalar,volVector}Field | $FOAM_CASE/\<time\>/\<inpField\>
      output file      | -                    | -
      output field     | volScalarField       | $FOAM_CASE/\<time\>/\<outField\>
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    div1
    {
        // Mandatory entries (unmodifiable)
        type            div;
        libs            (fieldFunctionObjects);

        // Mandatory (inherited) entry (runtime modifiable)
        field           \<field\>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Reqd  | Deflt
      type         | Type name: div                     | word |  yes  | -
      libs         | Library name: fieldFunctionObjects | word |  yes  | -
      field        | Name of the operand field          | word |  yes  | -
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link fieldExpression.H \endlink
     - \link zoneSubSet.H \endlink

    Minimal example by using the \c postProcess utility:
    \verbatim
        postProcess -func "div(\<field\>)"
    \endverbatim

Note
  - To execute \c div function object on an input \<field\>, a numerical scheme
    should be defined for \c div(\<field\>) in \c system/fvSchemes.divSchemes.
  - Optionally the user can specify \c cellZones to create a sub-mesh for the
    \c div calculation.

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - Foam::functionObjects::fieldExpression
    - ExtendedCodeGuide::functionObjects::field::div

SourceFiles
    div.C
    divTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_div_H
#define functionObjects_div_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class div Declaration
\*---------------------------------------------------------------------------*/

class div
:
    public fieldExpression
{
    // Private Member Functions

        //- Calculate the divergence of either a
        //- volScalarField or a surfaceScalarField and register the result
        template<class FieldType>
        bool calcDiv();

        //- Calculate the divergence field and return true if successful
        virtual bool calc();

        //- Helper function for writing submesh fields
        template<class Type>
        bool writeField();


public:

    //- Runtime type information
    TypeName("div");


    // Constructors

        //- Construct from Time and dictionary
        div
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        div(const div&) = delete;

        //- No copy assignment
        void operator=(const div&) = delete;


    //- Destructor
    virtual ~div() = default;


    // Member Functions

        //- Write the result field
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "divTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
