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
    Foam::functionObjects::pow

Group
    grpFieldFunctionObjects

Description
    Computes the power of an input \c volScalarField.

    \f[
        f = s f_0^n + t
    \f]

    where
    \vartable
        f     | Output volScalarField
        f_0   | Input volScalarField
        n     | Exponent
        s     | Scaling factor
        t     | Offset factor
    \endvartable

    \table
      Operand       | Type           | Location
      input         | volScalarField | $FOAM_CASE/\<time\>/\<inpField\>
      output file   | -              | -
      output field  | volScalarField | $FOAM_CASE/\<time\>/\<outField\>
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    pow1
    {
        // Mandatory entries (unmodifiable)
        type            pow;
        libs            (fieldFunctionObjects);

        // Mandatory (inherited) entry (runtime modifiable)
        field           <inpField>;

        // Mandatory entry (runtime modifiable)
        n               0.25;

        // Optional entries (runtime modifiable)
        checkDimensions false;
        scale           1.0;
        offset          0.0;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Req'd | Dflt
      type         | Type name: pow                     | word |  yes  | -
      libs         | Library name: fieldFunctionObjects | word |  yes  | -
      field        | Name of the operand field          | word |  yes  | -
      n            | Exponent                           | scalar | yes | -
      checkDimensions | Flag to check dimensions of the operand field <!--
                                                    --> | bool   | no  | true
      scale        | Scaling factor - \c s above        | scalar | no  | 1.0
      offset       | Offset factor - \c t above         | scalar | no  | 0.0
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link fieldExpression.H \endlink

    Minimal example by using the \c postProcess utility:
    \verbatim
        postProcess -func "pow(<inpField>)" -scale 1.0 -offset 0.0
    \endverbatim

See also
    - Foam::functionObject
    - Foam::functionObjects::fieldExpression
    - Foam::functionObjects::fvMeshFunctionObject
    - ExtendedCodeGuide::functionObjects::field::pow

SourceFiles
    pow.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_pow_H
#define functionObjects_pow_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class pow Declaration
\*---------------------------------------------------------------------------*/

class pow
:
    public fieldExpression
{
    // Private Data

        //- Flag to check dimensions of the operand
        Switch checkDimensions_;

        //- Exponent
        scalar n_;

        //- Scaling factor
        scalar scale_;

        //- Offset factor
        scalar offset_;


    // Private Member Functions

        //- Calculate the pow field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("pow");


    // Constructors

        //- Construct from Time and dictionary
        pow
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        pow(const pow&) = delete;

        //- No copy assignment
        void operator=(const pow&) = delete;


    //- Destructor
    virtual ~pow() = default;


    // Member Functions

        //- Read the randomise data
        virtual bool read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
