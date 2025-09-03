/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
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
    Foam::functionObjects::Lambda2

Group
    grpFieldFunctionObjects

Description
    Computes the second largest eigenvalue of the sum of the square of
    the symmetrical and anti-symmetrical parts of the velocity gradient tensor.

    Operands:
    \table
      Operand          | Type           | Location
      input            | volVectorField | \<case\>/\<time\>/\<inpField\>
      output file      | -              | -
      output field     | volScalarField | \<case\>/\<time\>/\<outField\>
    \endtable

    References:
    \verbatim
        Governing equation (tag:JH):
            Jeong, J., & Hussain, F. (1995).
            On the identification of a vortex.
            Journal of Fluid Mechanics, 285, 69-94.
            DOI:10.1017/S0022112095000462
    \endverbatim

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    Lambda21
    {
        // Mandatory entries
        type            Lambda2;
        libs            (fieldFunctionObjects);

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Reqd | Deflt
      type         | Type name: Lambda2                 | word | yes  | -
      libs         | Library name: fieldFunctionObjects | word | yes  | -
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link fieldExpression.H \endlink

    Minimal example by using the \c postProcess utility:
    \verbatim
        postProcess -func Lambda2
    \endverbatim

SourceFiles
    Lambda2.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_Lambda2_H
#define functionObjects_Lambda2_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class Lambda2 Declaration
\*---------------------------------------------------------------------------*/

class Lambda2
:
    public fieldExpression
{
    // Private Member Functions

        //- Calculate the Lambda2 field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("Lambda2");


    // Constructors

        //- Construct from Time and dictionary
        Lambda2
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        Lambda2(const Lambda2&) = delete;

        //- No copy assignment
        void operator=(const Lambda2&) = delete;


    //- Destructor
    virtual ~Lambda2() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
