/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    Foam::functionObjects::LambVector

Group
    grpFieldFunctionObjects

Description
    Computes Lamb vector, i.e. the cross product of vorticity and velocity.
    Calculation of the divergence of the Lamb vector can be performed by using
    \c div function object on this \c LambVector function object.

    The motivation of the function object is the literature-reported
    quantitative connection between the Lamb vector (divergence) and
    the spatially localised instantaneous fluid motions, e.g. high- and
    low-momentum fluid parcels, which possess considerable level of capacity
    to affect the rate of change of momentum, and to generate forces such as
    drag.

    Operands:
    \table
      Operand          | Type           | Location
      input            | volVectorField | $FOAM_CASE/\<time\>/\<inpField\>
      output file      | -              | -
      output field     | volVectorField | $FOAM_CASE/\<time\>/\<outField\>
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    LambVector1
    {
        // Mandatory entries (unmodifiable)
        type            LambVector;
        libs            (fieldFunctionObjects);

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Req'd | Dflt
      type         | Type name: LambVector              | word |  yes  | -
      libs         | Library name: fieldFunctionObjects | word |  yes  | -
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link fieldExpression.H \endlink

    Minimal example by using the \c postProcess utility:
    \verbatim
        postProcess -func LambVector
    \endverbatim

Note
    To execute \c LambVector function object on an input \<field\>, a numerical
    scheme should be defined for \c div(LambVector) in
    \c system/fvSchemes.divSchemes.

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - Foam::functionObjects::fieldExpression
    - ExtendedCodeGuide::functionObjects::field::LambVector

SourceFiles
    LambVector.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_LambVector_H
#define functionObjects_LambVector_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class LambVector Declaration
\*---------------------------------------------------------------------------*/

class LambVector
:
    public fieldExpression
{
    // Private Member Functions

        //- Calculate the LambVector field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("LambVector");


    // Constructors

        //- Construct from Time and dictionary
        LambVector
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        LambVector(const LambVector&) = delete;

        //- No copy assignment
        void operator=(const LambVector&) = delete;


    //- Destructor
    virtual ~LambVector() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
