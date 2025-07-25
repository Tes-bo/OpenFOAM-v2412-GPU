/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
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
    Foam::functionObjects::flowType

Group
    grpFieldFunctionObjects

Description
    Computes the flow type indicator of an input velocity field.

    The flow type indicator is obtained according to the following equation:
    \f[
        \lambda = \frac{|D| - |\omega|}{|D| + |\omega|}
    \f]

    where
    \vartable
      \lambda | Flow type indicator
      D       | Symmetric part of the gradient tensor of velocity
      \omega  | Skew-symmetric part of the gradient tensor of velocity
    \endvartable

    The flow type indicator values mean:
    \verbatim
      -1 = rotational flow
       0 = simple shear flow
       1 = planar extensional flow
    \endverbatim

    Operands:
    \table
      Operand          | Type           | Location
      input            | volVectorField | $FOAM_CASE/\<time\>/\<inpField\>
      output file      | -              | -
      output field     | volScalarField | $FOAM_CASE/\<time\>/\<outField\>
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    flowType1
    {
        // Mandatory entries (unmodifiable)
        type            flowType;
        libs            (fieldFunctionObjects);

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Req'd | Dflt
      type         | Type name: flowType                | word |  yes  | -
      libs         | Library name: fieldFunctionObjects | word |  yes  | -
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link fieldExpression.H \endlink

    Minimal example by using the \c postProcess utility:
    \verbatim
        postProcess -func flowType
    \endverbatim

See also
    - Foam::functionObject
    - Foam::functionObjects::fieldExpression
    - Foam::functionObjects::fvMeshFunctionObject
    - ExtendedCodeGuide::functionObjects::field::flowType

SourceFiles
    flowType.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_flowType_H
#define functionObjects_flowType_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class flowType Declaration
\*---------------------------------------------------------------------------*/

class flowType
:
    public fieldExpression
{
    // Private Member Functions

        //- Calculate the flowType field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("flowType");


    // Constructors

        //- Construct from Time and dictionary
        flowType
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        flowType(const flowType&) = delete;

        //- No copy assignment
        void operator=(const flowType&) = delete;


    //- Destructor
    virtual ~flowType() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
