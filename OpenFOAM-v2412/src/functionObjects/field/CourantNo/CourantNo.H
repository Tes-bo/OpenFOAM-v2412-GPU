/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
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
    Foam::functionObjects::CourantNo

Group
    grpFieldFunctionObjects

Description
    Computes the Courant number field for time-variant simulations.

    Operands:
    \table
      Operand           | Type           | Location
      input             | -              | -
      output file       | -              | -
      output field      | volScalarField | $FOAM_CASE/\<time\>/\<outField\>
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    CourantNo1
    {
        // Mandatory entries (unmodifiable)
        type        CourantNo;
        libs        (fieldFunctionObjects);

        // Optional entries (runtime modifiable)
        rho         rho;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                        | Type | Req'd | Dflt
      type       | Type name: CourantNo               | word |  yes  | -
      libs       | Library name: fieldFunctionObjects | word |  yes  | -
      rho        | Name of density field              | word |  no   | rho
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link fieldExpression.H \endlink

    Minimal example by using the \c postProcess utility:
    \verbatim
        postProcess -func CourantNo
    \endverbatim

See also
    - Foam::functionObject
    - Foam::functionObjects::fieldExpression
    - Foam::functionObjects::fvMeshFunctionObject
    - ExtendedCodeGuide::functionObjects::field::CourantNo

SourceFiles
    CourantNo.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_CourantNo_H
#define functionObjects_CourantNo_H

#include "fieldExpression.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class CourantNo Declaration
\*---------------------------------------------------------------------------*/

class CourantNo
:
    public fieldExpression
{
    // Private Data

        //- Name of density field
        word rhoName_;


    // Private Member Functions

        //- Divide the Courant number by rho if required
        tmp<volScalarField::Internal> byRho
        (
            const tmp<volScalarField::Internal>& Co
        ) const;

        //- Calculate the Courant number field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("CourantNo");


    // Constructors

        //- Construct from Time and dictionary
        CourantNo
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        CourantNo(const CourantNo&) = delete;

        //- No copy assignment
        void operator=(const CourantNo&) = delete;


    //- Destructor
    virtual ~CourantNo() = default;


    // Member Functions

        //- Read the CourantNo data
        virtual bool read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
