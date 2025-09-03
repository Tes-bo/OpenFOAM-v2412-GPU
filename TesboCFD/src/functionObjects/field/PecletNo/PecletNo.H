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
    Foam::functionObjects::PecletNo

Group
    grpFieldFunctionObjects

Description
    Computes the Peclet number as a \c surfaceScalarField.

    Operands:
    \table
      Operand        | Type               | Location
      input          | surfaceScalarField | $FOAM_CASE/\<time\>/\<inpField\>
      output file    | -                  | -
      output field   | surfaceScalarField | $FOAM_CASE/\<time\>/\<outField\>
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    PecletNo1
    {
        // Mandatory entries (unmodifiable)
        type            PecletNo;
        libs            (fieldFunctionObjects);

        // Optional entries (runtime modifiable)
        rho             rho;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Req'd | Dflt
      type         | Type name: PecletNo                | word |  yes  | -
      libs         | Library name: fieldFunctionObjects | word |  yes  | -
      rho          | Name of density field              | word |  no   | rho
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link fieldExpression.H \endlink

    Minimal example by using the \c postProcess utility:
    \verbatim
        postProcess -func PecletNo
    \endverbatim

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - Foam::functionObjects::fieldExpression
    - ExtendedCodeGuide::functionObjects::field::PecletNo

SourceFiles
    PecletNo.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_PecletNo_H
#define functionObjects_PecletNo_H

#include "fieldExpression.H"
#include "surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class PecletNo Declaration
\*---------------------------------------------------------------------------*/

class PecletNo
:
    public fieldExpression
{
    // Private Data

        //- Name of density field
        word rhoName_;


    // Private Member Functions

        //- Optionally scale the flux for compressible cases
        tmp<surfaceScalarField> rhoScale(const surfaceScalarField& phi) const;

        //- Calculate the Peclet number field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("PecletNo");


    // Constructors

        //- Construct from Time and dictionary
        PecletNo
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        PecletNo(const PecletNo&) = delete;

        //- No copy assignment
        void operator=(const PecletNo&) = delete;


    //- Destructor
    virtual ~PecletNo() = default;


    // Member Functions

        //- Read the PecletNo data
        virtual bool read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
