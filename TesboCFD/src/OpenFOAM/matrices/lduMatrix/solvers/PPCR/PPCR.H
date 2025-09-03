/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 Mattijs Janssens
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
    Foam::PPCR

Description
    Preconditioned pipelined conjugate residuals solver for symmetric
    lduMatrices using a run-time selectable preconditioner.

    Reference:
    \verbatim
        P. Ghysels, W. Vanroose.
        "Hiding global synchronization latency in the
         preconditioned Conjugate Gradient algorithm"
    \endverbatim
    and implementation details from
    \verbatim
        Paul Eller, William Gropp
        "Scalable Non-blocking Preconditioned Conjugate Gradient Methods"
    \endverbatim

See also
    PPCG

SourceFiles
    PPCR.C

\*---------------------------------------------------------------------------*/

#ifndef PPCR_H
#define PPCR_H

#include "PPCG.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class PPCR Declaration
\*---------------------------------------------------------------------------*/

class PPCR
:
    public PPCG
{
    // Private Member Functions

        //- No copy construct
        PPCR(const PPCR&) = delete;

        //- No copy assignment
        void operator=(const PPCR&) = delete;


public:

    //- Runtime type information
    TypeName("PPCR");


    // Constructors

        //- Construct from matrix components and solver controls
        PPCR
        (
            const word& fieldName,
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& interfaceBouCoeffs,
            const FieldField<Field, scalar>& interfaceIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const dictionary& solverControls
        );


    //- Destructor
    virtual ~PPCR() = default;


    // Member Functions

        //- Solve the matrix with this solver
        virtual solverPerformance solve
        (
            scalarField& psi,
            const scalarField& source,
            const direction cmpt=0
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
