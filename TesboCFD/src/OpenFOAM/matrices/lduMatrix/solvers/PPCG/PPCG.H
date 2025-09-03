/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 Mattijs Janssens
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
    Foam::PPCG

Description
    Preconditioned pipelined conjugate gradient solver for symmetric
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

SourceFiles
    PPCG.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_PPCG_H
#define Foam_PPCG_H

#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class PPCG Declaration
\*---------------------------------------------------------------------------*/

class PPCG
:
    public lduMatrix::solver
{
    // Private Member Data

        //- Cached preconditioner
        mutable autoPtr<lduMatrix::preconditioner> preconPtr_;


    // Private Member Functions

        //- Non-blocking version of sum(a*b), sum(a*c), sum(mag(sumMag))
        static void gSumMagProd
        (
            FixedList<solveScalar, 3>& globalSum,
            const solveScalarField& a,
            const solveScalarField& b,
            const solveScalarField& c,
            const solveScalarField& sumMag,
            UPstream::Request& request,
            const label comm
        );

        //- No copy construct
        PPCG(const PPCG&) = delete;

        //- No copy assignment
        void operator=(const PPCG&) = delete;


protected:

        //- CG solver. Operates either in conjugate-gradient mode or
        //- conjugate residual
        solverPerformance scalarSolveCG
        (
            solveScalarField& psi,
            const solveScalarField& source,
            const direction cmpt,
            const bool cgMode
        ) const;


public:

    //- Runtime type information
    TypeName("PPCG");


    // Constructors

        //- Construct from matrix components and solver controls
        PPCG
        (
            const word& fieldName,
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& interfaceBouCoeffs,
            const FieldField<Field, scalar>& interfaceIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const dictionary& solverControls
        );


    //- Destructor
    virtual ~PPCG() = default;


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
