/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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
    Foam::parProfilingSolver

Description
    Wrapper to switch on parProfiling around a linear solver.

    Used in combination with parProfiling functionObject.

Usage
    Example of linear solver specification in fvSolution:
    \verbatim
    solvers
    {
        p
        {
            solver          parProfiling;
            // Actual solver to use
            baseSolver      PCG;
            preconditioner  DIC;
            ..
        }
    }
    \endverbatim

SourceFiles
    parProfiling.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_parProfilingSolver_H
#define Foam_parProfilingSolver_H

#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class parProfilingSolver Declaration
\*---------------------------------------------------------------------------*/

class parProfilingSolver
:
    public lduMatrix::solver
{
    // Private Data

        //- Underlying solver
        autoPtr<lduMatrix::solver> solvePtr_;


    // Private Member Functions

        //- No copy construct
        parProfilingSolver(const parProfilingSolver&) = delete;

        //- No copy assignment
        void operator=(const parProfilingSolver&) = delete;


public:

    //- Runtime type information
    TypeName("parProfilingSolver");


    // Constructors

        //- Construct from matrix components and solver controls
        parProfilingSolver
        (
            const word& fieldName,
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& interfaceBouCoeffs,
            const FieldField<Field, scalar>& interfaceIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const dictionary& solverControls
        );


    //- Destructor
    virtual ~parProfilingSolver() = default;


    // Member Functions

        //- Solve the matrix with forwarding to the base solver
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
