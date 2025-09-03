/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::FDICPreconditioner

Group
    grpLduMatrixPreconditioners

Description
    Faster version of the DICPreconditioner diagonal-based incomplete
    Cholesky preconditioner for symmetric matrices
    (symmetric equivalent of DILU) in which the reciprocal of the
    preconditioned diagonal and the upper coefficients divided by the diagonal
    are calculated and stored.

SourceFiles
    FDICPreconditioner.C

\*---------------------------------------------------------------------------*/

#ifndef FDICPreconditioner_H
#define FDICPreconditioner_H

#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class FDICPreconditioner Declaration
\*---------------------------------------------------------------------------*/

class FDICPreconditioner
:
    public lduMatrix::preconditioner
{
    // Private data

        //- The reciprocal preconditioned diagonal
        solveScalarField rD_;
        solveScalarField rDuUpper_;
        solveScalarField rDlUpper_;


    // Private Member Functions

        //- No copy construct
        FDICPreconditioner(const FDICPreconditioner&) = delete;

        //- No copy assignment
        void operator=(const FDICPreconditioner&) = delete;


public:

    //- Runtime type information
    TypeName("FDIC");


    // Constructors

        //- Construct from matrix components and preconditioner solver controls
        FDICPreconditioner
        (
            const lduMatrix::solver&,
            const dictionary& solverControlsUnused
        );


    //- Destructor
    virtual ~FDICPreconditioner() = default;


    // Member Functions

        //- Return wA the preconditioned form of residual rA
        virtual void precondition
        (
            solveScalarField& wA,
            const solveScalarField& rA,
            const direction cmpt=0
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
