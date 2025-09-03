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
    Foam::distributedDICPreconditioner

Group
    grpLduMatrixPreconditioners

Description
    Version of DICpreconditioner that uses preconditioning across processor
    (and coupled) boundaries. Based on 'Parallel Preconditioners' chapter from
    'Iterative Methods for Sparse Linear Systems' by Yousef Saad.

    Leaves out the handling of boundary nodes after internal nodes since
    probably not beneficial (so no overlap of comms and internal calculation)

    By default preconditions across coupled boundaries (currently only
    tested for cyclicAMI). This can be disabled with the 'coupled' setting

        solver          PCG;
        preconditioner
        {
            preconditioner  distributedDIC;
            coupled         false;
        }

    The cyclicAMI boundary behaviour will only work if
    - running non-parallel or
    - different sides of cyclicAMI run on different processors i.e.
      there is no processor which has cells on both owner and neighbour
      of the patch pair.

See Also
    Foam::DICPreconditioner
    Foam::distributedDILUPreconditioner

SourceFiles
    distributedDICPreconditioner.C

\*---------------------------------------------------------------------------*/

#ifndef distributedDICPreconditioner_H
#define distributedDICPreconditioner_H

#include "distributedDILUPreconditioner.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class distributedDICPreconditioner Declaration
\*---------------------------------------------------------------------------*/

class distributedDICPreconditioner
:
    public distributedDILUPreconditioner
{
protected:

    // Protected member functions

        //- Update diagonal for colour
        virtual void forwardInternalDiag
        (
            solveScalarField& rD,
            const label colouri
        ) const;

        //- Update preconditioned variable walking forward on internal faces
        virtual void forwardInternal
        (
            solveScalarField& wA,
            const label colouri
        ) const;


public:

    //- Runtime type information
    TypeName("distributedDIC");


    // Constructors

        //- Construct from matrix components and preconditioner solver controls
        distributedDICPreconditioner
        (
            const lduMatrix::solver&,
            const dictionary& solverControlsUnused
        );


    //- Destructor
    virtual ~distributedDICPreconditioner() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
