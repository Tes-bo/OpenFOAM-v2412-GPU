/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::PBiCCCG

Description
    Preconditioned bi-conjugate gradient solver for asymmetric lduMatrices
    using a run-time selectable preconditioner.

SourceFiles
    PBiCCCG.C

\*---------------------------------------------------------------------------*/

#ifndef PBiCCCG_H
#define PBiCCCG_H

#include "LduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class PBiCCCG Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class DType, class LUType>
class PBiCCCG
:
    public LduMatrix<Type, DType, LUType>::solver
{
    // Private Member Functions

        //- No copy construct
        PBiCCCG(const PBiCCCG&) = delete;

        //- No copy assignment
        void operator=(const PBiCCCG&) = delete;


public:

    //- Runtime type information
    TypeName("PBiCCCG");


    // Constructors

        //- Construct from matrix components and solver data dictionary
        PBiCCCG
        (
            const word& fieldName,
            const LduMatrix<Type, DType, LUType>& matrix,
            const dictionary& solverDict
        );


    // Destructor

        virtual ~PBiCCCG() = default;


    // Member Functions

        //- Solve the matrix with this solver
        virtual SolverPerformance<Type> solve(Field<Type>& psi) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PBiCCCG.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
