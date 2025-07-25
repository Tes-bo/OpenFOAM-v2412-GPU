/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
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
    Foam::LBFGS

Description
    The quasi-Newton Limited-memory BFGS formula. Keeps nPrevSteps_ of the
    y and s vectors and approximates the inverse areas through them.
    Values of 3 < nPrevSteps_ < 20 are suggested.

SourceFiles
    LBFGS.C

\*---------------------------------------------------------------------------*/

#ifndef LBFGS_H
#define LBFGS_H

#include "quasiNewton.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class LBFGS Declaration
\*---------------------------------------------------------------------------*/

class LBFGS
:
    public quasiNewton
{
private:

    // Private Member Functions

        //- No copy construct
        LBFGS(const LBFGS&) = delete;

        //- No copy assignment
        void operator=(const LBFGS&) = delete;


protected:

    // Protected data

        //- Number of old corrections and grad differences kept
        label nPrevSteps_;

        //- The previous differences of derivatives. Holds nPrevSteps_ fields
        PtrList<scalarField> y_;

        //- The previous correction. Holds nPrevSteps_ fields
        PtrList<scalarField> s_;

        //- Use damping for s to ensure positive-definitiveness
        bool useSDamping_;

        //- Use damping for s to ensure positive-definitiveness
        bool useYDamping_;


    // Protected Member Functions

        //- Allocate matrices in the first optimisation cycle
        void allocateVectors();

        //- Move pointers in PtrList to the left and replace last element with
        //- given field
        void pivotFields(PtrList<scalarField>& list, const scalarField& f);

        //- Update y and s vectors
        void updateVectors
        (
            const scalarField& derivatives,
            const scalarField& derivativesOld
        );

        //- Use the damped version of s to ensure positive-definitiveness
        //  Usefull in conjunction with SQP
        void applyDamping
        (
            scalarField& y,
            scalarField& s
        );

        //- Update the Hessian matrix by updating the base vectors
        virtual void updateHessian();

        //- Update design variables
        virtual void update();

        //- Compute the inverse Hessian - vector product.
        //  Input should have the size of all design variables or the active
        //  ones, output is the size of the active design variables.
        //  Can optionally received a seed diagonal matrix
        virtual tmp<scalarField> invHessianVectorProduct
        (
            const scalarField& vector
        );

        //- Same as previous one, but with an explicit counter provided
        tmp<scalarField> invHessianVectorProduct
        (
            const scalarField& vector,
            const label counter,
            const refPtr<scalarField> diag = nullptr
        );

        //- Compute the Hessian - vector product
        //  Requires the solution of a system of equations twice the size of
        //  the bases. This should be OK since the latter is small.
        //  Input should have the size of all design variables or the active
        //  ones, output is the size of the active design variables
        virtual tmp<scalarField> HessianVectorProduct
        (
            const scalarField& vector
        );

        //- Same as previous one, but with an explicit counter provided
        tmp<scalarField> HessianVectorProduct
        (
            const scalarField& vector,
            const label counter
        );

        //- Return the diagonal of the Hessian.
        //  Requires the solution of a system of equations twice the size of
        //  the bases. This should be OK since the latter is small.
        //  Useful for preconditioning
        tmp<scalarField> HessianDiag();


        // Similar functions, but using the SR1 formula instead of BFGS.
        // Should become a separate class at some point

            //- Compute the Hessian - vector product
            //  Requires the solution of a system of equations twice the size
            //  of the bases. This should be OK since the latter is small.
            //  Input should have the size of all design variables or the
            //  active ones, output is the size of the active design variables
            virtual tmp<scalarField> SR1HessianVectorProduct
            (
                const scalarField& vector
            );

            //- Same as previous one, but with an explicit counter provided
            tmp<scalarField> SR1HessianVectorProduct
            (
                const scalarField& vector,
                const label counter
            );

            //- Return the diagonal of the Hessian.
            //  Requires the solution of a system of equations twice the size
            //  of the bases. This should be OK since the latter is small.
            //  Useful for preconditioning
            tmp<scalarField> SR1HessianDiag();


public:

    //- Runtime type information
    TypeName("LBFGS");


    // Constructors

        //- Construct from components
        LBFGS
        (
            const fvMesh& mesh,
            const dictionary& dict,
            autoPtr<designVariables>& designVars,
            const label nConstraints,
            const word& type
        );


    //- Destructor
    virtual ~LBFGS() = default;


    // Member Functions

       //- Write useful quantities to files
       virtual bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
