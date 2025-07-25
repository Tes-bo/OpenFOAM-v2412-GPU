/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2023 PCOpt/NTUA
    Copyright (C) 2021-2023 FOSS GP
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
    Foam::quasiNewton

Description
    Base class for quasi-Newton methods

SourceFiles
    quasiNewton.C

\*---------------------------------------------------------------------------*/

#ifndef quasiNewton_H
#define quasiNewton_H

#include "updateMethod.H"
#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class quasiNewton Declaration
\*---------------------------------------------------------------------------*/

class quasiNewton
:
    public updateMethod
{
protected:

    // Protected data

        //- Step for the Newton method
        scalar etaHessian_;

        //- Number of first steepest descent steps
        label nSteepestDescent_;

        //- Scale the initial unitary Hessian approximation
        bool scaleFirstHessian_;

        //- The Hessian or its inverse, depending on the deriving class.
        //  Has the size of the active design variables.
        //  autoPtr to avoid allocation by limited memory variants
        autoPtr<SquareMatrix<scalar>> Hessian_;

        //- The previous derivatives
        scalarField derivativesOld_;

        //- The previous correction
        scalarField correctionOld_;


    // Protected Member Functions

        //- Update approximation of the inverse Hessian
        virtual void updateHessian() = 0;

        //- Update design variables
        virtual void update() = 0;

        //- Allocate the Hessian matrix
        void allocateHessian();


private:

    // Private Member Functions

        //- No copy construct
        quasiNewton(const quasiNewton&) = delete;

        //- No copy assignment
        void operator=(const quasiNewton&) = delete;


public:

    //- Runtime type information
    TypeName("quasiNewton");


    // Constructors

        //- Construct from components
        quasiNewton
        (
            const fvMesh& mesh,
            const dictionary& dict,
            autoPtr<designVariables>& designVars,
            const label nConstraints,
            const word& type
        );


    //- Destructor
    virtual ~quasiNewton() = default;


    // Member Functions

       //- Compute design variables correction
       void computeCorrection();

       //- Update old correction. Useful for quasi-Newton methods coupled with
       //- line search
       virtual void updateOldCorrection(const scalarField& oldCorrection);

       //- Write useful quantities to files
       virtual bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
