/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2021 PCOpt/NTUA
    Copyright (C) 2013-2021 FOSS GP
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
    Foam::SQP

Description
    The quasi-Newton SQP formula for constrained optimisation

SourceFiles
    SQP.C

\*---------------------------------------------------------------------------*/

#ifndef SQP_H
#define SQP_H

#include "quasiNewton.H"
#include "SQPBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                             Class SQP Declaration
\*---------------------------------------------------------------------------*/

class SQP
:
    public quasiNewton,
    public SQPBase
{
protected:

    // Protected data

        //- Curvature threshold
        scalar dumpingThreshold_;


private:

    // Private Member Functions

        //- No copy construct
        SQP(const SQP&) = delete;

        //- No copy assignment
        void operator=(const SQP&) = delete;

        //- Update the Hessian matrix
        virtual void updateHessian();

        //- Compute new lamdas and update correction
        virtual void update();

        //- Store old fields needed for the next iter
        void storeOldFields();

        //- Get the part the merit function that depends on the constraints
        virtual scalar meritFunctionConstraintPart() const;


public:

    //- Runtime type information
    TypeName("SQP");


    // Constructors

        //- Construct from components
        SQP
        (
            const fvMesh& mesh,
            const dictionary& dict,
            autoPtr<designVariables>& designVars,
            const label nConstraints,
            const word& type
        );


    //- Destructor
    virtual ~SQP() = default;


    // Member Functions

        //- Compute design variables correction
        void computeCorrection();

        //- Compute merit function. Could be different than the objective
        //- in the presence of constraints
        virtual scalar computeMeritFunction();

        //- Derivative of the merit function. Could be different than the
        //- objective derivative in the presence of constraints
        virtual scalar meritFunctionDirectionalDerivative();

        //- Write useful quantities to files
        virtual bool writeData(Ostream& os) const;

        //- Write merit function information
        virtual bool writeAuxiliaryData();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
