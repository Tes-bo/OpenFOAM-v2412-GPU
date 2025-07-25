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
    Foam::BFGS

Description
    The quasi-Newton BFGS formula.
    quasiNewton::Hessian corresponds to the inverse Hessian matrix in BFGS

SourceFiles
    BFGS.C

\*---------------------------------------------------------------------------*/

#ifndef BFGS_H
#define BFGS_H

#include "quasiNewton.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                             Class BFGS Declaration
\*---------------------------------------------------------------------------*/

class BFGS
:
    public quasiNewton
{
private:

    // Private Member Functions

        //- No copy construct
        BFGS(const BFGS&) = delete;

        //- No copy assignment
        void operator=(const BFGS&) = delete;


protected:

    // Protected data

        //- Curvature threshold
        scalar curvatureThreshold_;


    // Protected Member Functions

        //- Update approximation of the inverse Hessian
        virtual void updateHessian();

        //- Update design variables
        virtual void update();


public:

    //- Runtime type information
    TypeName("BFGS");


    // Constructors

        //- Construct from components
        BFGS
        (
            const fvMesh& mesh,
            const dictionary& dict,
            autoPtr<designVariables>& designVars,
            const label nConstraints,
            const word& type
        );


    //- Destructor
    virtual ~BFGS() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
