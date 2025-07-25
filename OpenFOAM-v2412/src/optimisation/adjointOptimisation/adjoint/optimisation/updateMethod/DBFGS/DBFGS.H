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
    Foam::DBFGS

Description
    The quasi-Newton BFGS formula with the dampening proposed by Powell

SourceFiles
    DBFGS.C

\*---------------------------------------------------------------------------*/

#ifndef DBFGS_H
#define DBFGS_H

#include "quasiNewton.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class DBFGS Declaration
\*---------------------------------------------------------------------------*/

class DBFGS
:
    public quasiNewton
{
protected:

    // Protected data

        //- Curvature threshold
        scalar curvatureThreshold_;

        //- Threshold for damping
        scalar gamma_;


    // Protected Member Functions

        //- Update approximation of the inverse Hessian
        void updateHessian();

        //- Update design variables
        void update();


private:

    // Private Member Functions

        //- No copy construct
        DBFGS(const DBFGS&) = delete;

        //- No copy assignment
        void operator=(const DBFGS&) = delete;


public:

    //- Runtime type information
    TypeName("DBFGS");


    // Constructors

        //- Construct from components
        DBFGS
        (
            const fvMesh& mesh,
            const dictionary& dict,
            autoPtr<designVariables>& designVars,
            const label nConstraints,
            const word& type
        );


    //- Destructor
    virtual ~DBFGS() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
