/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
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
    Foam::displacementMethodvolumetricBSplinesMotionSolver

Description
    Wrapper class for the volumetricBSplinesMotionSolver motion solver

SourceFiles
    displacementMethodvolumetricBSplinesMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef displacementMethodvolumetricBSplinesMotionSolver_H
#define displacementMethodvolumetricBSplinesMotionSolver_H

#include "displacementMethod.H"
#include "volumetricBSplinesMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
       Class displacementMethodvolumetricBSplinesMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class displacementMethodvolumetricBSplinesMotionSolver
:
    public displacementMethod
{
private:

    // Private Member Functions

        //- No copy construct
        displacementMethodvolumetricBSplinesMotionSolver
        (
            const displacementMethodvolumetricBSplinesMotionSolver&
        ) = delete;

        //- No copy assignment
        void operator=
        (
            const displacementMethodvolumetricBSplinesMotionSolver&
        ) = delete;


public:

    //- Runtime type information
    TypeName("volumetricBSplinesMotionSolver");


    // Constructors

        //- Construct from components
        displacementMethodvolumetricBSplinesMotionSolver
        (
            fvMesh& mesh,
            const labelList& patchIDs
        );


    //- Destructor
    virtual ~displacementMethodvolumetricBSplinesMotionSolver() = default;


    // Member Functions

       //- Set motion filed related to model based on given motion
       void setMotionField(const pointVectorField& pointMovement);

       //- Set motion filed related to model based on given motion
       void setMotionField(const volVectorField& cellMovement);

       //- Set control field as a vectorField. For methods working with
       //- parameters (RBF etc)
       void setControlField(const vectorField& controlField);

       //- Set control field as a vectorField. For methods working with
       //- parameters (RBF etc)
       void setControlField(const scalarField& controlField);

       //- Bound control field in certain directions etc. For methods working
       //- with parameters (RBF etc)
       void boundControlField(vectorField& controlField);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
