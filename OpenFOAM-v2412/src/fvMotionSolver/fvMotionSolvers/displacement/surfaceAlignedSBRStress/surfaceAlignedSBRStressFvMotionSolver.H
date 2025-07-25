/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2015-2016 OpenCFD Ltd.
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
    Foam::surfaceAlignedSBRStressFvMotionSolver

Group
    grpMeshMotionSolvers

Description
    Mesh motion solver for an fvMesh.  Based on solving the cell-centre
    solid-body rotation stress equations for the motion displacement.
    The model calculates the necessary rotation to align the stl surface
    with the closest mesh face normals and it calculates the respective source
    term for the SBRStress equation.

    Input parameters:

    \verbatim
    surfaceAlignedSBRStressCoeffs
    {
        diffusivity         uniform;
        surfaces            ("terrain.stl");
        maxAng              85;
        minAng              5;
        nNonOrthogonalCorr  4;
        accFactor           1.3;
        smoothFactor        0.9;
        minSigmaDiff        1e-4;
    }
    \endverbatim

    surfaces: name of the stl surfaces to which apply alignment.

    maxAng: Maximum angle (between surface and most aligned mesh face normal)
            to which rotation is applied. (default 80 degress)

    minAng: Minimum angle (between surface and most aligned mesh face normal)
            to which rotation is applied (default 20 degress)

    nNonOrthogonalCorr: Non-orthogonal correction of the SBRStress equation

    accFactor: Proportionality constant applied to the source of the stress
               equation.(default 1)

    smoothFactor: Displacement smooth factor (1 very smooth , 0 no smoothing)
                (default 0.9)

    minSigmaDiff: Minimum magnitude difference between old and new magnitudes
                 of the applied stress.

                 If the difference between old and new is smaller than
                 minSigmaDiff no further addition to the stress source term
                 will be made (default: 1e-3)

SourceFiles
    surfaceAlignedSBRStressFvMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceAlignedSBRStressFvMotionSolver_H
#define surfaceAlignedSBRStressFvMotionSolver_H

#include "triSurfaceMesh.H"
#include "vectorList.H"
#include "displacementSBRStressFvMotionSolver.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class surfaceAlignedSBRStressFvMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class surfaceAlignedSBRStressFvMotionSolver
:
    public displacementSBRStressFvMotionSolver
{
    // Private data

        //- Names of the surfaces
        wordList surfaceNames_;

        //- TriSurfaceMesh
        PtrList<triSurfaceMesh>  surfaceMesh_;

        //- Rotation vector field
        volVectorField cellRot_;

        //- Maximum angle (between surface and most aligned mesh face normal)
        //- to which rotation is applied
        scalar maxAng_;

        //- Minimum angle (between surface and most aligned mesh face normal)
        //- to which rotation is applied
        scalar minAng_;

        //- Proportional constant applied to the source to accelerate
        //- convergence
        scalar accFactor_;

        //- Displacement smooth factor
        scalar smoothFactor_;

        //- Non-orthogonal correction of the SBRStress equation
        label nNonOrthogonalCorr_;

        //- Point displacement field
        pointVectorField& pointDisplacement_;

        //- Stress tensor
        volSymmTensorField sigmaD_;

        //- Minimum magnitude difference between old and new magnitudes
        //  stress applied
        scalar minSigmaDiff_;


    // Private Member Functions

        //- No copy construct
        surfaceAlignedSBRStressFvMotionSolver
        (
            const surfaceAlignedSBRStressFvMotionSolver&
        ) = delete;

        //- No copy assignment
        void operator=(const surfaceAlignedSBRStressFvMotionSolver&) = delete;

        //- Calculate cellRot
        void calculateCellRot();


public:

    //- Runtime type information
    TypeName("surfaceAlignedSBRStress");


    // Constructors

        //- Construct from polyMesh and IOdictionary
        surfaceAlignedSBRStressFvMotionSolver
        (
            const polyMesh&,
            const IOdictionary&
        );


    //- Destructor
    ~surfaceAlignedSBRStressFvMotionSolver();


    // Member Functions

        //- Solve for motion
        virtual void solve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
