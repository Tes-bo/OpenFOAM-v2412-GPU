/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2020 PCOpt/NTUA
    Copyright (C) 2013-2020 FOSS GP
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
    Foam::volumetricBSplinesMotionSolver

Description
    A mesh motion solver based on volumetric B-Splines

SourceFiles
    volumetricBSplinesMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef volumetricBSplinesMotionSolver_H
#define volumetricBSplinesMotionSolver_H

#include "motionSolver.H"
#include "volBSplinesBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                Class volumetricBSplinesMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class volumetricBSplinesMotionSolver
:
    public motionSolver
{
protected:

    // Protected data

        //- Reference to the underlaying volumetric B-Splines base
        volBSplinesBase& volBSplinesBase_;

        //- Movement of the control points
        vectorField controlPointsMovement_;

        //- Number of control boxes
        wordList controlBoxes_;

        //- List with volumetric B-splines boxes. No overlapping is supported
        mutable List<autoPtr<NURBS3DVolume>> volume_;

        //- Number of control boxes
        label nBoxes_;


private:

    // Private Member Functions

        //- No copy construct
        volumetricBSplinesMotionSolver
        (
            const volumetricBSplinesMotionSolver&
        ) = delete;

        //- No copy assignment
        void operator=(const volumetricBSplinesMotionSolver&) = delete;


public:

    //- Runtime type information
    TypeName("volumetricBSplinesMotionSolver");


    // Constructors

        //- Construct from mesh and dictionary
        volumetricBSplinesMotionSolver
        (
            const polyMesh& mesh,
            const IOdictionary& dict
        );


    //- Destructor
    virtual ~volumetricBSplinesMotionSolver() = default;


    // Member Functions

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion/ Does nothing
        virtual void solve();

        //- Update local data for geometry changes
        virtual void movePoints(const pointField&);

        //- Update the mesh corresponding to given map
        virtual void updateMesh(const mapPolyMesh&);

        //- Set control points movement
        void setControlPointsMovement(const vectorField& controlPointsMovement);

        //- Bound control points movement
        void boundControlPointMovement(vectorField& controlPointsMovement);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
