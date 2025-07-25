/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016 OpenCFD Ltd.
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
    Foam::displacementLayeredMotionMotionSolver

Group
    grpMeshMotionSolvers

Description
    Mesh motion solver for an (multi-block) extruded fvMesh. Gets given the
    structure of the mesh blocks and boundary conditions on these blocks.

    The displacementLayeredMotionCoeffs subdict of dynamicMeshDict specifies
    per region (=cellZone) the boundary conditions on two opposing patches
    (=faceZones). It then interpolates the boundary values using topological
    walking so requires the cellZone to be a layered mesh.

    The boundary conditions on faceZones are currently:

    follow: the faceZone follows the overall mesh displacement.
            Use this for faceZones on boundary faces (so it uses the
            proper boundary conditions on the pointDisplacement).

    uniformFollow: like 'follow' but takes the average value of
            a specified 'patch' (which is not necessarily colocated)

    fixedValue: fixed value.

    timeVaryingUniformFixedValue: table-driven fixed value.

    slip: the second of a pair of faceZones follows the tangential movement
          specified by the first faceZone. (= removes the normal component).

SourceFiles
    displacementLayeredMotionMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef displacementLayeredMotionMotionSolver_H
#define displacementLayeredMotionMotionSolver_H

#include "displacementMotionSolver.H"
#include "cylindricalCS.H"
#include "PtrMap.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class bitSet;


/*---------------------------------------------------------------------------*\
             Class displacementLayeredMotionMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class displacementLayeredMotionMotionSolver
:
    public displacementMotionSolver
{
    // Private Data

       //- Cylindrical coordinate systems (per cellZone)
       PtrMap<coordSystem::cylindrical> cylSystems_;


    // Private Member Functions

        //- Get existing or create new cylindrical system
        const coordSystem::cylindrical& getCylindrical
        (
            const label cellZoneI,
            const dictionary& zoneDict
        );

        void calcZoneMask
        (
            const label cellZoneI,
            bitSet& isZonePoint,
            bitSet& isZoneEdge
        ) const;

        void walkStructured
        (
            const label cellZoneI,
            const bitSet& isZonePoint,
            const bitSet& isZoneEdge,
            const labelList& seedPoints,
            const vectorField& seedData,
            scalarField& distance,
            vectorField& data,
            labelField& patchPoints
        ) const;

        tmp<vectorField> faceZoneEvaluate
        (
            const faceZone& fz,
            const labelList& meshPoints,
            const dictionary& dict,
            const PtrList<pointVectorField>& patchDisp,
            const label patchi
        ) const;

        void cellZoneSolve
        (
            const label cellZoneI,
            const dictionary& zoneDict
        );


        //- No copy construct
        displacementLayeredMotionMotionSolver
        (
            const displacementLayeredMotionMotionSolver&
        ) = delete;

        //- No copy assignment
        void operator=(const displacementLayeredMotionMotionSolver&) = delete;


public:

    //- Runtime type information
    TypeName("displacementLayeredMotion");


    // Constructors

        //- Construct from polyMesh and IOdictionary
        displacementLayeredMotionMotionSolver
        (
            const polyMesh& mesh,
            const IOdictionary& dict
        );

        //- Construct from polyMesh, IOdictionary, displacement and
        //- zero-time fields
        displacementLayeredMotionMotionSolver
        (
            const polyMesh& mesh,
            const IOdictionary& dict,
            const pointVectorField& pointDisplacement,
            const pointIOField& points0
        );


    //- Destructor
    ~displacementLayeredMotionMotionSolver() = default;


    // Member Functions

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve();

        //- Update topology
        virtual void updateMesh(const mapPolyMesh& mpm);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
