/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 Zeljko Tukovic, FSB Zagreb.
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::interfaceTrackingFvMesh

Description
    The interfaceTrackingFvMesh

SourceFiles
    interfaceTrackingFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceTrackingFvMesh_H
#define interfaceTrackingFvMesh_H

#include "dynamicMotionSolverFvMesh.H"
#include "regIOobject.H"
#include "faCFD.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "surfactantProperties.H"
#include "singlePhaseTransportModel.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class interfaceTrackingFvMesh Declaration
\*---------------------------------------------------------------------------*/

class interfaceTrackingFvMesh
:
    public dynamicMotionSolverFvMesh
{
    // Private Data

        //- Finite area mesh
        autoPtr<faMesh> aMeshPtr_;

        //- Free surface patch index
        label fsPatchIndex_;

        //- Free surface faPatch-es which do not move
        wordList fixedFreeSurfacePatches_;

        //- Free surface faPatch-es where wave shuld not reflect
        wordList nonReflectingFreeSurfacePatches_;

        //- Free surface patches for which point normals must be corrected
        wordList pointNormalsCorrectionPatches_;

        //- Is it free-surface points displacement direction
        //  parallel with free-surface point normals
        Switch normalMotionDir_;

        //- Free-surface points displacement direction
        //  if not normal motion direction
        vector motionDir_;

        //- Interface smoothing at the beginning of time step
        Switch smoothing_;

        //- Pure free-surface
        Switch pureFreeSurface_;

        //- Rigid free-surface
        Switch rigidFreeSurface_;

        //- Correct contact line point normals
        Switch correctContactLineNormals_;

        //- Surface tension coefficient of pure free-surface
        dimensionedScalar sigma0_;

        //- Fluid density
        dimensionedScalar rho_;

        //- Current interface tracking time index
        label timeIndex_;

        //- Free-surface velocity field
        mutable areaVectorField* UsPtr_;

        //- Points which are attached to the free-surface A side faces
        //  and which defines the free-surface shape
        mutable vectorIOField* controlPointsPtr_;

        //- Field which additionally determines
        //  the motion of free-surface points
        mutable labelList* motionPointsMaskPtr_;

        //- Displacement direction of free-surface points
        mutable vectorField* pointsDisplacementDirPtr_;

        //- Displacement direction of free-surface control points
        mutable vectorField* facesDisplacementDirPtr_;

        //- Free-surface net flux
        mutable areaScalarField* fsNetPhiPtr_;

        //- Free-surface flux
        mutable edgeScalarField* phisPtr_;

        //- Free-surface surfactant concetration
        mutable areaScalarField* surfactConcPtr_;

        //- Volume surfactant concetration
        mutable volScalarField* bulkSurfactConcPtr_;

        //- Surface tension field
        mutable areaScalarField* surfaceTensionPtr_;

        //- Surfactant properties
        mutable surfactantProperties* surfactantPtr_;

        //- Contact angle
        mutable areaScalarField* contactAnglePtr_;

    // Private Member Functions

        // Make demand-driven data

            //- Create free surface velocity field
            void makeUs() const;

            //- Create control points
            void makeControlPoints();

            //- Create motion points mask
            void makeMotionPointsMask();

            //- Create points and control point motion direction
            void makeDirections();

            //- Create free surface net flux
            void makeFsNetPhi() const;

            //- Create free-surface flux
            void makePhis();

            //- Create surfactant volume concentration field
            void makeBulkSurfactConc() const;

            //- Create surfactant concentration field
            void makeSurfactConc() const;

            //- Create surface tension field
            void makeSurfaceTension() const;

            //- Create surfactant properties
            void makeSurfactant() const;

            //- Create contact angle
            void makeContactAngle();

        //- No copy construct
        interfaceTrackingFvMesh(const interfaceTrackingFvMesh&) = delete;

        //- No copy assignment
        void operator=(const interfaceTrackingFvMesh&) = delete;

        //- Initialize data
        void initializeData();

        //- Update control points end displacement directions
        void updateDisplacementDirections();

        //- Initialize control points position
        void initializeControlPointsPosition();

        //- Calculate free surface points displacement
        //  for given new control points position
        tmp<vectorField> pointDisplacement();

        //- Calc least sqare plane point and normal
        tmp<vectorField> lsPlanePointAndNormal
        (
            const vectorField& points,
            const vector& origin,
            const vector& axis
        ) const;

        //- Smooth free-surface mesh
        void smoothFreeSurfaceMesh();

        //- Update free-surface flux
        void updateSurfaceFlux();

        //- Update free-surface surfactant concentration
        void updateSurfactantConcentration();

        //- Calculate total pressure force
        vector totalPressureForce() const;

        //- Calculate total viscous force
        vector totalViscousForce() const;

        //- Calculate total surface tension force
        vector totalSurfaceTensionForce() const;

        //- Maximal surface tension based Courant number
        scalar maxCourantNumber();

        //- Update properties
        void updateProperties();

        //- Correct free surface point normals at contact line
        void correctContactLinePointNormals();

        //- Correct free surface point displacement next to fixed contact line
        void correctPointDisplacement
        (
            const scalarField& sweptVolCorr,
            vectorField& displacement
        );

public:

    //- Runtime type information
    TypeName("interfaceTrackingFvMesh");


    // Constructors

        //- Construct from IOobject
        interfaceTrackingFvMesh(const IOobject& io, const bool doInit=true);
/*
        //- Construct from components without boundary.
        //  Boundary is added using addFvPatches() member function
        interfaceTrackingFvMesh
        (
            const IOobject& io,
            pointField&& points,
            faceList&& faces,
            labelList&& allOwner,
            labelList&& allNeighbour,
            const bool syncPar = true
        );
*/

    //- Destructor
    ~interfaceTrackingFvMesh();


    // Member Functions

        //- Initialise all non-demand-driven data
        virtual bool init(const bool doInit);

        fvMesh& mesh()
        {
            return *this;
        }

        const fvMesh& mesh() const
        {
            return *this;
        }

        //- Return reference to finite area mesh
        faMesh& aMesh()
        {
            return aMeshPtr_();
        }

        //- Return reference to finite area mesh
        const faMesh& aMesh() const
        {
            return aMeshPtr_();
        }

        const label& fsPatchIndex() const
        {
            return fsPatchIndex_;
        }

        //- Pure free-surface
        Switch pureFreeSurface() const noexcept
        {
            return pureFreeSurface_;
        }

        //- Rigid free-surface
        Switch rigidFreeSurface() const noexcept
        {
            return rigidFreeSurface_;
        }

        //- Rigid free-surface
        Switch& rigidFreeSurface() noexcept
        {
            return rigidFreeSurface_;
        }

        //- Correct contact line point normals
        Switch correctContactLineNormals() const noexcept
        {
            return correctContactLineNormals_;
        }

        //- Correct contact line point normals
        Switch& correctContactLineNormals() noexcept
        {
            return correctContactLineNormals_;
        }

        //- Surface tension coefficient of pure free-surface
        const dimensionedScalar& sigma() const
        {
            return sigma0_;
        }

        //- Return free-surface velocity field
        areaVectorField& Us();

        //- Return free-surface velocity field
        const areaVectorField& Us() const;

        //- Return free-surface net flux
        const areaScalarField& fsNetPhi() const;

        //- Return free-surface net flux
        areaScalarField& fsNetPhi();

        //- Correct surface velocity boundary conditions
        void correctUsBoundaryConditions();

        //- Update free-surface velocity field
        void updateUs();

        //- Return constant reference to velocity field
        const volVectorField& U() const;

        //- Return constant reference to pressure field
        const volScalarField& p() const;

        //- Return constant reference to velocity field
        const surfaceScalarField& phi() const;

        //- Return free surface normal derivative of velocity
        tmp<vectorField> freeSurfaceSnGradU();

        //- Return free surface normal derivative of normal velocity comp
        tmp<scalarField> freeSurfaceSnGradUn();

        //- Return free surface pressure jump
        tmp<scalarField> freeSurfacePressureJump();

        //- Return control points
        vectorField& controlPoints();

        //- Return reference to motion points mask field
        labelList& motionPointsMask();

        //- Return reference to point displacement direction field
        vectorField& pointsDisplacementDir();

        //- Return reference to control points displacement direction field
        vectorField& facesDisplacementDir();

        //- Motion direction swithc
        bool normalMotionDir() const
        {
            return normalMotionDir_;
        }

        //- Return free-surface fluid flux field
        edgeScalarField& Phis();

        //- Return free-surface surfactant concentration field
        areaScalarField& surfactantConcentration();

        //- Return free-surface surfactant concentration field
        const areaScalarField& surfactantConcentration() const;

        //- Return volume surfactant concentration field
        volScalarField& bulkSurfactantConcentration();

        //- Return volume surfactant concentration field
        const volScalarField& bulkSurfactantConcentration() const;

        //- Return surface tension field
        areaScalarField& surfaceTension();

        //- Return surface tension field
        const areaScalarField& surfaceTension() const;

        //- Return surface tension gradient
        tmp<areaVectorField> surfaceTensionGrad();

        //- Return surfactant properties
        const surfactantProperties& surfactant() const;

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();

        //- Clear control points
        void clearControlPoints()
        {
            deleteDemandDrivenData(controlPointsPtr_);
        }

        //- Write VTK freeSurface mesh
        void writeVTK() const;

        //- Write VTK freeSurface control points
        void writeVTKControlPoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
