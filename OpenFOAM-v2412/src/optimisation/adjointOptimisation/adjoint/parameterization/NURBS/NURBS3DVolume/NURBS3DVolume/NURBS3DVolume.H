/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    Foam::NURBS3DVolume

Description
    NURBS3DVolume morpher. Includes support functions for gradient computations
    Base class providing support for different coordinate systems

    Reference:
    \verbatim
        For a short introduction to a volumetric B-Splines morpher and its use
        in shape optimisation

            Papoutsis-Kiachagias, E., Magoulas, N., Mueller, J.,
            Othmer, C., & Giannakoglou, K. (2015).
            Noise reduction in car aerodynamics using a surrogate objective
            function and the continuous adjoint method with wall functions.
            Computers & Fluids, 122, 223-232.
            http://doi.org/10.1016/j.compfluid.2015.09.002
    \endverbatim

SourceFiles
    NURBS3DVolume.C
    NURBS3DVolumeI.H

\*---------------------------------------------------------------------------*/

#ifndef NURBS3DVolume_H
#define NURBS3DVolume_H

#include "NURBSbasis.H"
#include "boolVector.H"
#include "localIOdictionary.H"
#include "pointMesh.H"
#include "pointPatchField.H"
#include "pointPatchFieldsFwd.H"
#include "boundaryFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class NURBS3DVolume Declaration
\*---------------------------------------------------------------------------*/

class NURBS3DVolume
:
    public localIOdictionary
{
protected:

    // Protected Data

        typedef List<boolVector> boolVectorList;

        const fvMesh& mesh_;
        dictionary dict_;
        word name_;
        //- NURBS basis functions
        NURBSbasis basisU_;
        NURBSbasis basisV_;
        NURBSbasis basisW_;

        //- Max iterations for Newton-Raphson
        label maxIter_;

        //- Tolerance for Newton-Raphson
        scalar tolerance_;

        //- How many times to bound parametric coordinates until deciding it is
        //- outside the box
        label nMaxBound_;

        //- The volumetric B-Splines control points
        vectorField cps_;

        //- Map of points-in-box to mesh points
        autoPtr<labelList> mapPtr_;

        //- Map of mesh points to points-in-box
        //  Return -1 if point not within the box
        autoPtr<labelList> reverseMapPtr_;

        //- Parametric coordinates of pointsInBox
        autoPtr<pointVectorField> parametricCoordinatesPtr_;

        //- Coordinates in the local system for which CPs are defined
        vectorField localSystemCoordinates_;

        //- Confine movement in certain directions and control points. Refers
        //- to the local system
        label confineUMovement_;

        label confineVMovement_;

        label confineWMovement_;

        label confineBoundaryControlPoints_;

        //- Which movement components to freeze in each plane
        boolVectorList confineUMinCPs_;

        boolVectorList confineUMaxCPs_;

        boolVectorList confineVMinCPs_;

        boolVectorList confineVMaxCPs_;

        boolVectorList confineWMinCPs_;

        boolVectorList confineWMaxCPs_;

        //- Which of the cps are moved in an optimisation
        boolList activeControlPoints_;

        //- Which design variables are changed in an optimisation
        boolList activeDesignVariables_;

        //- Folder to store control points
        string cpsFolder_;

        //- Read parametric coordinates from file if present in the folder
        bool readStoredData_;


    // Protected Member Functions

        //- Find points within control points box
        void findPointsInBox(const vectorField& meshPoints);

        //- Compute parametric coordinates in order to match a given set
        //- of coordinates, based on the cps of the class
        // Uses a Newton-Raphson loop.
        // Argument is the points residing in the box
        void computeParametricCoordinates(const vectorField& points);

        void computeParametricCoordinates(tmp<vectorField> tPoints);

        //- Bound components to certain limits
        bool bound
        (
            vector& vec,
            scalar minValue = 1e-7,
            scalar maxValue = 0.999999
        ) const;

        //- Create lists with active design variables and control points
        void determineActiveDesignVariablesAndPoints();

        //- Confine movement in boundary control points if necessary
        void confineBoundaryControlPoints();

        //- Confine control point movement to maintain user-defined continuity
        void continuityRealatedConfinement();

        //- Confine movement in all control points for user-defined directions
        void confineControlPointsDirections();

        //- Deactivate control points if they affect no mesh point
        void confineInertControlPoints();

        //- Confine all three movements for a prescribed control point
        void confineControlPoint(const label cpI);

        //- Confine specific movements for a prescribed control point
        void confineControlPoint(const label cpI, const boolVector&);

        //- Create folders to store cps and derivatives
        void makeFolders();

        //- Transform a point from its coordinate system to a cartesian system
        virtual vector transformPointToCartesian
        (
            const vector& localCoordinates
        ) const = 0;

        //- Transformation tensor for dxdb, from local coordinate system to
        //- cartesian
        virtual tensor transformationTensorDxDb(label globalPointIndex) = 0;

        //- Update coordinates in the local system based on the cartesian points
        virtual void updateLocalCoordinateSystem
        (
            const vectorField& cartesianPoints
        ) = 0;


public:

    //- Runtime type information
    TypeName("NURBS3DVolume");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            NURBS3DVolume,
            dictionary,
            (
                const dictionary& dict,
                const fvMesh& mesh,
                bool computeParamCoors
            ),
            (dict, mesh, computeParamCoors)
        );

    // Constructors

        //- Construct from dictionary
        NURBS3DVolume
        (
            const dictionary& dict,
            const fvMesh& mesh,
            bool computeParamCoors = true
        );

        //- Construct as copy
        NURBS3DVolume(const NURBS3DVolume&);


    // Selectors

        //- Return a reference to the selected NURBS model
        static autoPtr<NURBS3DVolume> New
        (
            const dictionary& dict,
            const fvMesh& mesh,
            bool computeParamCoors = true
        );


    //- Destructor
    virtual ~NURBS3DVolume() = default;


    // Member Functions

        //- Compute parametric coordinates for a given set of points
        //  (not necessarily the mesh ones)
        tmp<vectorField> computeParametricCoordinates
        (
            const vectorField& points
        ) const;

        // Derivatives wrt parametric coordinates

            //- Volume point derivative wrt u at point u,v,w
            vector volumeDerivativeU
            (
                const scalar u,
                const scalar v,
                const scalar w
            ) const;

            //- Volume point derivative wrt v at point u,v,w
            vector volumeDerivativeV
            (
                const scalar u,
                const scalar v,
                const scalar w
            ) const;

            //- Volume point derivative wrt w at point u,v,w
            vector volumeDerivativeW
            (
                const scalar u,
                const scalar v,
                const scalar w
            ) const;

            //- Jacobian matrix wrt to the volume parametric coordinates
            tensor JacobianUVW(const vector& u) const;


        // Derivatives wrt to control points

            //- Volume point derivative wrt to control point cpI at point u,v,w
            //  Scalar since in the local system!
            scalar volumeDerivativeCP(const vector& u, const label cpI) const;

            //- Control point sensitivities computed using point-based surface
            //- sensitivities
            vectorField computeControlPointSensitivities
            (
                const pointVectorField& pointSens,
                const labelList& sensitivityPatchIDs
            );

            //- Control point sensitivities computed using face-based surface
            //- sensitivities
            vectorField computeControlPointSensitivities
            (
                const volVectorField& faceSens,
                const labelList& sensitivityPatchIDs
            );

            //- Control point sensitivities computed using face-based surface
            //- sensitivities
            vectorField computeControlPointSensitivities
            (
                const boundaryVectorField& faceSens,
                const labelList& sensitivityPatchIDs
            );

            //- Control point sensitivities computed using face-based surface
            //- sensitivities
            vector computeControlPointSensitivities
            (
                const vectorField& faceSens,
                const label patchI,
                const label cpI
            );

            //- Part of control point sensitivities related to the face normal
            //- variations
            tmp<tensorField> dndbBasedSensitivities
            (
                const label patchI,
                const label cpI,
                bool DimensionedNormalSens = true
            );

            //- Get patch dx/db
            tmp<tensorField> patchDxDb
            (
                const label patchI,
                const label cpI
            );

            //- Get patch dx/db
            tmp<tensorField> patchDxDbFace
            (
                const label patchI,
                const label cpI
            );


        // Cartesian coordinates

            //- Compute cartesian coordinates based on control points
            //- and parametric coordinates
            tmp<vectorField> coordinates(const vectorField& uVector) const;

            //- The same, for a specific point
            vector coordinates(const vector& uVector) const;

            //- Mesh movement based on given control point movement
            tmp<vectorField> computeNewPoints
            (
                const vectorField& controlPointsMovement
            );

            //- Boundary mesh movement based on given control point movement
            tmp<vectorField> computeNewBoundaryPoints
            (
                const vectorField& controlPointsMovement,
                const labelList& patchesToBeMoved,
                const bool moveCPs = true
            );


        // Control points management

            //- Get control point ID from its I-J-K coordinates
            label getCPID(const label i, const label j, const label k) const;

            //- Get I-J-K coordinates from control point ID
            void getIJK(label& i, label& j, label& k, const label cpID) const;

            //- Set new control points
            //  New values should be on the coordinates system original CPs
            //  were defined
            void setControlPoints(const vectorField& newCps);


            //- Bound control points movement in the boundary control points
            //- and in certain directions if needed
            void boundControlPointMovement
            (
                vectorField& controlPointsMovement
            ) const;

            //- Compute max. displacement at the boundary
            scalar computeMaxBoundaryDisplacement
            (
                const vectorField& controlPointsMovement,
                const labelList& patchesToBeMoved
            );


        // Access Functions

            // Functions triggering calculations

                //- Get mesh points that reside within the control points box
                tmp<vectorField> getPointsInBox();

                //- Get map of points in box to mesh points
                const labelList& getMap();

                //- Get map of mesh points to points in box.
                //- Return -1 if point is outside the box
                const labelList& getReverseMap();

                //- Get parametric coordinates
                const pointVectorField& getParametricCoordinates();

                //- Get dxCartesiandb for a certain control point
                tmp<pointTensorField> getDxDb(const label cpI);

                //- Get dxCartesiandb for a certain control point on cells
                tmp<volTensorField> getDxCellsDb(const label cpI);

                //- Get number of variables if CPs are moved symmetrically in U
                label nUSymmetry() const;

                //- Get number of variables if CPs are moved symmetrically in V
                label nVSymmetry() const;

                //- Get number of variables if CPs are moved symmetrically in W
                label nWSymmetry() const;

                //- Get number of variables per direction,
                //- if CPs are moved symmetrically
                Vector<label> nSymmetry() const;


            // Inline access functions

                //- Get box name
                inline const word& name() const;

                //- Which control points are active?
                //  A control point is active if at least one component can move
                inline const boolList& getActiveCPs() const;

                //- Which design variables are active?
                //  Numbering is (X1,Y1,Z1), (X2,Y2,Z2) ...
                inline const boolList& getActiveDesignVariables() const;

                //- Get control points
                inline const vectorField& getControlPoints() const;

                inline vectorField& getControlPoints();

                //- Get confine movements
                inline bool confineUMovement() const;

                inline bool confineVMovement() const;

                inline bool confineWMovement() const;

                //- Get basis functions
                inline const NURBSbasis& basisU() const;

                inline const NURBSbasis& basisV() const;

                inline const NURBSbasis& basisW() const;

                //- Get number of control points per direction
                inline Vector<label> nCPsPerDirection() const;

                //- Get mesh
                inline const fvMesh& mesh() const;

                //- Get dictionary
                inline const dictionary& dict() const;


        // Write Functions

            //- Write control points on a cartesian coordinates system for
            //- visualization
            void writeCps
            (
                const fileName& = "cpsFile",
                const bool transform = true
            ) const;

            //- Write parametric coordinates
            void writeParamCoordinates() const;

            //- Write the control points to support restart
            virtual bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "NURBS3DVolumeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
