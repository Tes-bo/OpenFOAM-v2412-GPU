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
    Foam::volBSplinesBase

Description
    Class constructing a number of volumetric B-Splines boxes,
    read from dynamicMeshDict. Useful for various sensitivities and
    optMeshMovement classes.

    Derives from MeshObject so that all instances know and update the same
    control points and parametric coordinates are computed only once

SourceFiles
    volBSplinesBase.C

\*---------------------------------------------------------------------------*/

#ifndef volBSplinesBase_H
#define volBSplinesBase_H

#include "NURBS3DVolume.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class volBSplinesBase Declaration
\*---------------------------------------------------------------------------*/

class volBSplinesBase
:
    public MeshObject<fvMesh, UpdateableMeshObject, volBSplinesBase>
{
    // Private Typedefs

        typedef MeshObject
        <
            fvMesh,
            UpdateableMeshObject,
            volBSplinesBase
        > MeshObject_type;


    // Private Member Functions

        //- No copy construct
        volBSplinesBase(const volBSplinesBase&) = delete;

        //- No copy assignment
        void operator=(const volBSplinesBase&) = delete;


protected:

    // Protected Data

        //- List with volumetric B-splines boxes.
        //  No overlapping is supported
        PtrList<NURBS3DVolume> volume_;

        //- Active design variables numbering for all boxes
        labelList activeDesignVariables_;


public:

    //- Runtime type information
    TypeName("volBSplinesBase");


    // Constructors

        //- Construct from components
        explicit volBSplinesBase(const fvMesh& mesh);


    //- Destructor
    virtual ~volBSplinesBase() = default;


    // Member Functions

        //- Get const reference to the vol. B-splines boxes
        const PtrList<NURBS3DVolume>& boxes() const;

        //- Get non-const reference to the vol. B-splines boxes
        PtrList<NURBS3DVolume>& boxesRef();

        //- Get const reference to a specific box
        const NURBS3DVolume& box(const label boxI) const;

        //- Get non-const reference to a specific box
        NURBS3DVolume& boxRef(const label boxI);

        //- Get reference to control points
        const vectorField& getControlPoints(const label& iNURB) const;

        //- Get control points from all boxes
        vectorField getAllControlPoints() const;

        //- Get cumulative number of control points from all boxes
        label getTotalControlPointsNumber() const;

        //- Get number of boxes
        label getNumberOfBoxes() const;

        //- Get start CP ID for each box
        labelList getStartCpID() const;

        //- Get start CP ID for each box
        labelList getStartVarID() const;

        //- Find box of certain control point
        label findBoxID(const label cpI) const;

        //- From design variable ID, return boxID, cpID and direction
        Vector<label> decomposeDV(const label dvI) const;

        //- Get active design variables
        const labelList& getActiveDesignVariables() const;

        //- Get max boundary displacement for a given control-points
        //- movement
        scalar computeMaxBoundaryDisplacement
        (
            const vectorField& controlPointsMovement,
            const labelList& patchesToBeMoved
        );

        //- Get the updated boundary points only
        tmp<vectorField> computeBoundaryDisplacement
        (
            const vectorField& controlPointsMovement,
            const labelList& patchesToBeMoved
        );

        //- Bound control points movement
        void boundControlPointMovement
        (
            vectorField& controlPointsMovement
        ) const;

        //- Move control points. No effect on mesh
        void moveControlPoints(const vectorField& controlPointsMovement);

        //- Write control points to constant and optimisation folders
        void writeControlPoints() const;

        //- Dummy function required by MeshObject.
        //  Since this class is going to initiate the mesh movement,
        //  there is nothing more to be done when the mesh points change
        virtual bool movePoints();

        //- Dummy function required by MeshObject.
        //  Since this class is going to initiate the mesh movement,
        //  there is nothing more to be done when the mesh points change
        virtual void updateMesh(const mapPolyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
