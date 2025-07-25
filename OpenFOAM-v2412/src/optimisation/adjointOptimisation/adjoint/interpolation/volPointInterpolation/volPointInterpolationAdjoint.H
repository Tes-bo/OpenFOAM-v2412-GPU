/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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
    Foam::volPointInterpolationAdjoint

Description
    Interpolate from cell centres to points (vertices) using inverse distance
    weighting

SourceFiles
    volPointInterpolationAdjoint.C
    volPointInterpolate.C

\*---------------------------------------------------------------------------*/

#ifndef volPointInterpolationAdjoint_H
#define volPointInterpolationAdjoint_H

#include "MeshObject.H"
#include "scalarList.H"
#include "volFields.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class fvMesh;
class pointMesh;

/*---------------------------------------------------------------------------*\
                       Class volPointInterpolationAdjoint Declaration
\*---------------------------------------------------------------------------*/

class volPointInterpolationAdjoint
:
    public MeshObject
    <
        fvMesh,
        UpdateableMeshObject,
        volPointInterpolationAdjoint
    >
{
    // Private Typedefs

        typedef MeshObject
        <
            fvMesh,
            UpdateableMeshObject,
            volPointInterpolationAdjoint
        > MeshObject_type;


protected:

    // Protected data

        //- Boundary addressing
        autoPtr<primitivePatch> boundaryPtr_;

        //- Per boundary face whether is on non-coupled, non-empty patch
        bitSet boundaryIsPatchFace_;

        //- Per mesh(!) point whether is on non-coupled, non-empty patch (on
        //  any processor)
        bitSet isPatchPoint_;

        //- Per mesh(!) point whether is on symmetry plane
        //  any processor)
        bitSet isSymmetryPoint_;

        //- Per boundary point the weights per pointFaces.
        scalarListList boundaryPointWeights_;


    // Protected Member Functions

        //- Construct addressing over all boundary faces
        void calcBoundaryAddressing();

        //- Make weights for points on uncoupled patches
        void makeBoundaryWeights(scalarField& sumWeights);

        //- Construct all point weighting factors
        void makeWeights();

        //- Helper: push master point data to collocated points
        template<class Type>
        void pushUntransformedData(List<Type>&) const;

        //- Get boundary field in same order as boundary faces. Field is
        //  zero on all coupled and empty patches
        template<class Type>
        tmp<Field<Type>> flatBoundaryField
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;

        //- Add separated contributions
        template<class Type>
        void addSeparated
        (
            GeometricField<Type, pointPatchField, pointMesh>&
        ) const;

        //- No copy construct
        volPointInterpolationAdjoint(const volPointInterpolationAdjoint&) = delete;

        //- No copy assignment
        void operator=(const volPointInterpolationAdjoint&) = delete;


public:

    // Declare name of the class and its debug switch
    ClassName("volPointInterpolationAdjoint");


    // Constructors

        //- Construct given fvMesh
        explicit volPointInterpolationAdjoint(const fvMesh&);


    //- Destructor
    ~volPointInterpolationAdjoint();


    // Member functions

    // Edit

        //- Update mesh topology using the morph engine
        void updateMesh(const mapPolyMesh&);

        //- Correct weighting factors for moving mesh.
        bool movePoints();


    // Interpolation Functions

        //- Interpolate sensitivties from points to faces
        //  conditions
        template<class Type>
        void interpolateSensitivitiesField
        (
            const GeometricField<Type, pointPatchField, pointMesh>& pf,
            typename GeometricField<Type, fvPatchField, volMesh>::Boundary& vf,
            const labelHashSet& patchIDs
        ) const;

        template<class Type>
        void interpolateBoundaryField
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf,
            GeometricField<Type, pointPatchField, pointMesh>& pf
        ) const;

        //- Interpolate sensitivties from points to faces
        void interpolateSensitivitiesField
        (
            const vectorField& pf,
            vectorField& vf,
            const labelHashSet& patchIDs
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "volPointInterpolateAdjoint.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
