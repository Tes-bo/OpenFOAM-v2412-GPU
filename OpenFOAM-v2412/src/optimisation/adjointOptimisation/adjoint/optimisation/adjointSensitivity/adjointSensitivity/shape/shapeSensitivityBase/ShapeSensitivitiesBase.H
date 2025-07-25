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
    Foam::ShapeSensitivitiesBase

Description
    Base class supporting Shape sensitivity derivatives.

    Reference:
    \verbatim
        For the FI formulation see
            Kavvadias, I., Papoutsis-Kiachagias, E., & Giannakoglou, K. (2015).
            On the proper treatment of grid sensitivities in continuous adjoint
            methods for shape optimization.
            Journal of Computational Physics, 301, 1–18.
            http://doi.org/10.1016/j.jcp.2015.08.012

        The ESI formulation is derived in a slightly different way than the
        one described in this paper, to provide a common mathematical
        formulation for both low- and high-Re meshes and to produce numerically
        identical results as the FI formulation. In brief, the boundary-bound
        part of the sensitivities is the patchInternalField of the tensor
        multiplying grad(dxdb) in the FI formulation.

    \endverbatim

SourceFiles
    ShapeSensitivitiesBase.C

\*---------------------------------------------------------------------------*/

#ifndef ShapeSensitivitiesBase_H
#define ShapeSensitivitiesBase_H

#include "adjointSensitivity.H"
#include "objective.H"
#include "volPointInterpolation.H"
#include "pointMesh.H"
#include "pointPatchField.H"
#include "pointPatchFieldsFwd.H"
#include "fixedValuePointPatchField.H"
#include "boundaryFieldsFwd.H"
#include "createZeroField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class ShapeSensitivitiesBase Declaration
\*---------------------------------------------------------------------------*/

class ShapeSensitivitiesBase
    :
    public adjointSensitivity
{

protected:

    // Protected data

        //- Patches on which to compute shape sensitivities
        labelHashSet sensitivityPatchIDs_;

        //- Whether to write all surface sensitivity fields
        bool writeAllSurfaceFiles_;

        // autoPtrs for fields holding sensitivities.
        // Not all of them are required for each case

        // Boundary sensitivities at faces. Shape opt & flow control
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        //- Wall face sens w.r.t. (x,y.z)
        autoPtr<boundaryVectorField> wallFaceSensVecPtr_;

        //- Wall face sens projected to normal
        autoPtr<boundaryScalarField> wallFaceSensNormalPtr_;

        //- Normal sens as vectors
        autoPtr<boundaryVectorField> wallFaceSensNormalVecPtr_;

        // Boundary sensitivities at points. Shape opt
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        //- Wall point sens w.r.t. (x,y.z)
        autoPtr<pointBoundaryVectorField> wallPointSensVecPtr_;

        //- Wall point sens projected to normal
        autoPtr<pointBoundaryScalarField> wallPointSensNormalPtr_;

        //- Normal sens as vectors
        autoPtr<pointBoundaryVectorField> wallPointSensNormalVecPtr_;


    // Protected Member Functions

        //- Allocate the adjoint eikonal solver
        void allocateEikonalSolver();

        //- Check if any of the available objective has a certain multiplier,
        //- provided through a function object
        bool hasMultiplier(bool (objective::*hasFunction)() const);

        //- Constructs volField based on boundaryField and writes it
        template<class Type>
        void constructAndWriteSensitivityField
        (
            const autoPtr
            <
                typename GeometricField<Type, fvPatchField, volMesh>::Boundary
            >& sensFieldPtr,
            const word& name
        ) const;

        //- Constructs pointField based on boundaryField and writes it
        template<class Type>
        void constructAndWriteSensitivtyPointField
        (
            const autoPtr<List<Field<Type>>>& sensFieldPtr,
            const word& name
        ) const;

        //- Constructs volField based on boundaryField and writes it
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh>>
        constructVolSensitivtyField
        (
            const autoPtr
            <
                typename GeometricField<Type, fvPatchField, volMesh>::Boundary
            >& sensFieldPtr,
            const word& name
        ) const;

        //- Write face-based sensitivities, if present
        void writeFaceBasedSens() const;

        //- Write point-based sensitivities, if present
        void writePointBasedSens() const;

        //- Clear surface/point fields
        void clearSurfaceFields();

        //- Allocate multiplier fields
        void allocateMultipliers();

        //- Clear multipliers
        void clearMultipliers();


private:

    // Private Member Functions

        //- No copy construct
        ShapeSensitivitiesBase(const ShapeSensitivitiesBase&) = delete;

        //- No copy assignment
        void operator=(const ShapeSensitivitiesBase&) = delete;


public:

    //- Runtime type information
    TypeName("ShapeSensitivitiesBase");


    // Constructors

        //- Construct from components
        ShapeSensitivitiesBase
        (
            const fvMesh& mesh,
            const dictionary& dict,
            adjointSolver& adjointSolver
        );


    //- Destructor
    virtual ~ShapeSensitivitiesBase() = default;


    // Member Functions

        //- Read dict if changed
        virtual bool readDict(const dictionary& dict);

        //- Get patch IDs on which sensitivities are computed
        inline const labelHashSet& sensitivityPatchIDs() const
        {
            return sensitivityPatchIDs_;
        }

        //- Overwrite sensitivityPatchIDs
        inline void setSensitivityPatchIDs(const labelHashSet& sensPatchIDs)
        {
            sensitivityPatchIDs_ = sensPatchIDs;
        }

        //- Return set of patches on which to compute direct sensitivities
        virtual const labelHashSet& geometryVariationIntegrationPatches() const;

        //- Accumulate sensitivity integrands
        //  Common function for the FI and E-SI approaches
        virtual void accumulateIntegrand(const scalar dt);

        //- Zero sensitivity fields and their constituents
        void clearSensitivities();

        //- Write sensitivity fields.
        //  If valid, copies boundaryFields to volFields and writes them.
        virtual void write(const word& baseName = word::null);

        //- Get wall face sensitivity vectors field
        tmp<volVectorField> getWallFaceSensVec();

        //- Get wall face sensitivity projected to normal field
        tmp<volScalarField> getWallFaceSensNormal();

        //- Get wall face normal sens as vectors field
        tmp<volVectorField> getWallFaceSensNormalVec();

        //- Get wall point sensitivity vectors field
        //  Uses volPointInterpolation
        tmp<pointVectorField> getWallPointSensVec();

        //- Get wall point sensitivity projected to normal field
        //  Uses volPointInterpolation
        tmp<pointScalarField> getWallPointSensNormal();

        //- Get wall point sens as vectors field
        //  Uses volPointInterpolation
        tmp<pointVectorField> getWallPointSensNormalVec();

        //- Get wall face sensitivity vectors field
        virtual const boundaryVectorField& getWallFaceSensVecBoundary() const;

        //- Get wall face sensitivity projected to normal field
        virtual const boundaryScalarField&
            getWallFaceSensNormalBoundary() const;

        //- Get wall face normal sens as vectors field
        virtual const boundaryVectorField&
            getWallFaceSensNormalVecBoundary() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ShapeSensitivitiesBaseTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
