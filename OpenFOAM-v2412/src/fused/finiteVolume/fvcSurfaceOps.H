/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 M.Janssens
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

InNamespace
    Foam::fvc

Description
    Surface integrate surfaceField creating a volField.
    Surface sum a surfaceField creating a volField.

SourceFiles
    fvcSurfaceOps.C

\*---------------------------------------------------------------------------*/


#ifndef fvcSurfaceOps_H
#define fvcSurfaceOps_H

#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- In-place operations on Fields. Add to FieldFunctions.C ?
#define INPLACE_PRODUCT_OPERATOR(product, CombineOp, Op, OpFunc)               \
                                                                               \
template<class Type1, class Type2>                                             \
void OpFunc                                                                    \
(                                                                              \
    Field<typename product<Type1, Type2>::type>& result,                       \
    const UList<Type1>& f1,                                                    \
    const UList<Type2>& f2                                                     \
)                                                                              \
{                                                                              \
    typedef typename product<Type1, Type2>::type resultType;                   \
    TFOR_ALL_F_OP_F_OP_F                                                       \
    (resultType, result, CombineOp, Type1, f1, Op, Type2, f2)                  \
}


INPLACE_PRODUCT_OPERATOR(outerProduct, +=, *, multiplyAdd)
INPLACE_PRODUCT_OPERATOR(outerProduct, -=, *, multiplySubtract)

#undef INPLACE_PRODUCT_OPERATOR



/*---------------------------------------------------------------------------*\
                     Namespace fvc functions Declaration
\*---------------------------------------------------------------------------*/

namespace fvc
{
    // Interpolation

        //- Interpolate to face (using cop) and additional face field
        template<class Type, class FType, class ResultType, class CellToFaceOp>
        void interpolate
        (
            const surfaceScalarField& lambdas,
            const GeometricField<Type, fvPatchField, volMesh>& vf,
            const GeometricField<FType, fvsPatchField, surfaceMesh>& sf,
            const CellToFaceOp& cop,
            GeometricField<ResultType, fvsPatchField, surfaceMesh>& result
        );

        //- Interpolate to face (using cop)
        template
        <
            class Type0,
            class Type1,
            class ResultType,
            class CellToFaceOp
        >
        void interpolate
        (
            const surfaceScalarField& weights,
            const GeometricField<Type0, fvPatchField, volMesh>& vf0,
            const GeometricField<Type1, fvPatchField, volMesh>& vf1,
            const CellToFaceOp& cop,
            GeometricField<ResultType, fvsPatchField, surfaceMesh>& result
        );


    // Interpolation and accumulation

        //- Interpolate to face (using cop) and accumulate.
        template<class Type, class ResultType, class CellToFaceOp>
        void surfaceSum
        (
            const surfaceScalarField& lambdas,
            const GeometricField<Type, fvPatchField, volMesh>& vf,
            const CellToFaceOp& cop,
            GeometricField<ResultType, fvPatchField, volMesh>& result,
            const bool doCorrectBoundaryConditions = true
        );

        //- Interpolate to face (using cop) and accumulate. Additional
        //- face field
        template<class Type, class FType, class ResultType, class CellToFaceOp>
        void surfaceSum
        (
            const surfaceScalarField& lambdas,
            const GeometricField<Type, fvPatchField, volMesh>& vf,
            const GeometricField<FType, fvsPatchField, surfaceMesh>& sf,
            const CellToFaceOp& cop,
            GeometricField<ResultType, fvPatchField, volMesh>& result,
            const bool doCorrectBoundaryConditions = true
        );

        //- Interpolate to face (using cop) and accumulate. Additional
        //- face fields
        template
        <
            class Type,
            class FType0,
            class FType1,
            class ResultType,
            class CellToFaceOp
        >
        void surfaceSum
        (
            const surfaceScalarField& lambdas,
            const GeometricField<Type, fvPatchField, volMesh>& vf,

            const GeometricField<FType0, fvsPatchField, surfaceMesh>& sf0,
            const GeometricField<FType1, fvsPatchField, surfaceMesh>& sf1,

            const CellToFaceOp& cop,

            GeometricField<ResultType, fvPatchField, volMesh>& result,
            const bool doCorrectBoundaryConditions = true
        );

        //- Interpolate to face (using cop) and apply Gauss. Note: uses V(),
        //  not Vsc()
        template<class Type, class ResultType, class CellToFaceOp>
        void GaussOp
        (
            const surfaceScalarField& lambdas,
            const GeometricField<Type, fvPatchField, volMesh>& vf,
            const CellToFaceOp& cop,
            GeometricField<ResultType, fvPatchField, volMesh>& result
        );


    // Difference and accumulation

        //- sum of snGrad
        template<class Type, class ResultType, class CellToFaceOp>
        void surfaceSnSum
        (
            const surfaceScalarField& deltaCoeffs,
            const GeometricField<Type, fvPatchField, volMesh>& vf,

            const CellToFaceOp& cop,

            GeometricField<ResultType, fvPatchField, volMesh>& result,
            const bool doCorrectBoundaryConditions
        );
        //- sum of snGrad with additional surface field
        template<class Type, class ResultType, class CellToFaceOp>
        void surfaceSnSum
        (
            const surfaceScalarField& deltaCoeffs,
            const GeometricField<Type, fvPatchField, volMesh>& vf,
            const GeometricField<Type, fvsPatchField, surfaceMesh>& sadd,

            const CellToFaceOp& cop,

            GeometricField<ResultType, fvPatchField, volMesh>& result,
            const bool doCorrectBoundaryConditions
        );

        //- sum of snGrad with additional (interpolated) volField
        template<class Type, class GType, class ResultType, class CellToFaceOp>
        void surfaceSnSum
        (
            const surfaceScalarField& gammaWeights,
            const GeometricField<GType, fvPatchField, volMesh>& gamma,

            const surfaceScalarField& deltaCoeffs,
            const GeometricField<Type, fvPatchField, volMesh>& vf,

            const CellToFaceOp& cop,

            GeometricField<ResultType, fvPatchField, volMesh>& result,
            const bool doCorrectBoundaryConditions
        );

        //- sum of snGrad with additional (interpolated) volfields
        template
        <
            class Type,
            class GType0,
            class GType1,
            class ResultType,
            class CellToFaceOp
        >
        void surfaceSnSum
        (
            const surfaceScalarField& weights,

            const GeometricField<GType0, fvPatchField, volMesh>& gamma0,
            const GeometricField<GType1, fvPatchField, volMesh>& gamma1,

            const surfaceScalarField& deltaCoeffs,
            const GeometricField<Type, fvPatchField, volMesh>& vf,

            const CellToFaceOp& cop,

            GeometricField<ResultType, fvPatchField, volMesh>& result,
            const bool doCorrectBoundaryConditions
        );

        //- sum of snGrad with additional surface field
        template<class Type, class GType, class ResultType, class CellToFaceOp>
        void surfaceSnSum
        (
            const surfaceScalarField& gammaWeights,
            const GeometricField<GType, fvPatchField, volMesh>& gamma,

            const surfaceScalarField& deltaCoeffs,
            const GeometricField<Type, fvPatchField, volMesh>& vf,

            const GeometricField<Type, fvsPatchField, surfaceMesh>& sadd,

            const CellToFaceOp& cop,

            GeometricField<ResultType, fvPatchField, volMesh>& result,
            const bool doCorrectBoundaryConditions
        );


    // Other

        //- Interpolate to face (using cop) and apply distribution vectors
        template<class Type, class ResultType, class CellToFaceOp>
        void surfaceOp
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf,
            const surfaceVectorField& ownLs,
            const surfaceVectorField& neiLs,
            const CellToFaceOp& cop,
            GeometricField<ResultType, fvPatchField, volMesh>& result
        );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvcSurfaceOps.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
