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

\*---------------------------------------------------------------------------*/

#include "fvcSurfaceOps.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class ResultType, class CellToFaceOp>
void surfaceSum
(
    const surfaceScalarField& lambdas,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const CellToFaceOp& cop,
    GeometricField<ResultType, fvPatchField, volMesh>& result,
    const bool doCorrectBoundaryConditions
)
{
    const fvMesh& mesh = vf.mesh();
    const auto& Sf = mesh.Sf();
    const auto& P = mesh.owner();
    const auto& N = mesh.neighbour();

    const auto& vfi = vf.primitiveField();
    auto& sfi = result.primitiveFieldRef();

    // See e.g. surfaceInterpolationScheme<Type>::dotInterpolate

    // Internal field
    {
        const auto& Sfi = Sf.primitiveField();
        const auto& lambda = lambdas.primitiveField();

        for (label facei=0; facei<P.size(); facei++)
        {
            const label ownCelli = P[facei];
            const label neiCelli = N[facei];

            const ResultType faceVal
            (
                cop
                (
                    Sfi[facei],
                    lambda[facei],
                    vfi[ownCelli],
                    vfi[neiCelli]
                )
            );
            sfi[ownCelli] += faceVal;
            sfi[neiCelli] -= faceVal;
        }
    }


    // Boundary field
    {
        forAll(mesh.boundary(), patchi)
        {
            const auto& pFaceCells = mesh.boundary()[patchi].faceCells();
            const auto& pSf = Sf.boundaryField()[patchi];
            const auto& pvf = vf.boundaryField()[patchi];
            const auto& pLambda = lambdas.boundaryField()[patchi];

            if (pvf.coupled())
            {
                auto tpnf(pvf.patchNeighbourField());
                auto& pnf = tpnf();

                for (label facei=0; facei<pFaceCells.size(); facei++)
                {
                    // Interpolate between owner-side and neighbour-side values
                    const ResultType faceVal
                    (
                        cop
                        (
                            pSf[facei],
                            pLambda[facei],
                            vfi[pFaceCells[facei]],
                            pnf[facei]
                        )
                    );

                    sfi[pFaceCells[facei]] += faceVal;
                }
            }
            else
            {
                for (label facei=0; facei<pFaceCells.size(); facei++)
                {
                    // Use patch value only
                    const ResultType faceVal
                    (
                        cop
                        (
                            pSf[facei],
                            scalar(1.0),
                            pvf[facei],
                            pTraits<Type>::zero  // not used
                        )
                    );
                    sfi[pFaceCells[facei]] += faceVal;
                }
            }
        }
    }

    if (doCorrectBoundaryConditions)
    {
        result.correctBoundaryConditions();
    }
}


template<class Type, class FType, class ResultType, class CellToFaceOp>
void surfaceSum
(
    const surfaceScalarField& lambdas,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const GeometricField<FType, fvsPatchField, surfaceMesh>& sadd,
    const CellToFaceOp& cop,
    GeometricField<ResultType, fvPatchField, volMesh>& result,
    const bool doCorrectBoundaryConditions
)
{
    const fvMesh& mesh = vf.mesh();
    const auto& Sf = mesh.Sf();
    const auto& P = mesh.owner();
    const auto& N = mesh.neighbour();

    const auto& vfi = vf.primitiveField();
    auto& sfi = result.primitiveFieldRef();

    // See e.g. surfaceInterpolationScheme<Type>::dotInterpolate

    // Internal field
    {
        const auto& Sfi = Sf.primitiveField();
        const auto& lambda = lambdas.primitiveField();
        const auto& saddi = sadd.primitiveField();

        for (label facei=0; facei<P.size(); facei++)
        {
            const label ownCelli = P[facei];
            const label neiCelli = N[facei];

            const ResultType faceVal
            (
                cop
                (
                    Sfi[facei],
                    lambda[facei],
                    vfi[ownCelli],
                    vfi[neiCelli],

                    saddi[facei]        // additional face value
                )
            );
            sfi[ownCelli] += faceVal;
            sfi[neiCelli] -= faceVal;
        }
    }


    // Boundary field
    {
        forAll(mesh.boundary(), patchi)
        {
            const auto& pFaceCells = mesh.boundary()[patchi].faceCells();
            const auto& pSf = Sf.boundaryField()[patchi];
            const auto& pvf = vf.boundaryField()[patchi];
            const auto& pLambda = lambdas.boundaryField()[patchi];
            const auto& psadd = sadd.boundaryField()[patchi];

            if (pvf.coupled())
            {
                auto tpnf(pvf.patchNeighbourField());
                auto& pnf = tpnf();

                for (label facei=0; facei<pFaceCells.size(); facei++)
                {
                    // Interpolate between owner-side and neighbour-side values
                    const ResultType faceVal
                    (
                        cop
                        (
                            pSf[facei],
                            pLambda[facei],
                            vfi[pFaceCells[facei]],
                            pnf[facei],
                            psadd[facei]
                        )
                    );

                    sfi[pFaceCells[facei]] += faceVal;
                }
            }
            else
            {
                for (label facei=0; facei<pFaceCells.size(); facei++)
                {
                    // Use patch value only
                    const ResultType faceVal
                    (
                        cop
                        (
                            pSf[facei],
                            scalar(1.0),
                            pvf[facei],
                            pTraits<Type>::zero,  // not used

                            psadd[facei]
                        )
                    );
                    sfi[pFaceCells[facei]] += faceVal;
                }
            }
        }
    }

    if (doCorrectBoundaryConditions)
    {
        result.correctBoundaryConditions();
    }
}


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
    const bool doCorrectBoundaryConditions
)
{
    const fvMesh& mesh = vf.mesh();
    const auto& Sf = mesh.Sf();
    const auto& P = mesh.owner();
    const auto& N = mesh.neighbour();

    const auto& vfi = vf.primitiveField();
    auto& resulti = result.primitiveFieldRef();

    // See e.g. surfaceInterpolationScheme<Type>::dotInterpolate

    // Internal field
    {
        const auto& Sfi = Sf.primitiveField();
        const auto& lambda = lambdas.primitiveField();
        const auto& sf0i = sf0.primitiveField();
        const auto& sf1i = sf1.primitiveField();

        for (label facei=0; facei<P.size(); facei++)
        {
            const label ownCelli = P[facei];
            const label neiCelli = N[facei];

            const ResultType faceVal
            (
                cop
                (
                    Sfi[facei],

                    lambda[facei],
                    vfi[ownCelli],
                    vfi[neiCelli],

                    sf0i[facei],        // additional face value
                    sf1i[facei]         // additional face value
                )
            );
            resulti[ownCelli] += faceVal;
            resulti[neiCelli] -= faceVal;
        }
    }


    // Boundary field
    {
        forAll(mesh.boundary(), patchi)
        {
            const auto& pFaceCells = mesh.boundary()[patchi].faceCells();
            const auto& pSf = Sf.boundaryField()[patchi];
            const auto& pvf = vf.boundaryField()[patchi];
            const auto& pLambda = lambdas.boundaryField()[patchi];
            const auto& psf0 = sf0.boundaryField()[patchi];
            const auto& psf1 = sf1.boundaryField()[patchi];

            if (pvf.coupled())
            {
                auto tpnf(pvf.patchNeighbourField());
                auto& pnf = tpnf();

                for (label facei=0; facei<pFaceCells.size(); facei++)
                {
                    // Interpolate between owner-side and neighbour-side values
                    const ResultType faceVal
                    (
                        cop
                        (
                            pSf[facei],

                            pLambda[facei],
                            vfi[pFaceCells[facei]],
                            pnf[facei],

                            psf0[facei],
                            psf1[facei]
                        )
                    );

                    resulti[pFaceCells[facei]] += faceVal;
                }
            }
            else
            {
                for (label facei=0; facei<pFaceCells.size(); facei++)
                {
                    // Use patch value only
                    const ResultType faceVal
                    (
                        cop
                        (
                            pSf[facei],
                            scalar(1.0),
                            pvf[facei],
                            pTraits<Type>::zero,  // not used

                            psf0[facei],
                            psf1[facei]
                        )
                    );
                    resulti[pFaceCells[facei]] += faceVal;
                }
            }
        }
    }

    if (doCorrectBoundaryConditions)
    {
        result.correctBoundaryConditions();
    }
}


template<class Type, class ResultType, class CombineOp>
void GaussOp
(
    const surfaceScalarField& lambdas,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const CombineOp& cop,
    GeometricField<ResultType, fvPatchField, volMesh>& result
)
{
    const fvMesh& mesh = vf.mesh();

    // Sum contributions
    surfaceSum(lambdas, vf, cop, result, false);

    auto& sfi = result.primitiveFieldRef();
    sfi /= mesh.V();

    result.correctBoundaryConditions();
}


template<class Type, class ResultType, class CombineOp>
void surfaceOp
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const surfaceVectorField& ownLs,
    const surfaceVectorField& neiLs,
    const CombineOp& cop,
    GeometricField<ResultType, fvPatchField, volMesh>& result
)
{
    const fvMesh& mesh = vf.mesh();
    const auto& Sf = mesh.Sf();
    const auto& P = mesh.owner();
    const auto& N = mesh.neighbour();

    const auto& vfi = vf.primitiveField();
    auto& sfi = result.primitiveFieldRef();

    // Internal field
    {
        const auto& Sfi = Sf.primitiveField();

        for (label facei=0; facei<P.size(); facei++)
        {
            const label ownCelli = P[facei];
            const label neiCelli = N[facei];

            const Type faceVal
            (
                cop
                (
                    Sfi[facei],     // needed?
                    vfi[ownCelli],
                    vfi[neiCelli]
                )
            );
            sfi[ownCelli] += ownLs[facei]*faceVal;
            sfi[neiCelli] -= neiLs[facei]*faceVal;
        }
    }


    // Boundary field
    {
        forAll(mesh.boundary(), patchi)
        {
            const auto& pFaceCells = mesh.boundary()[patchi].faceCells();
            const auto& pSf = Sf.boundaryField()[patchi];
            const auto& pvf = vf.boundaryField()[patchi];
            const auto& pOwnLs = ownLs.boundaryField()[patchi];

            if (pvf.coupled())
            {
                auto tpnf(pvf.patchNeighbourField());
                auto& pnf = tpnf();

                for (label facei=0; facei<pFaceCells.size(); facei++)
                {
                    const Type faceVal
                    (
                        cop
                        (
                            pSf[facei], // needed?
                            vfi[pFaceCells[facei]],
                            pnf[facei]
                        )
                    );

                    sfi[pFaceCells[facei]] += pOwnLs[facei]*faceVal;
                }
            }
            else
            {
                for (label facei=0; facei<pFaceCells.size(); facei++)
                {
                    const Type faceVal
                    (
                        cop
                        (
                            pSf[facei],
                            vfi[pFaceCells[facei]],
                            pvf[facei]
                        )
                    );
                    sfi[pFaceCells[facei]] += pOwnLs[facei]*faceVal;
                }
            }
        }
    }

    result.correctBoundaryConditions();
}


template<class Type, class ResultType, class CellToFaceOp>
void surfaceSnSum
(
    const surfaceScalarField& deltaCoeffs,
    const GeometricField<Type, fvPatchField, volMesh>& vf,

    const CellToFaceOp& cop,

    GeometricField<ResultType, fvPatchField, volMesh>& result,
    const bool doCorrectBoundaryConditions
)
{
    const fvMesh& mesh = vf.mesh();
    const auto& Sf = mesh.Sf();
    const auto& P = mesh.owner();
    const auto& N = mesh.neighbour();

    const auto& vfi = vf.primitiveField();
    auto& sfi = result.primitiveFieldRef();

    // Internal field
    {
        const auto& Sfi = Sf.primitiveField();
        const auto& dc = deltaCoeffs.primitiveField();

        for (label facei=0; facei<P.size(); facei++)
        {
            const label ownCelli = P[facei];
            const label neiCelli = N[facei];

            const ResultType faceVal
            (
                cop
                (
                    Sfi[facei],         // area vector
                    dc[facei],          // delta coefficients
                    vfi[ownCelli],
                    vfi[neiCelli]
                )
            );
            sfi[ownCelli] += faceVal;
            sfi[neiCelli] -= faceVal;
        }
    }


    // Boundary field
    {
        forAll(mesh.boundary(), patchi)
        {
            const auto& pFaceCells = mesh.boundary()[patchi].faceCells();
            const auto& pSf = Sf.boundaryField()[patchi];
            const auto& pvf = vf.boundaryField()[patchi];
            const auto& pdc = deltaCoeffs.boundaryField()[patchi];

            if (pvf.coupled())
            {
                auto tpnf(pvf.patchNeighbourField());
                auto& pnf = tpnf();

                for (label facei=0; facei<pFaceCells.size(); facei++)
                {
                    const label ownCelli = pFaceCells[facei];

                    // Interpolate between owner-side and neighbour-side values
                    const ResultType faceVal
                    (
                        cop
                        (
                            pSf[facei],         // area vector
                            pdc[facei],
                            vfi[ownCelli],
                            pnf[facei]
                        )
                    );

                    sfi[ownCelli] += faceVal;
                }
            }
            else
            {
                auto tpnf(pvf.snGrad());
                auto& pnf = tpnf();

                for (label facei=0; facei<pFaceCells.size(); facei++)
                {
                    const label ownCelli = pFaceCells[facei];

                    // Use patch value only
                    const ResultType faceVal
                    (
                        cop
                        (
                            pSf[facei],             // area vector
                            scalar(1.0),            // use 100% of pnf
                            pTraits<Type>::zero,
                            pnf[facei]
                        )
                    );
                    sfi[ownCelli] += faceVal;
                }
            }
        }
    }

    if (doCorrectBoundaryConditions)
    {
        result.correctBoundaryConditions();
    }
}


template<class Type, class ResultType, class CellToFaceOp>
void surfaceSnSum
(
    const surfaceScalarField& deltaCoeffs,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sadd,

    const CellToFaceOp& cop,

    GeometricField<ResultType, fvPatchField, volMesh>& result,
    const bool doCorrectBoundaryConditions
)
{
    const fvMesh& mesh = vf.mesh();
    const auto& Sf = mesh.Sf();
    const auto& P = mesh.owner();
    const auto& N = mesh.neighbour();

    const auto& vfi = vf.primitiveField();
    auto& sfi = result.primitiveFieldRef();

    // Internal field
    {
        const auto& Sfi = Sf.primitiveField();
        const auto& dc = deltaCoeffs.primitiveField();
        const auto& saddi = sadd.primitiveField();

        for (label facei=0; facei<P.size(); facei++)
        {
            const label ownCelli = P[facei];
            const label neiCelli = N[facei];

            const ResultType faceVal
            (
                cop
                (
                    Sfi[facei],         // area vector
                    dc[facei],          // delta coefficients
                    vfi[ownCelli],
                    vfi[neiCelli],
                    saddi[facei]        // face addition value
                )
            );
            sfi[ownCelli] += faceVal;
            sfi[neiCelli] -= faceVal;
        }
    }


    // Boundary field
    {
        forAll(mesh.boundary(), patchi)
        {
            const auto& pFaceCells = mesh.boundary()[patchi].faceCells();
            const auto& pSf = Sf.boundaryField()[patchi];
            const auto& pvf = vf.boundaryField()[patchi];
            const auto& pdc = deltaCoeffs.boundaryField()[patchi];
            const auto& psadd = sadd.boundaryField()[patchi];

            if (pvf.coupled())
            {
                auto tpnf(pvf.patchNeighbourField());
                auto& pnf = tpnf();

                for (label facei=0; facei<pFaceCells.size(); facei++)
                {
                    const label ownCelli = pFaceCells[facei];

                    // Interpolate between owner-side and neighbour-side values
                    const ResultType faceVal
                    (
                        cop
                        (
                            pSf[facei],         // area vector
                            pdc[facei],
                            vfi[ownCelli],
                            pnf[facei],
                            psadd[facei]
                        )
                    );

                    sfi[ownCelli] += faceVal;
                }
            }
            else
            {
                auto tpnf(pvf.snGrad());
                auto& pnf = tpnf();

                for (label facei=0; facei<pFaceCells.size(); facei++)
                {
                    const label ownCelli = pFaceCells[facei];

                    // Use patch value only
                    const ResultType faceVal
                    (
                        cop
                        (
                            pSf[facei],             // area vector
                            scalar(1.0),            // use 100% of pnf
                            pTraits<Type>::zero,
                            pnf[facei],
                            psadd[facei]
                        )
                    );
                    sfi[ownCelli] += faceVal;
                }
            }
        }
    }

    if (doCorrectBoundaryConditions)
    {
        result.correctBoundaryConditions();
    }
}


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
)
{
    const fvMesh& mesh = vf.mesh();
    const auto& Sf = mesh.Sf();
    const auto& P = mesh.owner();
    const auto& N = mesh.neighbour();

    const auto& gammai = gamma.primitiveField();
    const auto& vfi = vf.primitiveField();
    auto& sfi = result.primitiveFieldRef();

    // Internal field
    {
        const auto& Sfi = Sf.primitiveField();
        const auto& weights = gammaWeights.primitiveField();
        const auto& dc = deltaCoeffs.primitiveField();

        for (label facei=0; facei<P.size(); facei++)
        {
            const label ownCelli = P[facei];
            const label neiCelli = N[facei];

            const ResultType faceVal
            (
                cop
                (
                    Sfi[facei],         // area vector

                    weights[facei],     // interpolation weights
                    gammai[ownCelli],
                    gammai[neiCelli],

                    dc[facei],          // delta coefficients
                    vfi[ownCelli],
                    vfi[neiCelli]
                )
            );
            sfi[ownCelli] += faceVal;
            sfi[neiCelli] -= faceVal;
        }
    }


    // Boundary field
    {
        forAll(mesh.boundary(), patchi)
        {
            const auto& pFaceCells = mesh.boundary()[patchi].faceCells();
            const auto& pSf = Sf.boundaryField()[patchi];
            const auto& pvf = vf.boundaryField()[patchi];
            const auto& pdc = deltaCoeffs.boundaryField()[patchi];
            const auto& pweights = gammaWeights.boundaryField()[patchi];
            const auto& pgamma = gamma.boundaryField()[patchi];

            if (pvf.coupled())
            {
                auto tpnf(pvf.patchNeighbourField());
                auto& pnf = tpnf();
                auto tgammanf(pgamma.patchNeighbourField());
                auto& gammanf = tgammanf();

                for (label facei=0; facei<pFaceCells.size(); facei++)
                {
                    const label ownCelli = pFaceCells[facei];

                    // Interpolate between owner-side and neighbour-side values
                    const ResultType faceVal
                    (
                        cop
                        (
                            pSf[facei],         // area vector

                            pweights[facei],
                            gammai[ownCelli],
                            gammanf[facei],

                            pdc[facei],
                            vfi[ownCelli],
                            pnf[facei]
                        )
                    );

                    sfi[ownCelli] += faceVal;
                }
            }
            else
            {
                auto tpnf(pvf.snGrad());
                auto& pnf = tpnf();

                for (label facei=0; facei<pFaceCells.size(); facei++)
                {
                    const label ownCelli = pFaceCells[facei];

                    // Use patch value only
                    const ResultType faceVal
                    (
                        cop
                        (
                            pSf[facei],             // area vector

                            scalar(1.0),
                            pgamma[facei],
                            pTraits<GType>::zero,   // not used


                            scalar(1.0),            // use 100% of pnf
                            pTraits<Type>::zero,
                            pnf[facei]
                        )
                    );
                    sfi[ownCelli] += faceVal;
                }
            }
        }
    }

    if (doCorrectBoundaryConditions)
    {
        result.correctBoundaryConditions();
    }
}


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
)
{
    const fvMesh& mesh = vf.mesh();
    const auto& Sf = mesh.Sf();
    const auto& P = mesh.owner();
    const auto& N = mesh.neighbour();

    const auto& gammai = gamma.primitiveField();
    const auto& vfi = vf.primitiveField();
    auto& sfi = result.primitiveFieldRef();

    // Internal field
    {
        const auto& Sfi = Sf.primitiveField();
        const auto& weights = gammaWeights.primitiveField();
        const auto& dc = deltaCoeffs.primitiveField();
        const auto& saddi = sadd.primitiveField();

        for (label facei=0; facei<P.size(); facei++)
        {
            const label ownCelli = P[facei];
            const label neiCelli = N[facei];

            const ResultType faceVal
            (
                cop
                (
                    Sfi[facei],         // area vector

                    weights[facei],     // interpolation weights
                    gammai[ownCelli],
                    gammai[neiCelli],

                    dc[facei],          // delta coefficients
                    vfi[ownCelli],
                    vfi[neiCelli],

                    saddi[facei]        // face addition value
                )
            );
            sfi[ownCelli] += faceVal;
            sfi[neiCelli] -= faceVal;
        }
    }


    // Boundary field
    {
        forAll(mesh.boundary(), patchi)
        {
            const auto& pFaceCells = mesh.boundary()[patchi].faceCells();
            const auto& pSf = Sf.boundaryField()[patchi];
            const auto& pvf = vf.boundaryField()[patchi];
            const auto& pdc = deltaCoeffs.boundaryField()[patchi];
            const auto& pweights = gammaWeights.boundaryField()[patchi];
            const auto& pgamma = gamma.boundaryField()[patchi];
            const auto& psadd = sadd.boundaryField()[patchi];

            if (pvf.coupled())
            {
                auto tpnf(pvf.patchNeighbourField());
                auto& pnf = tpnf();
                auto tgammanf(pgamma.patchNeighbourField());
                auto& gammanf = tgammanf();

                for (label facei=0; facei<pFaceCells.size(); facei++)
                {
                    const label ownCelli = pFaceCells[facei];

                    // Interpolate between owner-side and neighbour-side values
                    const ResultType faceVal
                    (
                        cop
                        (
                            pSf[facei],         // area vector

                            pweights[facei],
                            gammai[ownCelli],
                            gammanf[facei],

                            pdc[facei],
                            vfi[ownCelli],
                            pnf[facei],

                            psadd[facei]
                        )
                    );

                    sfi[ownCelli] += faceVal;
                }
            }
            else
            {
                auto tpnf(pvf.snGrad());
                auto& pnf = tpnf();

                for (label facei=0; facei<pFaceCells.size(); facei++)
                {
                    const label ownCelli = pFaceCells[facei];

                    // Use patch value only
                    const ResultType faceVal
                    (
                        cop
                        (
                            pSf[facei],             // area vector

                            scalar(1.0),
                            pgamma[facei],
                            pTraits<GType>::zero,   // not used


                            scalar(1.0),            // use 100% of pnf
                            pTraits<Type>::zero,
                            pnf[facei],

                            psadd[facei]
                        )
                    );
                    sfi[ownCelli] += faceVal;
                }
            }
        }
    }

    if (doCorrectBoundaryConditions)
    {
        result.correctBoundaryConditions();
    }
}


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
    const surfaceScalarField& gammaWeights,
    const GeometricField<GType0, fvPatchField, volMesh>& gamma0,
    const GeometricField<GType1, fvPatchField, volMesh>& gamma1,

    const surfaceScalarField& deltaCoeffs,
    const GeometricField<Type, fvPatchField, volMesh>& vf,

    const CellToFaceOp& cop,

    GeometricField<ResultType, fvPatchField, volMesh>& result,
    const bool doCorrectBoundaryConditions
)
{
    const fvMesh& mesh = vf.mesh();
    const auto& Sf = mesh.Sf();
    const auto& P = mesh.owner();
    const auto& N = mesh.neighbour();

    const auto& gamma0i = gamma0.primitiveField();
    const auto& gamma1i = gamma1.primitiveField();
    const auto& vfi = vf.primitiveField();
    auto& sfi = result.primitiveFieldRef();

    // Internal field
    {
        const auto& Sfi = Sf.primitiveField();
        const auto& weights = gammaWeights.primitiveField();
        const auto& dc = deltaCoeffs.primitiveField();

        for (label facei=0; facei<P.size(); facei++)
        {
            const label ownCelli = P[facei];
            const label neiCelli = N[facei];

            const ResultType faceVal
            (
                cop
                (
                    Sfi[facei],         // area vector

                    weights[facei],     // interpolation weights

                    gamma0i[ownCelli],
                    gamma0i[neiCelli],

                    gamma1i[ownCelli],
                    gamma1i[neiCelli],

                    dc[facei],          // delta coefficients
                    vfi[ownCelli],
                    vfi[neiCelli]
                )
            );
            sfi[ownCelli] += faceVal;
            sfi[neiCelli] -= faceVal;
        }
    }


    // Boundary field
    {
        forAll(mesh.boundary(), patchi)
        {
            const auto& pFaceCells = mesh.boundary()[patchi].faceCells();
            const auto& pSf = Sf.boundaryField()[patchi];
            const auto& pvf = vf.boundaryField()[patchi];
            const auto& pdc = deltaCoeffs.boundaryField()[patchi];
            const auto& pweights = gammaWeights.boundaryField()[patchi];
            const auto& pgamma0 = gamma0.boundaryField()[patchi];
            const auto& pgamma1 = gamma1.boundaryField()[patchi];

            if (pvf.coupled())
            {
                auto tpnf(pvf.patchNeighbourField());
                auto& pnf = tpnf();
                auto tgamma0nf(pgamma0.patchNeighbourField());
                auto& gamma0nf = tgamma0nf();
                auto tgamma1nf(pgamma1.patchNeighbourField());
                auto& gamma1nf = tgamma1nf();

                for (label facei=0; facei<pFaceCells.size(); facei++)
                {
                    const label ownCelli = pFaceCells[facei];

                    // Interpolate between owner-side and neighbour-side values
                    const ResultType faceVal
                    (
                        cop
                        (
                            pSf[facei],         // area vector

                            pweights[facei],

                            gamma0i[ownCelli],
                            gamma0nf[facei],

                            gamma1i[ownCelli],
                            gamma1nf[facei],

                            pdc[facei],
                            vfi[ownCelli],
                            pnf[facei]
                        )
                    );

                    sfi[ownCelli] += faceVal;
                }
            }
            else
            {
                auto tpnf(pvf.snGrad());
                auto& pnf = tpnf();

                for (label facei=0; facei<pFaceCells.size(); facei++)
                {
                    const label ownCelli = pFaceCells[facei];

                    // Use patch value only
                    const ResultType faceVal
                    (
                        cop
                        (
                            pSf[facei],             // area vector

                            scalar(1.0),

                            pgamma0[facei],
                            pTraits<GType0>::zero,  // not used

                            pgamma1[facei],
                            pTraits<GType1>::zero,  // not used

                            scalar(1.0),            // use 100% of pnf
                            pTraits<Type>::zero,
                            pnf[facei]
                        )
                    );
                    sfi[ownCelli] += faceVal;
                }
            }
        }
    }

    if (doCorrectBoundaryConditions)
    {
        result.correctBoundaryConditions();
    }
}


template<class Type, class FType, class ResultType, class CellToFaceOp>
void interpolate
(
    const surfaceScalarField& weights,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const GeometricField<FType, fvsPatchField, surfaceMesh>& sf,
    const CellToFaceOp& cop,
    GeometricField<ResultType, fvsPatchField, surfaceMesh>& result
)
{
    const fvMesh& mesh = vf.mesh();
    const auto& Sf = mesh.Sf();
    const auto& P = mesh.owner();
    const auto& N = mesh.neighbour();

    const auto& vfi = vf.primitiveField();

    // Internal field
    {
        const auto& Sfi = Sf.primitiveField();
        const auto& weight = weights.primitiveField();
        const auto& sfi = sf.primitiveField();

        auto& resulti = result.primitiveFieldRef();

        for (label facei=0; facei<P.size(); facei++)
        {
            const label ownCelli = P[facei];
            const label neiCelli = N[facei];

            cop
            (
                Sfi[facei],

                weight[facei],
                vfi[ownCelli],
                vfi[neiCelli],

                sfi[facei],

                resulti[facei]
            );
        }
    }


    // Boundary field
    {
        forAll(mesh.boundary(), patchi)
        {
            const auto& pFaceCells = mesh.boundary()[patchi].faceCells();
            const auto& pSf = Sf.boundaryField()[patchi];
            const auto& pvf = vf.boundaryField()[patchi];
            const auto& pweight = weights.boundaryField()[patchi];
            const auto& psf = sf.boundaryField()[patchi];
            auto& presult = result.boundaryFieldRef()[patchi];

            if (pvf.coupled())
            {
                auto tpnf(pvf.patchNeighbourField());
                auto& pnf = tpnf();

                for (label facei=0; facei<pFaceCells.size(); facei++)
                {
                    // Interpolate between owner-side and neighbour-side values
                    cop
                    (
                        pSf[facei],

                        pweight[facei],
                        vfi[pFaceCells[facei]],
                        pnf[facei],

                        psf[facei],

                        presult[facei]
                    );
                }
            }
            else
            {
                for (label facei=0; facei<pFaceCells.size(); facei++)
                {
                    // Use patch value only
                    cop
                    (
                        pSf[facei],

                        scalar(1.0),
                        pvf[facei],
                        pTraits<Type>::zero,    // not used

                        psf[facei],

                        presult[facei]
                    );
                }
            }
        }
    }
}
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
)
{
    const fvMesh& mesh = vf0.mesh();
    const auto& Sf = mesh.Sf();
    const auto& P = mesh.owner();
    const auto& N = mesh.neighbour();

    const auto& vf0i = vf0.primitiveField();
    const auto& vf1i = vf1.primitiveField();

    // Internal field
    {
        const auto& Sfi = Sf.primitiveField();
        const auto& weight = weights.primitiveField();

        auto& resulti = result.primitiveFieldRef();

        for (label facei=0; facei<P.size(); facei++)
        {
            const label ownCelli = P[facei];
            const label neiCelli = N[facei];

            cop
            (
                Sfi[facei],

                weight[facei],

                vf0i[ownCelli],
                vf0i[neiCelli],

                vf1i[ownCelli],
                vf1i[neiCelli],

                resulti[facei]
            );
        }
    }


    // Boundary field
    {
        forAll(mesh.boundary(), patchi)
        {
            const auto& pFaceCells = mesh.boundary()[patchi].faceCells();
            const auto& pSf = Sf.boundaryField()[patchi];
            const auto& pvf0 = vf0.boundaryField()[patchi];
            const auto& pvf1 = vf1.boundaryField()[patchi];
            const auto& pweight = weights.boundaryField()[patchi];
            auto& presult = result.boundaryFieldRef()[patchi];

            if (pvf0.coupled() || pvf1.coupled())
            {
                auto tpnf0(pvf0.patchNeighbourField());
                auto& pnf0 = tpnf0();

                auto tpnf1(pvf1.patchNeighbourField());
                auto& pnf1 = tpnf1();

                for (label facei=0; facei<pFaceCells.size(); facei++)
                {
                    // Interpolate between owner-side and neighbour-side values
                    cop
                    (
                        pSf[facei],

                        pweight[facei],

                        vf0i[pFaceCells[facei]],
                        pnf0[facei],

                        vf1i[pFaceCells[facei]],
                        pnf1[facei],

                        presult[facei]
                    );
                }
            }
            else
            {
                for (label facei=0; facei<pFaceCells.size(); facei++)
                {
                    // Use patch value only
                    cop
                    (
                        pSf[facei],

                        scalar(1.0),

                        pvf0[facei],
                        pTraits<Type0>::zero, // not used

                        pvf1[facei],
                        pTraits<Type1>::zero, // not used

                        presult[facei]
                    );
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
