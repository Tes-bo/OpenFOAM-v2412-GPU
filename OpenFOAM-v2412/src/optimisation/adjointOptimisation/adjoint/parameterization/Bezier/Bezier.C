/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019, 2022-2023 PCOpt/NTUA
    Copyright (C) 2013-2019, 2022-2023 FOSS GP
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

\*---------------------------------------------------------------------------*/

#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "Bezier.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(Bezier, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Bezier::Bezier(const fvMesh& mesh, const dictionary& dict)
:
    mesh_(mesh),
    dict_(dict),
    nBezier_(dict_.subDict("Bezier").get<label>("nBezier")),
    dxidXj_(nBezier_),
    confineXmovement_
    (
        dict_.subDict("Bezier").get<boolList>("confineXmovement")
    ),
    confineYmovement_
    (
        dict_.subDict("Bezier").get<boolList>("confineYmovement")
    ),
    confineZmovement_
    (
        dict_.subDict("Bezier").get<boolList>("confineZmovement")
    ),
    confineMovement_(3, boolList(nBezier_, false)),
    activeDesignVariables_(3*nBezier_)
{
    if
    (
        confineXmovement_.size() != nBezier_
     || confineYmovement_.size() != nBezier_
     || confineZmovement_.size() != nBezier_
    )
    {
        FatalErrorInFunction
            << "confineMovement lists sizes "
            << confineXmovement_.size()  << " "
            << confineYmovement_.size()  << " "
            << confineZmovement_.size()  << " "
            << "are incompatible with nBezier " << nBezier_
            << endl << endl
            << exit(FatalError);
    }
    confineMovement_[0] = confineXmovement_;
    confineMovement_[1] = confineYmovement_;
    confineMovement_[2] = confineZmovement_;

    // Determine active design variables
    label iActive(0);
    for (label iDir = 0; iDir < 3; ++iDir)
    {
        for (label iCP = 0; iCP < nBezier_; ++iCP)
        {
            if (!confineMovement_[iDir][iCP])
            {
                activeDesignVariables_[iActive++] = iDir*nBezier_ + iCP;
            }
        }
    }
    activeDesignVariables_.setSize(iActive);

    // Read bezier parameterization
    for (label iCP = 0; iCP < nBezier_; ++iCP)
    {
        dxidXj_.set
        (
            iCP,
            new pointTensorField
            (
                IOobject
                (
                    "dxidXj_"+name(iCP),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                pointMesh::New(mesh_)
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

label Bezier::nBezier() const
{
    return nBezier_;
}


PtrList<pointTensorField>& Bezier::dxidXj()
{
    return dxidXj_;
}


const boolList& Bezier::confineXmovement() const
{
    return confineXmovement_;
}


const boolList& Bezier::confineYmovement() const
{
    return confineYmovement_;
}


const boolList& Bezier::confineZmovement() const
{
    return confineZmovement_;
}


const boolListList& Bezier::confineMovement() const
{
    return confineMovement_;
}


tmp<tensorField> Bezier::dndbBasedSensitivities
(
    const label patchI,
    const label cpI,
    bool returnDimensionedNormalSens
) const
{
    const fvPatch& patch = mesh_.boundary()[patchI];
    const polyPatch& ppatch = patch.patch();

    // Return field
    auto tdndbSens = tmp<tensorField>::New(patch.size(), Zero);
    auto& dndbSens = tdndbSens.ref();

    // Auxiliary quantities
    deltaBoundary deltaBoundary(mesh_);
    const label patchStart = ppatch.start();
    const tensorField& dxdbInt = dxidXj_[cpI].primitiveField();

    // Loop over patch faces
    forAll(patch, fI)
    {
        const face& fGlobal = mesh_.faces()[fI + patchStart];
        const pointField facePoints = fGlobal.points(mesh_.points());

        // Loop over face points
        tensorField facePointDerivs(facePoints.size(), Zero);
        forAll(fGlobal, pI)
        {
            facePointDerivs[pI] = dxdbInt[fGlobal[pI]];
        }

        // Determine whether to return variance of dimensioned or unit normal
        if (returnDimensionedNormalSens)
        {
            dndbSens[fI] =
                deltaBoundary.makeFaceCentresAndAreas_d
                (
                    facePoints,
                    facePointDerivs
                )[1];
        }
        else
        {
            dndbSens[fI] =
                deltaBoundary.makeFaceCentresAndAreas_d
                (
                    facePoints,
                    facePointDerivs
                )[2];
        }
    }
    return tdndbSens;
}


tmp<vectorField> Bezier::dndbBasedSensitivities
(
    const label patchI,
    const label cpI,
    const label idir,
    bool returnDimensionedNormalSens
) const
{
    const fvPatch& patch = mesh_.boundary()[patchI];
    const polyPatch& ppatch = patch.patch();

    // Return field
    auto tdndbSens = tmp<vectorField>::New(patch.size(), Zero);
    auto& dndbSens = tdndbSens.ref();

    // Auxiliary quantities
    deltaBoundary deltaBoundary(mesh_);
    const label patchStart = ppatch.start();
    const tensorField& dxdbInt = dxidXj_[cpI].primitiveField();
    vectorField dxdbDir(dxdbInt.size(), Zero);
    unzipRow(dxdbInt, vector::components(idir), dxdbDir);

    // Loop over patch faces
    forAll(patch, fI)
    {
        const face& fGlobal = mesh_.faces()[fI + patchStart];
        const pointField facePoints = fGlobal.points(mesh_.points());

        // Loop over face points
        vectorField facePointDerivs(facePoints.size(), Zero);
        forAll(fGlobal, pI)
        {
            facePointDerivs[pI] = dxdbDir[fGlobal[pI]];
        }

        // Determine whether to return variance of dimensioned or unit normal
        if (returnDimensionedNormalSens)
        {
            dndbSens[fI] =
                deltaBoundary.makeFaceCentresAndAreas_d
                (
                    facePoints, facePointDerivs
                )[1];
        }
        else
        {
            dndbSens[fI] =
                deltaBoundary.makeFaceCentresAndAreas_d
                (
                    facePoints,
                    facePointDerivs
                )[2];
        }
    }

    return tdndbSens;
}


tmp<tensorField> Bezier::dxdbFace
(
    const label patchI,
    const label cpI,
    bool useChainRule
) const
{
    const polyPatch& patch = mesh_.boundary()[patchI].patch();

    // Return field
    auto tdxdbFace = tmp<tensorField>::New(patch.size(), Zero);
    auto& dxdbFace = tdxdbFace.ref();
    if (useChainRule)
    {
        // Auxiliary quantities
        deltaBoundary deltaBoundary(mesh_);
        const label patchStart = patch.start();
        const tensorField& dxdbInt = dxidXj_[cpI].primitiveField();

        // Loop over patch faces
        forAll(patch, fI)
        {
            const face& fGlobal = mesh_.faces()[fI + patchStart];
            const pointField facePoints = fGlobal.points(mesh_.points());

            // Loop over face points
            tensorField facePointDerivs(facePoints.size(), Zero);
            forAll(fGlobal, pI)
            {
                facePointDerivs[pI] = dxdbInt[fGlobal[pI]];
            }
            dxdbFace[fI] =
                deltaBoundary.makeFaceCentresAndAreas_d
                (
                    facePoints,
                    facePointDerivs
                )[0];
        }
    }
    // Simple average of dxdb in points. Older and less accurate
    else
    {
        PrimitivePatchInterpolation<polyPatch> patchInter(patch);
        dxdbFace = patchInter.pointToFaceInterpolate
        (
            dxidXj_[cpI].boundaryField()[patchI].patchInternalField()()
        )();
    }
    return tdxdbFace;
}


tmp<vectorField> Bezier::dxdbFace
(
    const label patchI,
    const label cpI,
    const label idir,
    bool useChainRule
) const
{
    const polyPatch& patch = mesh_.boundary()[patchI].patch();

    // Return field
    auto tdxdbFace = tmp<vectorField>::New(patch.size(), Zero);
    auto& dxdbFace = tdxdbFace.ref();

    if (useChainRule)
    {
        // Auxiliary quantities
        deltaBoundary deltaBoundary(mesh_);
        const label patchStart = patch.start();
        const tensorField& dxdbInt = dxidXj_[cpI].primitiveField();
        vectorField dxdbDir(dxdbInt.size(), Zero);
        dxdbDir.replace(0, dxdbInt.component(3*idir));
        dxdbDir.replace(1, dxdbInt.component(3*idir + 1));
        dxdbDir.replace(2, dxdbInt.component(3*idir + 2));

        // Loop over patch faces
        forAll(patch, fI)
        {
            const face& fGlobal = mesh_.faces()[fI + patchStart];
            const pointField facePoints = fGlobal.points(mesh_.points());

            // Loop over face points
            vectorField facePointDerivs(facePoints.size(), Zero);
            forAll(fGlobal, pI)
            {
                facePointDerivs[pI] = dxdbDir[fGlobal[pI]];
            }
            dxdbFace[fI] =
                deltaBoundary.makeFaceCentresAndAreas_d
                (
                    facePoints,
                    facePointDerivs
                )[0];
        }
    }
    // Simple average of dxdb in points. Older and less accurate
    else
    {
        PrimitivePatchInterpolation<polyPatch> patchInter(patch);
        vectorField dxdb(patch.nPoints(), Zero);
        dxdb.replace
        (
            idir,
            dxidXj_[cpI].boundaryField()[patchI].patchInternalField()().
                component(0)
        );
        dxdbFace = patchInter.pointToFaceInterpolate(dxdb)();
    }
    return tdxdbFace;
}


tensorField Bezier::facePoints_d
(
   const label globalFaceI,
   const label cpI
) const
{
    const face& faceI(mesh_.faces()[globalFaceI]);
    tensorField fPoints_d(faceI.size(), Zero);
    forAll(faceI, fpI)
    {
        fPoints_d[fpI] = dxidXj_[cpI].primitiveField()[faceI[fpI]];
    }
    return fPoints_d;
}


vectorField Bezier::facePoints_d
(
    const label globalFaceI,
    const label cpI,
    const label idir
) const
{
    const face& faceI(mesh_.faces()[globalFaceI]);
    vectorField fPoints_d(faceI.size(), Zero);
    forAll(faceI, fpI)
    {
        const tensor& dxdbTensor = dxidXj_[cpI].primitiveField()[faceI[fpI]];
        fPoints_d[fpI].x() = dxdbTensor.component(3*idir);
        fPoints_d[fpI].y() = dxdbTensor.component(3*idir + 1);
        fPoints_d[fpI].z() = dxdbTensor.component(3*idir + 2);
    }
    return fPoints_d;
}


const labelList& Bezier::getActiveDesignVariables() const
{
    return activeDesignVariables_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
