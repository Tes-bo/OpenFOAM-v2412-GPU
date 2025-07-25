/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
    Copyright (C) 2020 Henning Scheufler
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

#include "fvcSmooth.H"
#include "volFields.H"
#include "FaceCellWave.H"
#include "smoothData.H"
#include "sweepData.H"
#include "fvMatrices.H"
#include "fvcVolumeIntegrate.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"
#include "zeroGradientFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::fvc::smooth
(
    volScalarField& field,
    const scalar coeff
)
{
    const fvMesh& mesh = field.mesh();
    scalar maxRatio = 1 + coeff;

    DynamicList<label> changedFaces(mesh.nFaces()/100 + 100);
    DynamicList<smoothData> changedFacesInfo(changedFaces.size());

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nbr = neighbour[facei];

        // Check if owner value much larger than neighbour value or vice versa
        if (field[own] > maxRatio*field[nbr])
        {
            changedFaces.append(facei);
            changedFacesInfo.append(smoothData(field[own]));
        }
        else if (field[nbr] > maxRatio*field[own])
        {
            changedFaces.append(facei);
            changedFacesInfo.append(smoothData(field[nbr]));
        }
    }

    // Insert all faces of coupled patches - FaceCellWave will correct them
    forAll(mesh.boundaryMesh(), patchi)
    {
        const polyPatch& patch = mesh.boundaryMesh()[patchi];

        if (patch.coupled())
        {
            forAll(patch, patchFacei)
            {
                label facei = patch.start() + patchFacei;
                label own = mesh.faceOwner()[facei];

                changedFaces.append(facei);
                changedFacesInfo.append(smoothData(field[own]));
            }
        }
    }

    changedFaces.shrink();
    changedFacesInfo.shrink();

    // Set initial field on cells
    List<smoothData> cellData(mesh.nCells());

    forAll(field, celli)
    {
        cellData[celli] = field[celli];
    }

    // Set initial field on faces
    List<smoothData> faceData(mesh.nFaces());

    smoothData::trackData td;
    td.maxRatio = maxRatio;

    // Propagate information over whole domain
    FaceCellWave<smoothData, smoothData::trackData> smoothData
    (
        mesh,
        changedFaces,
        changedFacesInfo,
        faceData,
        cellData,
        mesh.globalData().nTotalCells(),   // max iterations
        td
    );

    forAll(field, celli)
    {
        field[celli] = cellData[celli].value();
    }

    field.correctBoundaryConditions();
}


void Foam::fvc::spread
(
    volScalarField& field,
    const volScalarField& alpha,
    const label nLayers,
    const scalar alphaDiff,
    const scalar alphaMax,
    const scalar alphaMin
)
{
    const fvMesh& mesh = field.mesh();

    DynamicList<label> changedFaces(mesh.nFaces()/100 + 100);
    DynamicList<smoothData> changedFacesInfo(changedFaces.size());

    // Set initial field on cells
    List<smoothData> cellData(mesh.nCells());

    forAll(field, celli)
    {
        cellData[celli] = field[celli];
    }

    // Set initial field on faces
    List<smoothData> faceData(mesh.nFaces());

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nbr = neighbour[facei];

        if (mag(alpha[own] - alpha[nbr]) > alphaDiff)
        {
            changedFaces.append(facei);
            changedFacesInfo.append
            (
                smoothData(max(field[own], field[nbr]))
            );
        }
    }

    // Insert all faces of coupled patches - FaceCellWave will correct them
    forAll(mesh.boundaryMesh(), patchi)
    {
        const polyPatch& patch = mesh.boundaryMesh()[patchi];

        if (patch.coupled())
        {
            forAll(patch, patchFacei)
            {
                label facei = patch.start() + patchFacei;
                label own = mesh.faceOwner()[facei];

                const scalarField alphapn
                (
                    alpha.boundaryField()[patchi].patchNeighbourField()
                );

                if (mag(alpha[own] - alphapn[patchFacei]) > alphaDiff)
                {
                    changedFaces.append(facei);
                    changedFacesInfo.append(smoothData(field[own]));
                }
            }
        }
    }

    changedFaces.shrink();
    changedFacesInfo.shrink();

    smoothData::trackData td;
    td.maxRatio = 1.0;

    // Propagate information over whole domain
    FaceCellWave<smoothData, smoothData::trackData> smoothData
    (
        mesh,
        faceData,
        cellData,
        td
    );

    smoothData.setFaceInfo(changedFaces, changedFacesInfo);

    smoothData.iterate(nLayers);

    forAll(field, celli)
    {
        field[celli] = cellData[celli].value();
    }

    field.correctBoundaryConditions();
}


void Foam::fvc::sweep
(
    volScalarField& field,
    const volScalarField& alpha,
    const label nLayers,
    const scalar alphaDiff
)
{
    const fvMesh& mesh = field.mesh();

    DynamicList<label> changedFaces(mesh.nFaces()/100 + 100);
    DynamicList<sweepData> changedFacesInfo(changedFaces.size());

    // Set initial field on cells
    List<sweepData> cellData(mesh.nCells());

    // Set initial field on faces
    List<sweepData> faceData(mesh.nFaces());

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const vectorField& Cf = mesh.faceCentres();

    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nbr = neighbour[facei];

        if (mag(alpha[own] - alpha[nbr]) > alphaDiff)
        {
            changedFaces.append(facei);
            changedFacesInfo.append
            (
                sweepData(max(field[own], field[nbr]), Cf[facei])
            );
        }
    }

    // Insert all faces of coupled patches - FaceCellWave will correct them
    forAll(mesh.boundaryMesh(), patchi)
    {
        const polyPatch& patch = mesh.boundaryMesh()[patchi];

        if (patch.coupled())
        {
            forAll(patch, patchFacei)
            {
                label facei = patch.start() + patchFacei;
                label own = mesh.faceOwner()[facei];

                const scalarField alphapn
                (
                    alpha.boundaryField()[patchi].patchNeighbourField()
                );

                if (mag(alpha[own] - alphapn[patchFacei]) > alphaDiff)
                {
                    changedFaces.append(facei);
                    changedFacesInfo.append
                    (
                        sweepData(field[own], Cf[facei])
                    );
                }
            }
        }
    }

    changedFaces.shrink();
    changedFacesInfo.shrink();

    // Propagate information over whole domain
    FaceCellWave<sweepData> sweepData
    (
        mesh,
        faceData,
        cellData
    );

    sweepData.setFaceInfo(changedFaces, changedFacesInfo);

    sweepData.iterate(nLayers);

    forAll(field, celli)
    {
        if (cellData[celli].valid(sweepData.data()))
        {
            field[celli] = max(field[celli], cellData[celli].value());
        }
    }

    field.correctBoundaryConditions();
}


void Foam::fvc::spreadSource
(
    volScalarField& mDotOut,
    const volScalarField& mDotIn,
    const volScalarField& alpha1,
    const volScalarField& alpha2,
    const dimensionedScalar& D,
    const scalar cutoff
)
{
    const fvMesh& mesh = alpha1.mesh();

    volScalarField mDotSmear
    (
        IOobject
        (
            "mDotSmear",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(mDotOut.dimensions(), Zero),
        fvPatchFieldBase::zeroGradientType()
    );

    // Smearing of source term field
    fvScalarMatrix mSourceEqn
    (
        fvm::Sp(scalar(1), mDotSmear)
      - fvm::laplacian(D, mDotSmear)
      ==
        mDotIn
    );

    mSourceEqn.solve();

    // Cut cells with cutoff < alpha1 < 1-cutoff and rescale remaining
    // source term field

    dimensionedScalar intvDotLiquid("intvDotLiquid", dimMass/dimTime, 0.0);
    dimensionedScalar intvDotVapor ("intvDotVapor", dimMass/dimTime, 0.0);

    const scalarField& Vol = mesh.V();

    forAll(mesh.C(), celli)
    {
        if (alpha1[celli] < cutoff)
        {
            intvDotVapor.value() +=
                alpha2[celli]*mDotSmear[celli]*Vol[celli];
        }
        else if (alpha1[celli] > 1.0 - cutoff)
        {
            intvDotLiquid.value() +=
                alpha1[celli]*mDotSmear[celli]*Vol[celli];
        }
    }

    reduce(intvDotVapor.value(), sumOp<scalar>());
    reduce(intvDotLiquid.value(), sumOp<scalar>());

    // Calculate Nl and Nv
    dimensionedScalar Nl ("Nl", dimless, Zero);
    dimensionedScalar Nv ("Nv", dimless, Zero);

    const dimensionedScalar intmSource0(fvc::domainIntegrate(mDotIn));

    if (intvDotVapor.value() > VSMALL)
    {
        Nv = intmSource0/intvDotVapor;
    }
    if (intvDotLiquid.value() > VSMALL)
    {
        Nl = intmSource0/intvDotLiquid;
    }

    // Set source terms in cells with alpha1 < cutoff or alpha1 > 1-cutoff
    forAll(mesh.C(), celli)
    {
        if (alpha1[celli] < cutoff)
        {
            mDotOut[celli] = Nv.value()*(1 - alpha1[celli])*mDotSmear[celli];
        }
        else if (alpha1[celli] > 1.0 - cutoff)
        {
            //mDotOut[celli] = 0;
            mDotOut[celli] = -Nl.value()*alpha1[celli]*mDotSmear[celli];
        }
        else
        {
            mDotOut[celli] = 0;
        }
    }
}


// ************************************************************************* //
