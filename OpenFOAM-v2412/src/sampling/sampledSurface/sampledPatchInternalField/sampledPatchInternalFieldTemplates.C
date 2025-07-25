/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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

#include "sampledPatchInternalField.H"
#include "interpolationCellPoint.H"
#include "primitivePatchInterpolation.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledPatchInternalField::sampleOnFaces
(
    const interpolation<Type>& sampler
) const
{
    const auto& vField = sampler.psi();

    // One value per face
    auto tvalues = tmp<Field<Type>>::New(patchFaceLabels().size());
    auto& values = tvalues.ref();

    forAll(patchStart(), i)
    {
        // Get patchface wise data by sampling internal field
        Field<Type> interpVals = vField.primitiveField();
        mappers_[i].map().distribute(interpVals);

        // Store at correct position in values
        label end =
        (
            i < patchStart().size()-1
          ? patchStart()[i+1]
          : patchFaceLabels().size()
        );

        for (label triI = patchStart()[i]; triI < end; triI++)
        {
            values[triI] = interpVals[patchFaceLabels()[triI]];
        }
    }

    return tvalues;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledPatchInternalField::sampleOnPoints
(
    const interpolation<Type>& interpolator
) const
{
    label sz = 0;
    forAll(patchIDs(), i)
    {
        sz += mesh().boundaryMesh()[patchIDs()[i]].size();
    }

    Field<Type> allPatchVals(sz);
    sz = 0;

    forAll(patchIDs(), i)
    {
        // See mappedFixedValueFvPatchField
        const mapDistribute& distMap = mappers_[i].map();

        // Send back sample points to processor that holds the cell.
        // Mark cells with point::max so we know which ones we need
        // to interpolate (since expensive).
        vectorField samples(mappers_[i].samplePoints());
        distMap.reverseDistribute(mesh().nCells(), point::max, samples);

        Field<Type> patchVals(mesh().nCells());

        forAll(samples, celli)
        {
            if (samples[celli] != point::max)
            {
                patchVals[celli] = interpolator.interpolate
                (
                    samples[celli],
                    celli
                );
            }
        }

        distMap.distribute(patchVals);

        // Now patchVals holds the interpolated data in patch face order.
        // Collect.
        SubList<Type>(allPatchVals, patchVals.size(), sz) = patchVals;
        sz += patchVals.size();
    }

    // Interpolate to points. Reconstruct the patch of all faces to aid
    // interpolation.

    labelList meshFaceLabels(allPatchVals.size());
    sz = 0;
    for (const label patchId : patchIDs())
    {
        const polyPatch& pp = mesh().boundaryMesh()[patchId];
        forAll(pp, i)
        {
            meshFaceLabels[sz++] = pp.start()+i;
        }
    }

    indirectPrimitivePatch allPatches
    (
        IndirectList<face>(mesh().faces(), meshFaceLabels),
        mesh().points()
    );

    return PrimitivePatchInterpolation<indirectPrimitivePatch>
    (
        allPatches
    ).faceToPointInterpolate(allPatchVals);
}


// ************************************************************************* //
