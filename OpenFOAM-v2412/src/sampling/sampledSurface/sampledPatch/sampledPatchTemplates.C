/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

#include "sampledPatch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledPatch::sampleOnFaces
(
    const interpolation<Type>& sampler
) const
{
    const auto& vField = sampler.psi();

    // One value per face
    auto tvalues = tmp<Field<Type>>::New(patchFaceLabels_.size());
    auto& values = tvalues.ref();

    forAll(patchFaceLabels_, i)
    {
        const label patchi = patchIDs_[patchIndex_[i]];
        const label patchFacei = patchFaceLabels_[i];

        values[i] = vField.boundaryField()[patchi][patchFacei];
    }

    return tvalues;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledPatch::sampleOnFaces
(
    const SurfaceField<Type>& sField
) const
{
    // One value per face
    auto tvalues = tmp<Field<Type>>::New(patchFaceLabels_.size());
    auto& values = tvalues.ref();

    forAll(patchFaceLabels_, i)
    {
        const label patchi = patchIDs_[patchIndex_[i]];
        const label patchFacei = patchFaceLabels_[i];

        values[i] = sField.boundaryField()[patchi][patchFacei];
    }

    return tvalues;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledPatch::sampleOnPoints
(
    const interpolation<Type>& interpolator
) const
{
    // One value per vertex
    auto tvalues = tmp<Field<Type>>::New(points().size());
    auto& values = tvalues.ref();

    const labelList& own = mesh().faceOwner();

    bitSet pointDone(points().size());

    forAll(faces(), i)
    {
        const face& f = faces()[i];
        const label patchi = patchIDs_[patchIndex_[i]];
        const label patchFacei = patchFaceLabels_[i];

        const polyPatch& pp = mesh().boundaryMesh()[patchi];

        const label facei = patchFacei + pp.start();
        const label celli = own[facei];

        for (const label pointi : f)
        {
            if (pointDone.set(pointi))
            {
                values[pointi] = interpolator.interpolate
                (
                    points()[pointi],
                    celli,
                    facei
                );
            }
        }
    }

    return tvalues;
}


// ************************************************************************* //
