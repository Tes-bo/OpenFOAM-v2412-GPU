/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2022 OpenCFD Ltd.
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

#include "sampledFaceZone.H"
#include "volFieldsFwd.H"
#include "pointFields.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledFaceZone::sampleOnFaces
(
    const interpolation<Type>& sampler
) const
{
    const auto& vField = sampler.psi();
    const labelList& own = mesh().faceOwner();
    const labelList& nei = mesh().faceNeighbour();

    // One value per face
    auto tvalues = tmp<Field<Type>>::New(faceId_.size());
    auto& values = tvalues.ref();

    forAll(faceId_, i)
    {
        const label facei = faceId_[i];
        const label patchi = facePatchId_[i];

        if (patchi != -1)
        {
            // Boundary face - face id is the patch-local face id
            values[i] = vField.boundaryField()[patchi][facei];
        }
        else
        {
            // Internal face
            values[i] = 0.5*(vField[own[facei]] + vField[nei[facei]]);
        }
    }

    return tvalues;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledFaceZone::sampleOnFaces
(
    const SurfaceField<Type>& sField
) const
{
    // One value per face
    auto tvalues = tmp<Field<Type>>::New(faceId_.size());
    auto& values = tvalues.ref();

    forAll(faceId_, i)
    {
        const label facei = faceId_[i];
        const label patchi = facePatchId_[i];

        if (patchi != -1)
        {
            // Boundary face - face id is the patch-local face id
            values[i] = sField.boundaryField()[patchi][facei];
        }
        else
        {
            // Internal face
            values[i] = sField[facei];
        }
    }

    return tvalues;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledFaceZone::sampleOnPoints
(
    const interpolation<Type>& interpolator
) const
{
    // One value per vertex
    auto tvalues = tmp<Field<Type>>::New(points().size(), Zero);
    auto& values = tvalues.ref();

    const labelList& own = mesh().faceOwner();

    bitSet pointDone(points().size());

    forAll(faces(), i)
    {
        const face& f = faces()[i];
        label facei = faceId_[i];
        const label patchi = facePatchId_[i];

        if (patchi != -1)
        {
            // Boundary face. patch-local face id -> mesh face id
            const polyPatch& pp = mesh().boundaryMesh()[patchi];

            facei += pp.start();
        }


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
