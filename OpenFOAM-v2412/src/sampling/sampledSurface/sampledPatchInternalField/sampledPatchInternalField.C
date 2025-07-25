/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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
#include "dictionary.H"
#include "polyMesh.H"
#include "polyPatch.H"
#include "volFields.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sampledPatchInternalField, 0);
    addNamedToRunTimeSelectionTable
    (
        sampledSurface,
        sampledPatchInternalField,
        word,
        patchInternalField
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledPatchInternalField::sampledPatchInternalField
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledPatch(name, mesh, dict),
    mappers_(patchIDs().size())
{
    mappedPatchBase::offsetMode mode =
        mappedPatchBase::offsetModeNames_.getOrDefault
        (
            "offsetMode", dict, mappedPatchBase::NORMAL
        );

    switch (mode)
    {
        case mappedPatchBase::NORMAL:
        {
            const scalar distance(dict.get<scalar>("distance"));
            forAll(patchIDs(), i)
            {
                mappers_.set
                (
                    i,
                    new mappedPatchBase
                    (
                        mesh.boundaryMesh()[patchIDs()[i]],
                        mesh.name(),                        // sampleRegion
                        mappedPatchBase::NEARESTCELL,       // sampleMode
                        word::null,                         // samplePatch
                        -distance                  // sample inside my domain
                    )
                );
            }
        }
        break;

        case mappedPatchBase::UNIFORM:
        {
            const point offset(dict.get<point>("offset"));
            forAll(patchIDs(), i)
            {
                mappers_.set
                (
                    i,
                    new mappedPatchBase
                    (
                        mesh.boundaryMesh()[patchIDs()[i]],
                        mesh.name(),                        // sampleRegion
                        mappedPatchBase::NEARESTCELL,       // sampleMode
                        word::null,                         // samplePatch
                        offset                  // sample inside my domain
                    )
                );
            }
        }
        break;

        case mappedPatchBase::NONUNIFORM:
        {
            const pointField offsets(dict.get<pointField>("offsets"));
            forAll(patchIDs(), i)
            {
                mappers_.set
                (
                    i,
                    new mappedPatchBase
                    (
                        mesh.boundaryMesh()[patchIDs()[i]],
                        mesh.name(),                        // sampleRegion
                        mappedPatchBase::NEARESTCELL,       // sampleMode
                        word::null,                         // samplePatch
                        offsets                  // sample inside my domain
                    )
                );
            }
        }
        break;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::sampledPatchInternalField::sample
(
    const interpolation<scalar>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::vectorField> Foam::sampledPatchInternalField::sample
(
    const interpolation<vector>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::sphericalTensorField> Foam::sampledPatchInternalField::sample
(
    const interpolation<sphericalTensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::symmTensorField> Foam::sampledPatchInternalField::sample
(
    const interpolation<symmTensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::tensorField> Foam::sampledPatchInternalField::sample
(
    const interpolation<tensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::scalarField> Foam::sampledPatchInternalField::interpolate
(
    const interpolation<scalar>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::vectorField> Foam::sampledPatchInternalField::interpolate
(
    const interpolation<vector>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::sphericalTensorField>
Foam::sampledPatchInternalField::interpolate
(
    const interpolation<sphericalTensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::symmTensorField> Foam::sampledPatchInternalField::interpolate
(
    const interpolation<symmTensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::tensorField> Foam::sampledPatchInternalField::interpolate
(
    const interpolation<tensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


void Foam::sampledPatchInternalField::print(Ostream& os, int level) const
{
    os  << "sampledPatchInternalField: " << name() << " :"
        << " patches:" << patchNames();

    if (level)
    {
        os  << "  faces:" << faces().size()
            << "  points:" << points().size();
    }
}


// ************************************************************************* //
