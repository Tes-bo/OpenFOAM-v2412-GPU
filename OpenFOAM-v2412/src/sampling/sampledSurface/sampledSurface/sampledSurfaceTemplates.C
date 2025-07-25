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

#include "sampledSurface.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledSurface::sampleOnFaces
(
    const interpolation<Type>& sampler,
    const labelUList& elements,
    const faceList& fcs,
    const pointField& pts,
    const Type& defaultValue
)
{
    const label len = elements.size();

    if (len != fcs.size())
    {
        FatalErrorInFunction
            << "size mismatch: "
            << "sampled elements (" << len
            << ") != faces (" << fcs.size() << ')'
            << exit(FatalError);
    }

    auto tvalues = tmp<Field<Type>>::New(len);
    auto& values = tvalues.ref();

    for (label i=0; i < len; ++i)
    {
        const label celli = elements[i];
        if (celli < 0)
        {
            values[i] = defaultValue;
        }
        else
        {
            const point pt = fcs[i].centre(pts);

            values[i] = sampler.interpolate(pt, celli);
        }
    }

    return tvalues;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledSurface::sampleOnPoints
(
    const interpolation<Type>& interpolator,
    const labelUList& elements,
    const faceList& fcs,
    const pointField& pts
)
{
    const label len = elements.size();

    if (len != fcs.size())
    {
        FatalErrorInFunction
            << "size mismatch: "
            << "sampled elements (" << len
            << ") != faces (" << fcs.size() << ')'
            << exit(FatalError);
    }

    // One value per point
    // Initialize with Zero to handle missed/degenerate faces
    auto tvalues = tmp<Field<Type>>::New(pts.size(), Zero);
    auto& values = tvalues.ref();

    bitSet pointDone(pts.size());

    forAll(fcs, facei)
    {
        const face& f = fcs[facei];
        const label celli = elements[facei];

        for (const label pointi : f)
        {
            if (pointDone.set(pointi))
            {
                values[pointi] = interpolator.interpolate
                (
                    pts[pointi],
                    celli
                );
            }
        }
    }

    return tvalues;
}


template<class Type>
Foam::tmp<Foam::VolumeField<Type>>
Foam::sampledSurface::pointAverage
(
    const PointField<Type>& pfld
)
{
    const fvMesh& mesh = dynamic_cast<const fvMesh&>(pfld.mesh()());

    auto tcellAvg = tmp<VolumeField<Type>>::New
    (
        IOobject
        (
            "cellAvg",
            mesh.time().timeName(),
            pfld.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh,
        dimensioned<Type>(dimless, Zero)
    );
    auto& cellAvg = tcellAvg.ref();

    labelField nPointCells(mesh.nCells(), Zero);

    for (label pointi = 0; pointi < mesh.nPoints(); ++pointi)
    {
        const Type& val = pfld[pointi];
        const labelList& pCells = mesh.pointCells(pointi);

        for (const label celli : pCells)
        {
            cellAvg[celli] += val;
            ++nPointCells[celli];
        }
    }

    forAll(cellAvg, celli)
    {
        cellAvg[celli] /= nPointCells[celli];
    }

    // Give value to calculatedFvPatchFields
    cellAvg.correctBoundaryConditions();

    return tcellAvg;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class GeoMeshType>
bool Foam::sampledSurface::storeRegistryField
(
    const objectRegistry& obr,
    const word& fieldName,
    const dimensionSet& dims,
    const Field<Type>& values,
    word lookupName
) const
{
    polySurface* surfptr = this->getRegistrySurface(obr, lookupName);

    if (surfptr)
    {
        surfptr->storeField<Type, GeoMeshType>
        (
            fieldName, dims, values
        );
    }

    return surfptr;
}


template<class Type, class GeoMeshType>
bool Foam::sampledSurface::storeRegistryField
(
    const objectRegistry& obr,
    const word& fieldName,
    const dimensionSet& dims,
    Field<Type>&& values,
    word lookupName
) const
{
    polySurface* surfptr = this->getRegistrySurface(obr, lookupName);

    if (surfptr)
    {
        surfptr->storeField<Type, GeoMeshType>
        (
            fieldName, dims, std::move(values)
        );
    }

    return surfptr;
}


// ************************************************************************* //
