/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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

#include "fvMeshSubsetProxy.H"
#include "volFields.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Type>
Foam::tmp
<
    Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>
>
Foam::fvMeshSubsetProxy::zeroGradientField
(
    const DimensionedField<Type, volMesh>& df
)
{
    IOobject io(df);
    io.readOpt(IOobject::NO_READ);
    io.writeOpt(IOobject::NO_WRITE);
    io.registerObject(false);

    auto tfield = tmp<GeometricField<Type, fvPatchField, volMesh>>::New
    (
        io,
        df.mesh(),
        dimensioned<Type>(df.dimensions(), Zero),
        fvPatchFieldBase::zeroGradientType()
    );
    tfield.ref().primitiveFieldRef() = df;
    tfield.ref().oriented() = df.oriented();
    tfield.ref().correctBoundaryConditions();

    return tfield;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::fvMeshSubsetProxy::interpolateInternal
(
    const fvMeshSubset& subsetter,
    const DimensionedField<Type, volMesh>& df
)
{
    auto tfield = zeroGradientField<Type>(df);

    if (subsetter.hasSubMesh())
    {
        return interpolate(subsetter, tfield());
    }

    return tfield;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::fvMeshSubsetProxy::interpolateInternal
(
    const fvMeshSubset& subsetter,
    const tmp<DimensionedField<Type, volMesh>>& tdf
)
{
    // TODO - move dimensioned mesh into internal,
    // but needs different GeometricField constructors

    if (tdf.valid())
    {
        if (subsetter.hasSubMesh())
        {
            auto tproxied = interpolate(subsetter, tdf);
            auto tfield = zeroGradientField<Type>(tproxied());

            tdf.clear();
            tproxied.clear();
            return tfield;
        }
        else
        {
            auto tfield = zeroGradientField<Type>(tdf());

            tdf.clear();
            return tfield;
        }
    }

    return nullptr;
}


template<class GeoField>
Foam::tmp<GeoField>
Foam::fvMeshSubsetProxy::interpolate
(
    const fvMeshSubset& subsetter,
    const GeoField& fld
)
{
    if (subsetter.hasSubMesh())
    {
        auto tfield = subsetter.interpolate(fld);

        tfield.ref().checkOut();
        tfield.ref().rename(fld.name());
        return tfield;
    }

    return fld;
}


template<class GeoField>
Foam::tmp<GeoField>
Foam::fvMeshSubsetProxy::interpolate
(
    const fvMeshSubset& subsetter,
    const tmp<GeoField>& tfield
)
{
    if (tfield.valid() && subsetter.hasSubMesh())
    {
        auto tproxied = interpolate(subsetter, tfield());
        tfield.clear();

        return tproxied;
    }

    // Nothing to be done
    return tfield;
}


// * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::fvMeshSubsetProxy::interpolateInternal
(
    const DimensionedField<Type, volMesh>& df
) const
{
    return interpolateInternal(subsetter_, df);
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::fvMeshSubsetProxy::interpolateInternal
(
    const tmp<DimensionedField<Type, volMesh>>& tdf
) const
{
    return interpolateInternal(subsetter_, tdf);
}


template<class GeoField>
Foam::tmp<GeoField>
Foam::fvMeshSubsetProxy::interpolate(const GeoField& fld) const
{
    return interpolate(subsetter_, fld);
}


template<class GeoField>
Foam::tmp<GeoField>
Foam::fvMeshSubsetProxy::interpolate(const tmp<GeoField>& tfield) const
{
    return interpolate(subsetter_, tfield);
}


// ************************************************************************* //
