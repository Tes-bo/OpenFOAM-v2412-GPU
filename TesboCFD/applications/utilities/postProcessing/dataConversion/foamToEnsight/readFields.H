/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2023 OpenCFD Ltd.
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
    Foam

Description
    Helper routines for reading a field or fields,
    for foamToEnsight

SourceFiles
    readFields.C

\*---------------------------------------------------------------------------*/

#ifndef ensight_readFields_H
#define ensight_readFields_H

#include "instantList.H"
#include "IOobjectList.H"
#include "fvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Get the field or FatalError
template<class GeoField>
tmp<GeoField> getField
(
    const IOobject& io,
    const typename GeoField::Mesh& mesh
)
{
    return tmp<GeoField>::New(io, mesh);
}


//- Get the field or return nullptr
template<class GeoField>
tmp<GeoField> getField
(
    const IOobject* io,
    const typename GeoField::Mesh& mesh
)
{
    if (io)
    {
        return tmp<GeoField>::New(*io, mesh);
    }

    return nullptr;
}


//- Get the named field from the objects, or return nullptr.
template<class GeoField>
tmp<GeoField> getField
(
    const typename GeoField::Mesh& mesh,
    const IOobjectList& objects,
    const word& fieldName
)
{
    // Can do something with syncPar on failure ...

    return getField<GeoField>(objects.findObject(fieldName), mesh);
}


//- Convert an internal field to zero-gradient volume field
template<class Type>
tmp<VolumeField<Type>> makeZeroGradientField
(
    const tmp<VolumeInternalField<Type>>& tdf
)
{
    if (tdf)
    {
        auto& df = tdf.ref();

        auto tfield = VolumeField<Type>::New
        (
            df.name(),
            df.mesh(),
            df.dimensions(),
            std::move(df.field()),
            fvPatchFieldBase::zeroGradientType()
        );

        tfield.ref().oriented() = df.oriented();
        tfield.ref().correctBoundaryConditions();

        tdf.clear();

        return tfield;
    }

    tdf.clear();

    return nullptr;
}


//- Convert a volume field to zero-gradient volume field
template<class Type>
tmp<VolumeField<Type>> makeZeroGradientField
(
    const tmp<VolumeField<Type>>& tdf
)
{
    if (tdf)
    {
        auto& df = tdf.ref();

        auto tfield = VolumeField<Type>::New
        (
            df.name(),
            df.mesh(),
            df.dimensions(),
            std::move(df.primitiveFieldRef(false)),  // No update accessTime
            fvPatchFieldBase::zeroGradientType()
        );

        tfield.ref().oriented() = df.oriented();
        tfield.ref().correctBoundaryConditions();

        tdf.clear();

        return tfield;
    }

    tdf.clear();

    return nullptr;
}


//- Check if fields are good to use (available at all times)
//  ignore special fields (_0 fields),
//  ignore fields that are not available for all time-steps
label checkData
(
    const objectRegistry& obr,
    const instantList& timeDirs,
    wordList& objectNames,
    const fileName& local = fileName::null
);


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
