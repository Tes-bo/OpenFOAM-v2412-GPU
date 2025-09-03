/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2023 OpenCFD Ltd.
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

#include "ensightOutputVolField.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::label Foam::functionObjects::ensightWrite::writeVolFieldsImpl
(
    ensightOutput::floatBufferType& scratch,
    const fvMeshSubset& proxy,
    const wordHashSet& candidateNames
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> GeoField;

    const fvMesh& baseMesh = proxy.baseMesh();

    label count = 0;

    for
    (
        const GeoField& origField
      : baseMesh.csorted<GeoField>(candidateNames)
    )
    {
        const word& fieldName = origField.name();

        auto tfield = fvMeshSubsetProxy::interpolate(proxy, origField);
        const auto& field = tfield();

        autoPtr<ensightFile> os = ensCase().newData<Type>(fieldName);

        ensightOutput::writeVolField<Type>
        (
            scratch,
            os.ref(),
            field,
            ensMesh(),
            caseOpts_.nodeValues()
        );

        Log << "    " << fieldName << "\n";

        ++count;
    }

    return count;
}


// ************************************************************************* //
