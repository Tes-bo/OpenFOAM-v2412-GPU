/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class PatchType>
template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::vtk::GenericPatchGeoFieldsWriter<PatchType>::getFaceField
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sfld,
    const labelUList& faceAddr
) const
{
    if (this->patch().size() != faceAddr.size())
    {
        FatalErrorInFunction
            << "Inconsistent sizing: patch has "
            << this->patch().size() << " faces, addressing has "
            << faceAddr.size() << " faces!" << nl
            << Foam::exit(FatalError);
    }

    const polyBoundaryMesh& patches = sfld.mesh().boundaryMesh();

    auto tfld = tmp<Field<Type>>::New(faceAddr.size());
    auto iter = tfld.ref().begin();

    for (const label facei : faceAddr)
    {
        const label patchi = patches.whichPatch(facei);

        if (patchi == -1)
        {
            *iter = sfld[facei];
        }
        else
        {
            const label localFacei = facei - patches[patchi].start();
            *iter = sfld.boundaryField()[patchi][localFacei];
        }

        ++iter;
    }

    return tfld;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class PatchType>
template<class Type>
void Foam::vtk::GenericPatchGeoFieldsWriter<PatchType>::write
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& field,
    const labelUList& faceAddr
)
{
    this->GenericPatchWriter<PatchType>::writeCellData
    (
        field.name(),
        getFaceField(field, faceAddr)()
    );
}


template<class PatchType>
template<class Type>
void Foam::vtk::GenericPatchGeoFieldsWriter<PatchType>::write
(
    const GeometricField<Type, faPatchField, areaMesh>& field
)
{
    this->GenericPatchWriter<PatchType>::writeCellData
    (
        field.name(),
        field.primitiveField()
    );
}


// ************************************************************************* //
