/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
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

#include "volSurfaceMapping.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::volSurfaceMapping::mapToSurface
(
    const GeometricBoundaryField<Type, fvPatchField, volMesh>& bfld,
    Field<Type>& result
) const
{
    // The polyPatch/local-face for each of the faceLabels
    const List<labelPair>& patchFaces = mesh_.whichPatchFaces();

    // FULLDEBUG: checkSize ?
    // or simply result.resize(mesh_.nFaces());

    forAll(patchFaces, i)
    {
        const labelPair& patchAndFace = patchFaces[i];

        if (patchAndFace.first() >= 0)
        {
            result[i] = bfld[patchAndFace];
        }
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::volSurfaceMapping::mapToSurface
(
    const GeometricBoundaryField<Type, fvPatchField, volMesh>& bfld
) const
{
    auto tresult = tmp<Field<Type>>::New(mesh_.nFaces(), Zero);
    mapToSurface(bfld, tresult.ref());
    return tresult;
}


template<class Type>
void Foam::volSurfaceMapping::mapToSurface
(
    const GeometricField<Type, fvPatchField, volMesh>& vfld,
    Field<Type>& result
) const
{
    mapToSurface(vfld.boundaryField(), result);
}


template<class Type>
void Foam::volSurfaceMapping::mapToSurface
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvf,
    Field<Type>& result
) const
{
    mapToSurface(tvf().boundaryField(), result);
    tvf.clear();
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::volSurfaceMapping::mapToSurface
(
    const GeometricField<Type, fvPatchField, volMesh>& vfld
) const
{
    return mapToSurface(vfld.boundaryField());
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::volSurfaceMapping::mapToSurface
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvf
) const
{
    tmp<Field<Type>> tresult(mapToSurface(tvf().boundaryField()));
    tvf.clear();
    return tresult;
}


template<class Type>
void Foam::volSurfaceMapping::mapToSurface
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& fld,
    Field<Type>& result
) const
{
    const auto& bfld = fld.boundaryField();

    // The polyPatch/local-face for each of the faceLabels
    const List<labelPair>& patchFaces = mesh_.whichPatchFaces();

    // FULLDEBUG: checkSize ?
    // or simply result.resize(mesh_.nFaces());

    forAll(patchFaces, i)
    {
        const labelPair& patchAndFace = patchFaces[i];

        if (patchAndFace.first() >= 0)
        {
            // Value from boundary
            result[i] = bfld[patchAndFace];
        }
        else if (patchAndFace.second() >= 0)
        {
            // Value from internal
            result[i] = fld[patchAndFace.second()];
        }
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::volSurfaceMapping::mapToSurface
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& fld
) const
{
    auto tresult = tmp<Field<Type>>::New(mesh_.nFaces(), Zero);
    mapToSurface(fld, tresult.ref());
    return tresult;
}


template<class Type>
void Foam::volSurfaceMapping::mapToSurface
(
    const tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>& tsf,
    Field<Type>& result
) const
{
    mapToSurface(tsf(), result);
    tsf.clear();
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::volSurfaceMapping::mapToSurface
(
    const tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>& tsf
) const
{
    tmp<Field<Type>> tresult(mapToSurface(tsf()));
    tsf.clear();
    return tresult;
}


template<class Type>
void Foam::volSurfaceMapping::mapToSurface
(
    const UPtrList<Field<Type>>& patchFields,
    Field<Type>& result
) const
{
    // The polyPatch/local-face for each of the faceLabels
    const List<labelPair>& patchFaces = mesh_.whichPatchFaces();

    // FULLDEBUG: checkSize ?
    // or simply result.resize(mesh_.nFaces());

    forAll(patchFaces, i)
    {
        const label patchi = patchFaces[i].first();
        const label facei = patchFaces[i].second();

        const auto* pfld = patchFields.get(patchi);

        if (pfld)
        {
            result[i] = (*pfld)[facei];
        }
    }
}


template<class Type>
void Foam::volSurfaceMapping::mapToSurface
(
    const PtrMap<Field<Type>>& patchFields,
    Field<Type>& result
) const
{
    // The polyPatch/local-face for each of the faceLabels
    const List<labelPair>& patchFaces = mesh_.whichPatchFaces();

    // FULLDEBUG: checkSize ?
    // or simply result.resize(mesh_.nFaces());

    forAll(patchFaces, i)
    {
        const label patchi = patchFaces[i].first();
        const label facei = patchFaces[i].second();

        const auto* pfld = patchFields.get(patchi);

        if (pfld)
        {
            result[i] = (*pfld)[facei];
        }
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::volSurfaceMapping::mapToSurface
(
    const UPtrList<Field<Type>>& patchFields
) const
{
    auto tresult = tmp<Field<Type>>::New(mesh_.nFaces(), Zero);
    mapToSurface(patchFields, tresult.ref());
    return tresult;
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::volSurfaceMapping::mapToSurface
(
    const PtrMap<Field<Type>>& patchFields
) const
{
    auto tresult = tmp<Field<Type>>::New(mesh_.nFaces(), Zero);
    mapToSurface(patchFields, tresult.ref());
    return tresult;
}


template<class Type>
void Foam::volSurfaceMapping::mapInternalToSurface
(
    const GeometricBoundaryField<Type, fvPatchField, volMesh>& bfld,
    Field<Type>& result
) const
{
    PtrList<Field<Type>> patchFields;

    // All referenced polyPatches (sorted order)
    const labelList& patches = mesh_.whichPolyPatches();

    if (!patches.empty())
    {
        // maxPolyPatch+1
        patchFields.resize(patches.last()+1);
    }

    // Populate patchInternalField
    for (const label patchi : patches)
    {
        patchFields.set(patchi, bfld[patchi].patchInternalField());
    }

    mapToSurface(patchFields, result);
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::volSurfaceMapping::mapInternalToSurface
(
    const GeometricBoundaryField<Type, fvPatchField, volMesh>& bfld
) const
{
    auto tresult = tmp<Field<Type>>::New(mesh_.nFaces(), Zero);

    mapInternalToSurface(bfld, tresult.ref());

    return tresult;
}


template<class Type>
void Foam::volSurfaceMapping::mapInternalToSurface
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    Field<Type>& result
) const
{
    mapInternalToSurface(vf.boundaryField(), result);
}


template<class Type>
void Foam::volSurfaceMapping::mapInternalToSurface
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvf,
    Field<Type>& result
) const
{
    mapInternalToSurface(tvf().boundaryField(), result);
    tvf.clear();
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::volSurfaceMapping::mapInternalToSurface
(
    const GeometricField<Type, fvPatchField, volMesh>& vfld
) const
{
    return mapInternalToSurface(vfld.boundaryField());
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::volSurfaceMapping::mapInternalToSurface
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvf
) const
{
    tmp<Field<Type>> tresult(mapInternalToSurface(tvf().boundaryField()));
    tvf.clear();
    return tresult;
}


template<class Type>
void Foam::volSurfaceMapping::mapToVolume
(
    const DimensionedField<Type, areaMesh>& af,
    GeometricBoundaryField<Type, fvPatchField, volMesh>& dest,
    const label destPatchi
) const
{
    // The polyPatch/local-face for each of the faceLabels
    const List<labelPair>& patchFaces = mesh_.whichPatchFaces();

    forAll(patchFaces, i)
    {
        const labelPair& patchAndFace = patchFaces[i];

        if
        (
            patchAndFace.first() >= 0
         && (patchAndFace.first() == destPatchi || destPatchi < 0)
        )
        {
            dest[patchAndFace] = af[i];
        }
    }
}


template<class Type>
void Foam::volSurfaceMapping::mapToVolume
(
    const tmp<DimensionedField<Type, areaMesh>>& taf,
    GeometricBoundaryField<Type, fvPatchField, volMesh>& dest,
    const label destPatchi
) const
{
    mapToVolume(taf(), dest, destPatchi);
    taf.clear();
}


template<class Type>
void Foam::volSurfaceMapping::mapToVolume
(
    const DimensionedField<Type, areaMesh>& af,
    GeometricField<Type, fvPatchField, volMesh>& dest,
    const label destPatchi
) const
{
    mapToVolume(af, dest.boundaryFieldRef(), destPatchi);
}


template<class Type>
void Foam::volSurfaceMapping::mapToVolume
(
    const tmp<DimensionedField<Type, areaMesh>>& taf,
    GeometricField<Type, fvPatchField, volMesh>& dest,
    const label destPatchi
) const
{
    mapToVolume(taf, dest.boundaryFieldRef(), destPatchi);
}


template<class Type>
void Foam::volSurfaceMapping::mapToVolume
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& taf,
    GeometricBoundaryField<Type, fvPatchField, volMesh>& dest,
    const label destPatchi
) const
{
    mapToVolume(taf().internalField(), dest, destPatchi);
    taf.clear();
}


template<class Type>
void Foam::volSurfaceMapping::mapToVolume
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& taf,
    GeometricField<Type, fvPatchField, volMesh>& dest,
    const label destPatchi
) const
{
    mapToVolume(taf().internalField(), dest.boundaryFieldRef(), destPatchi);
    taf.clear();
}


template<class Type>
void Foam::volSurfaceMapping::mapToVolumePatch
(
    const DimensionedField<Type, areaMesh>& af,
    Field<Type>& dest,
    const label destPatchi
) const
{
    // The polyPatch/local-face for each of the faceLabels
    const List<labelPair>& patchFaces = mesh_.whichPatchFaces();

    forAll(patchFaces, i)
    {
        const label patchi = patchFaces[i].first();
        const label facei = patchFaces[i].second();

        if (patchi >= 0 && patchi == destPatchi)
        {
            dest[facei] = af[i];
        }
    }
}


template<class Type>
void Foam::volSurfaceMapping::mapToVolumePatch
(
    const tmp<DimensionedField<Type, areaMesh>>& taf,
    Field<Type>& dest,
    const label destPatchi
) const
{
    mapToVolumePatch(taf(), dest, destPatchi);
    taf.clear();
}


template<class Type>
void Foam::volSurfaceMapping::mapToVolumePatch
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& taf,
    Field<Type>& dest,
    const label destPatchi
) const
{
    mapToVolumePatch(taf().internalField(), dest, destPatchi);
    taf.clear();
}


// ************************************************************************* //
