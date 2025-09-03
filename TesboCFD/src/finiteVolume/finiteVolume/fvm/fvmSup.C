/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::zeroField
Foam::fvm::Su
(
    const Foam::zero,
    const GeometricField<Type, fvPatchField, volMesh>& fld
)
{
    return zeroField();
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::Su
(
    const dimensioned<Type>& su,
    const GeometricField<Type, fvPatchField, volMesh>& fld
)
{
    auto tmat = tmp<fvMatrix<Type>>::New
    (
        fld,
        dimVol*su.dimensions()
    );
    auto& mat = tmat.ref();
    const auto& domain = fld.mesh().V();

    if (magSqr(su.value()) > VSMALL)
    {
        mat.source() -= domain*su.value();
    }

    return tmat;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::Su
(
    const DimensionedField<Type, volMesh>& su,
    const GeometricField<Type, fvPatchField, volMesh>& fld
)
{
    auto tmat = tmp<fvMatrix<Type>>::New
    (
        fld,
        dimVol*su.dimensions()
    );
    auto& mat = tmat.ref();
    const auto& domain = fld.mesh().V();

    mat.source() -= domain*su.field();

    return tmat;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::Su
(
    const tmp<DimensionedField<Type, volMesh>>& tsu,
    const GeometricField<Type, fvPatchField, volMesh>& fld
)
{
    tmp<fvMatrix<Type>> tmat = fvm::Su(tsu(), fld);
    tsu.clear();
    return tmat;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::Su
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tsu,
    const GeometricField<Type, fvPatchField, volMesh>& fld
)
{
    tmp<fvMatrix<Type>> tmat = fvm::Su(tsu(), fld);
    tsu.clear();
    return tmat;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::zeroField
Foam::fvm::Sp
(
    const Foam::zero,
    const GeometricField<Type, fvPatchField, volMesh>&
)
{
    return zeroField();
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::Sp
(
    const dimensionedScalar& sp,
    const GeometricField<Type, fvPatchField, volMesh>& fld
)
{
    auto tmat = tmp<fvMatrix<Type>>::New
    (
        fld,
        dimVol*sp.dimensions()*fld.dimensions()
    );
    auto& mat = tmat.ref();
    const auto& domain = fld.mesh().V();

    if (mag(sp.value()) > ROOTVSMALL)
    {
        mat.diag() += domain*sp.value();
    }

    return tmat;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::Sp
(
    const DimensionedField<scalar, volMesh>& sp,
    const GeometricField<Type, fvPatchField, volMesh>& fld
)
{
    auto tmat = tmp<fvMatrix<Type>>::New
    (
        fld,
        dimVol*sp.dimensions()*fld.dimensions()
    );
    auto& mat = tmat.ref();
    const auto& domain = fld.mesh().V();

    mat.diag() += domain*sp.field();

    return tmat;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::Sp
(
    const tmp<DimensionedField<scalar, volMesh>>& tsp,
    const GeometricField<Type, fvPatchField, volMesh>& fld
)
{
    tmp<fvMatrix<Type>> tmat = fvm::Sp(tsp(), fld);
    tsp.clear();
    return tmat;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::Sp
(
    const tmp<volScalarField>& tsp,
    const GeometricField<Type, fvPatchField, volMesh>& fld
)
{
    tmp<fvMatrix<Type>> tmat = fvm::Sp(tsp(), fld);
    tsp.clear();
    return tmat;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::zeroField
Foam::fvm::SuSp
(
    const Foam::zero,
    const GeometricField<Type, fvPatchField, volMesh>& fld
)
{
    return zeroField();
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::SuSp
(
    const dimensionedScalar& susp,
    const GeometricField<Type, fvPatchField, volMesh>& fld
)
{
    auto tmat = tmp<fvMatrix<Type>>::New
    (
        fld,
        dimVol*susp.dimensions()*fld.dimensions()
    );
    auto& mat = tmat.ref();
    const auto& domain = fld.mesh().V();

    if (susp.value() > ROOTVSMALL)
    {
        mat.diag() += domain*susp.value();
    }
    else if (susp.value() < -ROOTVSMALL)
    {
        mat.source() -= domain*susp.value()*fld.primitiveField();
    }

    return tmat;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::SuSp
(
    const DimensionedField<scalar, volMesh>& susp,
    const GeometricField<Type, fvPatchField, volMesh>& fld
)
{
    auto tmat = tmp<fvMatrix<Type>>::New
    (
        fld,
        dimVol*susp.dimensions()*fld.dimensions()
    );
    auto& mat = tmat.ref();
    const auto& domain = fld.mesh().V();

    mat.diag() += domain*max(susp.field(), scalar(0));

    mat.source() -= domain*min(susp.field(), scalar(0))*fld.primitiveField();

    return tmat;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::SuSp
(
    const tmp<DimensionedField<scalar, volMesh>>& tsusp,
    const GeometricField<Type, fvPatchField, volMesh>& fld
)
{
    tmp<fvMatrix<Type>> tmat = fvm::SuSp(tsusp(), fld);
    tsusp.clear();
    return tmat;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::SuSp
(
    const tmp<volScalarField>& tsusp,
    const GeometricField<Type, fvPatchField, volMesh>& fld
)
{
    tmp<fvMatrix<Type>> tmat = fvm::SuSp(tsusp(), fld);
    tsusp.clear();
    return tmat;
}


// ************************************************************************* //
