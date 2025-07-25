/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
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

#include "areaFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::zeroField
Foam::fam::Su
(
    const Foam::zero,
    const GeometricField<Type, faPatchField, areaMesh>& fld
)
{
    return zeroField();
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>>
Foam::fam::Su
(
    const dimensioned<Type>& su,
    const GeometricField<Type, faPatchField, areaMesh>& fld
)
{
    auto tmat = tmp<faMatrix<Type>>::New
    (
        fld,
        dimArea*su.dimensions()
    );
    auto& mat = tmat.ref();
    const auto& domain = fld.mesh().S();

    if (magSqr(su.value()) > VSMALL)
    {
        mat.source() -= domain*su.value();
    }

    return tmat;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>>
Foam::fam::Su
(
    const DimensionedField<Type, areaMesh>& su,
    const GeometricField<Type, faPatchField, areaMesh>& fld
)
{
    auto tmat = tmp<faMatrix<Type>>::New
    (
        fld,
        dimArea*su.dimensions()
    );
    auto& mat = tmat.ref();
    const auto& domain = fld.mesh().S();

    mat.source() -= domain*su.field();

    return tmat;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>>
Foam::fam::Su
(
    const tmp<DimensionedField<Type, areaMesh>>& tsu,
    const GeometricField<Type, faPatchField, areaMesh>& fld
)
{
    tmp<faMatrix<Type>> tmat = fam::Su(tsu(), fld);
    tsu.clear();
    return tmat;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>>
Foam::fam::Su
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tsu,
    const GeometricField<Type, faPatchField, areaMesh>& fld
)
{
    tmp<faMatrix<Type>> tmat = fam::Su(tsu(), fld);
    tsu.clear();
    return tmat;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::zeroField
Foam::fam::Sp
(
    const Foam::zero,
    const GeometricField<Type, faPatchField, areaMesh>& fld
)
{
    return zeroField();
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>>
Foam::fam::Sp
(
    const dimensionedScalar& sp,
    const GeometricField<Type, faPatchField, areaMesh>& fld
)
{
    auto tmat = tmp<faMatrix<Type>>::New
    (
        fld,
        dimArea*sp.dimensions()*fld.dimensions()
    );
    auto& mat = tmat.ref();
    const auto& domain = fld.mesh().S();

    if (mag(sp.value()) > ROOTVSMALL)
    {
        mat.diag() += domain*sp.value();
    }

    return tmat;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>>
Foam::fam::Sp
(
    const DimensionedField<scalar, areaMesh>& sp,
    const GeometricField<Type, faPatchField, areaMesh>& fld
)
{
    auto tmat = tmp<faMatrix<Type>>::New
    (
        fld,
        dimArea*sp.dimensions()*fld.dimensions()
    );
    auto& mat = tmat.ref();
    const auto& domain = fld.mesh().S();

    mat.diag() += domain*sp.field();

    return tmat;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>>
Foam::fam::Sp
(
    const tmp<DimensionedField<scalar, areaMesh>>& tsp,
    const GeometricField<Type, faPatchField, areaMesh>& fld
)
{
    tmp<faMatrix<Type>> tmat = fam::Sp(tsp(), fld);
    tsp.clear();
    return tmat;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>>
Foam::fam::Sp
(
    const tmp<areaScalarField>& tsp,
    const GeometricField<Type, faPatchField, areaMesh>& fld
)
{
    tmp<faMatrix<Type>> tmat = fam::Sp(tsp(), fld);
    tsp.clear();
    return tmat;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::zeroField
Foam::fam::SuSp
(
    const Foam::zero,
    const GeometricField<Type, faPatchField, areaMesh>& fld
)
{
    return zeroField();
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>>
Foam::fam::SuSp
(
    const dimensionedScalar& susp,
    const GeometricField<Type, faPatchField, areaMesh>& fld
)
{
    auto tmat = tmp<faMatrix<Type>>::New
    (
        fld,
        dimArea*susp.dimensions()*fld.dimensions()
    );
    auto& mat = tmat.ref();
    const auto& domain = fld.mesh().S();

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
Foam::tmp<Foam::faMatrix<Type>>
Foam::fam::SuSp
(
    const DimensionedField<scalar, areaMesh>& susp,
    const GeometricField<Type, faPatchField, areaMesh>& fld
)
{
    auto tmat = tmp<faMatrix<Type>>::New
    (
        fld,
        dimArea*susp.dimensions()*fld.dimensions()
    );
    auto& mat = tmat.ref();
    const auto& domain = fld.mesh().S();

    mat.diag() += domain*max(susp.field(), scalar(0));

    mat.source() -= domain*min(susp.field(), scalar(0))*fld.primitiveField();

    return tmat;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>>
Foam::fam::SuSp
(
    const tmp<DimensionedField<scalar, areaMesh>>& tsusp,
    const GeometricField<Type, faPatchField, areaMesh>& fld
)
{
    tmp<faMatrix<Type>> tmat = fam::SuSp(tsusp(), fld);
    tsusp.clear();
    return tmat;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>>
Foam::fam::SuSp
(
    const tmp<areaScalarField>& tsusp,
    const GeometricField<Type, faPatchField, areaMesh>& fld
)
{
    tmp<faMatrix<Type>> tmat = fam::SuSp(tsusp(), fld);
    tsusp.clear();
    return tmat;
}


// ************************************************************************* //
