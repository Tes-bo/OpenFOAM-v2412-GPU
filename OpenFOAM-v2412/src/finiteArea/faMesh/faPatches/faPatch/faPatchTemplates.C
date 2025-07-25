/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

\*----------------------------------------------------------------------------*/

#include "faPatch.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::faPatch::patchInternalField
(
    const UList<Type>& internalData,
    const labelUList& addressing,
    Field<Type>& pfld
) const
{
    const label len = this->size();

    pfld.resize_nocopy(len);

    for (label i = 0; i < len; ++i)
    {
        pfld[i] = internalData[addressing[i]];
    }
}


template<class Type>
void Foam::faPatch::patchInternalField
(
    const UList<Type>& internalData,
    Field<Type>& pfld
) const
{
    patchInternalField(internalData, this->edgeFaces(), pfld);
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::faPatch::patchInternalField
(
    const UList<Type>& internalData
) const
{
    auto tpfld = tmp<Field<Type>>::New();
    patchInternalField(internalData, this->edgeFaces(), tpfld.ref());
    return tpfld;
}


template<class GeometricField, class AnyType>
const typename GeometricField::Patch& Foam::faPatch::patchField
(
    const GeometricField& gf
) const
{
    return gf.boundaryField()[this->index()];
}


// ************************************************************************* //
