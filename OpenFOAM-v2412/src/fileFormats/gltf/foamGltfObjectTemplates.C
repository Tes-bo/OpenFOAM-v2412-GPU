/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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

template<class Type>
void Foam::glTF::object::addData(const Type& fld)
{
    const direction nCmpts = pTraits<typename Type::value_type>::nComponents;

    label count = data_.size();
    data_.resize(data_.size() + fld.size()*nCmpts);

    forAll(fld, fieldi)
    {
        for (direction d = 0; d < nCmpts; ++d)
        {
            data_[count++] = component(fld[fieldi], d);
        }
    }
}


template<class Type1, class Type2>
void Foam::glTF::object::addData(const Type1& fld1, const Type2& fld2)
{
    if (fld1.size() != fld2.size())
    {
        FatalErrorInFunction
            << "Field lengths must be the same. Field1:"
            << fld1.size() << " Field2:" << fld2.size()
            << abort(FatalError);
    }

    const direction nCmpts1 = pTraits<typename Type1::value_type>::nComponents;
    const direction nCmpts2 = pTraits<typename Type2::value_type>::nComponents;

    label count = data_.size();
    data_.resize(data_.size() + fld1.size()*(nCmpts1 + nCmpts2));

    forAll(fld1, fieldi)
    {
        for (direction d = 0; d < nCmpts1; ++d)
        {
            data_[count++] = component(fld1[fieldi], d);
        }

        for (direction d = 0; d < nCmpts2; ++d)
        {
            data_[count++] = component(fld2[fieldi], d);
        }
    }
}


// ************************************************************************* //
