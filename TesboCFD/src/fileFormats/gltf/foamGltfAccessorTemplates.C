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

#include "foamGltfAccessor.H"
#include "exprTraits.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Type>
Foam::string Foam::glTF::accessor::getValueType()
{
    switch (exprTypeTraits<Type>::value)
    {
        case exprTypeTraits<label>::value :
        case exprTypeTraits<scalar>::value :
        {
            return "SCALAR";
        }

        case exprTypeTraits<vector>::value :
        {
            return "VEC3";
        }

        case exprTypeTraits<tensor>::value :
        {
            return "MAT3";
        }

        default:
        {
            FatalErrorInFunction
                << "Unable to process "
                << pTraits<Type>::typeName << " fields"
                << abort(FatalError);
        }
    }

    return string();
}


template<class Type>
Foam::string Foam::glTF::accessor::toString(const Type& val)
{
    OStringStream buf;
    buf << "[ ";
    for (direction dir = 0; dir < pTraits<Type>::nComponents; ++dir)
    {
        if (dir) buf << ", ";
        buf << float(component(val, dir));
    }
    buf << " ]";

    return buf.str();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::glTF::accessor::set(const Field<Type>& fld, bool calcMinMax)
{
    count_ = fld.size();

    type_ = accessor::getValueType<Type>();

    componentType_ = key(componentTypes::FLOAT);

    minMax_ = calcMinMax;

    if (minMax_)
    {
        Type minValue = min(fld);
        Type maxValue = max(fld);

        min_ = accessor::toString(minValue);
        max_ = accessor::toString(maxValue);
    }
}


// ************************************************************************* //
