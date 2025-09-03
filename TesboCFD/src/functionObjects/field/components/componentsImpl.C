/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class GeoFieldType>
bool Foam::functionObjects::components::calcComponents
(
    const GeoFieldType& field
)
{
    typedef typename GeoFieldType::value_type Type;

    resultNames_.resize(pTraits<Type>::nComponents);

    bool stored = true;

    for (direction i = 0; i < pTraits<Type>::nComponents; ++i)
    {
        resultName_ = fieldName_ + word(pTraits<Type>::componentNames[i]);
        resultNames_[i] = resultName_;

        stored = stored && store(resultName_, field.component(i));
    }

    return stored;
}


template<class Type>
bool Foam::functionObjects::components::calcComponents()
{
    const auto* vfield = cfindObject<VolumeField<Type>>(fieldName_);
    if (vfield)
    {
        return calcComponents(*vfield);
    }

    const auto* sfield = cfindObject<SurfaceField<Type>>(fieldName_);
    if (sfield)
    {
        return calcComponents(*sfield);
    }

    return false;
}


// ************************************************************************* //
