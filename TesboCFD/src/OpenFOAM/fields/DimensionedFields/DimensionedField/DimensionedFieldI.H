/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
inline const typename GeoMesh::Mesh&
Foam::DimensionedField<Type, GeoMesh>::mesh() const noexcept
{
    return mesh_;
}


template<class Type, class GeoMesh>
inline const Foam::dimensionSet&
Foam::DimensionedField<Type, GeoMesh>::dimensions() const noexcept
{
    return dimensions_;
}


template<class Type, class GeoMesh>
inline Foam::dimensionSet&
Foam::DimensionedField<Type, GeoMesh>::dimensions() noexcept
{
    return dimensions_;
}


template<class Type, class GeoMesh>
inline Foam::orientedType
Foam::DimensionedField<Type, GeoMesh>::oriented() const noexcept
{
    return oriented_;
}


template<class Type, class GeoMesh>
inline Foam::orientedType&
Foam::DimensionedField<Type, GeoMesh>::oriented() noexcept
{
    return oriented_;
}


template<class Type, class GeoMesh>
inline bool Foam::DimensionedField<Type, GeoMesh>::is_oriented() const noexcept
{
    return oriented_.is_oriented();
}


template<class Type, class GeoMesh>
inline void Foam::DimensionedField<Type, GeoMesh>::setOriented(bool on) noexcept
{
    oriented_.setOriented(on);
}


template<class Type, class GeoMesh>
inline const Foam::Field<Type>&
Foam::DimensionedField<Type, GeoMesh>::field() const noexcept
{
    return *this;
}


template<class Type, class GeoMesh>
inline Foam::Field<Type>&
Foam::DimensionedField<Type, GeoMesh>::field() noexcept
{
    return *this;
}


// ************************************************************************* //
