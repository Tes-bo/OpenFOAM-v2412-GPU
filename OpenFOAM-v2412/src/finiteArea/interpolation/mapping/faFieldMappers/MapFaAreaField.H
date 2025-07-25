/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
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

Class
    Foam::MapFaAreaField

Description
    Map area internal field on topology change.  This is a partial
    template specialisation, see MapGeometricFields.

Author
    Zeljko Tukovic, FMENA
    Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

#ifndef MapFaAreaField_H
#define MapFaAreaField_H

#include "Field.H"
#include "areaFaMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type, class MeshMapper>
class MapInternalField<Type, MeshMapper, areaMesh>
{
public:

    MapInternalField()
    {}

    void operator()
    (
        Field<Type>& field,
        const MeshMapper& mapper
    ) const;
};


template<class Type, class MeshMapper>
void MapInternalField<Type, MeshMapper, areaMesh>::operator()
(
    Field<Type>& field,
    const MeshMapper& mapper
) const
{
    if (field.size() != mapper.areaMap().sizeBeforeMapping())
    {
        FatalErrorInFunction
           << "Incompatible size before mapping.  Field size: " << field.size()
           << " map size: " << mapper.areaMap().sizeBeforeMapping()
           << abort(FatalError);
    }

    field.autoMap(mapper.areaMap());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
