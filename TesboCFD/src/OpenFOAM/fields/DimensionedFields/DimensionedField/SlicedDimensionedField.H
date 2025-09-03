/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenFOAM Foundation
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

Class
    Foam::SlicedDimensionedField

Description
    Specialisation of DimensionedField that holds a slice of a given
    field so that it acts as a DimensionedField.

\*---------------------------------------------------------------------------*/

#ifndef Foam_SlicedDimensionedField_H
#define Foam_SlicedDimensionedField_H

#include "DimensionedField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class SlicedDimensionedField Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class GeoMesh>
class SlicedDimensionedField
:
    public DimensionedField<Type, GeoMesh>
{
public:

    // Public Typedefs

        //- Type of mesh on which this SlicedDimensionedField is instantiated
        typedef typename GeoMesh::Mesh Mesh;


    // Constructors

        //- Construct from components and field to slice
        SlicedDimensionedField
        (
            const IOobject& io,
            const Mesh& mesh,
            const dimensionSet& dims,
            const Field<Type>& iField
        )
        :
            DimensionedField<Type, GeoMesh>(io, mesh, dims, Field<Type>())
        {
            // Set internalField to the slice of the complete field
            UList<Type>::shallowCopy
            (
                SubList<Type>(iField, GeoMesh::size(mesh))
            );
        }


    //- Destructor
    ~SlicedDimensionedField()
    {
        // Set internalField to nullptr to avoid deletion of underlying field
        UList<Type>::shallowCopy(nullptr);
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
