/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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
    Foam::MapInternalField

Description
    Generic internal field mapper.  For "real" mapping, add template
    specialisations for mapping of internal fields depending on mesh
    type.

\*---------------------------------------------------------------------------*/

#ifndef MapGeometricFields_H
#define MapGeometricFields_H

#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type, class MeshMapper, class GeoMesh>
class MapInternalField
{
public:

    MapInternalField() = default;

    void operator()
    (
        DimensionedField<Type, GeoMesh>& field,
        const MeshMapper& mapper
    ) const;
};


//- Generic Geometric field mapper.
//  For "real" mapping, add template specialisations
//  for mapping of internal fields depending on mesh type.
template
<
    class Type,
    template<class> class PatchField,
    class MeshMapper,
    class GeoMesh
>
void MapGeometricFields
(
    const MeshMapper& mapper
)
{
    typedef GeometricField<Type, PatchField, GeoMesh> FieldType;

    const UPtrList<const FieldType> fields
    (
        mapper.thisDb().objectRegistry::template csorted<FieldType>()
    );

    // It is necessary to enforce that all old-time fields are stored
    // before the mapping is performed.  Otherwise, if the
    // old-time-level field is mapped before the field itself, sizes
    // will not match.

    for (const auto& field : fields)
    {
        //Note: check can be removed once pointFields are actually stored on
        //      the pointMesh instead of now on the polyMesh!
        if (&field.mesh() == &mapper.mesh())
        {
            const_cast<FieldType&>(field).storeOldTimes();
        }
    }

    for (const auto& field : fields)
    {
        FieldType& fld = const_cast<FieldType&>(field);

        if (&field.mesh() == &mapper.mesh())
        {
            if (polyMesh::debug)
            {
                Info<< "Mapping "
                    << FieldType::typeName << ' ' << field.name() << endl;
            }

            // Map the internal field
            MapInternalField<Type, MeshMapper, GeoMesh>()
            (
                fld.internalFieldRef(),
                mapper
            );

            // Map the patch fields
            auto& bfield = fld.boundaryFieldRef();

            forAll(bfield, patchi)
            {
                // Cannot check sizes for patch fields because of
                // empty fields in FV and because point fields get their size
                // from the patch which has already been resized
                //

                bfield[patchi].autoMap(mapper.boundaryMap()[patchi]);
            }

            fld.instance() = fld.time().timeName();
        }
        else if (polyMesh::debug)
        {
            Info<< "Not mapping "
                << FieldType::typeName << ' ' << field.name()
                << " since originating mesh differs from that of mapper."
                << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
