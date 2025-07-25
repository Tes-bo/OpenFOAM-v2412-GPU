/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2021 OpenCFD Ltd.
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
bool Foam::functionObjects::DMD::getSnapshot()
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceFieldType;

    if (foundObject<VolFieldType>(fieldName_))
    {
        return storeSnapshot<VolFieldType>();
    }
    else if (foundObject<SurfaceFieldType>(fieldName_))
    {
        return storeSnapshot<SurfaceFieldType>();
    }

    return false;
}


template<class GeoFieldType>
bool Foam::functionObjects::DMD::storeSnapshot()
{
    if (step_ == 0)
    {
        initialise();
    }

    if (z_.size() == 1)
    {
        return true;
    }

    // Move previous-time snapshot into previous-time slot in "z"
    // Effectively moves the lower half of "z" to its upper half
    std::rotate(z_.begin(), z_.begin() + nSnap_, z_.end());

    // Copy new current-time snapshot into current-time slot in "z"
    // Effectively copies the new field elements into the lower half of "z"
    const label nComps =
        pTraits<typename GeoFieldType::value_type>::nComponents;

    const GeoFieldType& field = lookupObject<GeoFieldType>(fieldName_);

    label rowi = nSnap_;
    if (patches_.empty())
    {
        const label nField = field.size();

        for (direction dir = 0; dir < nComps; ++dir)
        {
            z_.subColumn(0, rowi, nField) = field.component(dir);
            rowi += nField;
        }
    }
    else
    {
        const labelList patchis
        (
            mesh_.boundaryMesh().patchSet(patches_).sortedToc()
        );

        for (const label patchi : patchis)
        {
            const typename GeoFieldType::Boundary& bf = field.boundaryField();

            const Field<typename GeoFieldType::value_type>& pbf = bf[patchi];

            const label nField = pbf.size();

            if (nField > 0)
            {
                for (direction dir = 0; dir < nComps; ++dir)
                {
                    z_.subColumn(0, rowi, nField) = pbf.component(dir);
                    rowi += nField;
                }
            }
        }
    }

    return true;
}


template<class Type>
bool Foam::functionObjects::DMD::nComponents
(
    const word& fieldName,
    label& nComps
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceFieldType;

    if (mesh_.foundObject<VolFieldType>(fieldName))
    {
        nComps = pTraits<typename VolFieldType::value_type>::nComponents;
        return true;
    }
    else if (mesh_.foundObject<SurfaceFieldType>(fieldName))
    {
        nComps = pTraits<typename SurfaceFieldType::value_type>::nComponents;
        return true;
    }

    return false;
}


// ************************************************************************* //
