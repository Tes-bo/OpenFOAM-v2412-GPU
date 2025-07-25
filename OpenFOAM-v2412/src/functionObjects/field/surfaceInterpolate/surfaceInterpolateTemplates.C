/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2019 OpenCFD Ltd.
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

#include "surfaceInterpolate.H"
#include "volFields.H"
#include "linear.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::functionObjects::surfaceInterpolate::interpolateFields()
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    // Convert field to map
    HashTable<word> fieldMap(2*fieldSet_.size());
    for (const auto& namePair : fieldSet_)
    {
        fieldMap.insert(namePair.first(), namePair.second());
    }


    for (const VolFieldType& fld : obr_.csorted<VolFieldType>())
    {
        if (fieldMap.found(fld.name()))
        {
            // const word sName = "interpolate(" + fld.name() + ')';
            word& sName = fieldMap[fld.name()];

            if (obr_.found(sName))
            {
                Log << "        updating field " << sName << endl;
            }
            else
            {
                Log << "        interpolating "
                    << fld.name() << " to create " << sName << endl;
            }

            store(sName, linearInterpolate(fld));
        }
    }
}


// ************************************************************************* //
