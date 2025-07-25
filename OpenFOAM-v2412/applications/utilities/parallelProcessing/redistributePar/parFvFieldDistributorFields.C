/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

\*---------------------------------------------------------------------------*/

#include "parFvFieldDistributor.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::parFvFieldDistributor::distributeAllFields
(
    const IOobjectList& objects,
    const wordRes& selectedFields
) const
{
    do
    {
        #undef  doLocalCode
        #define doLocalCode(Method)                                           \
        {                                                                     \
            this->Method <scalar> (objects, selectedFields);                  \
            this->Method <vector> (objects, selectedFields);                  \
            this->Method <sphericalTensor> (objects, selectedFields);         \
            this->Method <symmTensor> (objects, selectedFields);              \
            this->Method <tensor> (objects, selectedFields);                  \
        }

        doLocalCode(distributeInternalFields);
        doLocalCode(distributeVolumeFields);
        doLocalCode(distributeSurfaceFields);

        #undef doLocalCode
    }
    while (false);
}


// ************************************************************************* //
