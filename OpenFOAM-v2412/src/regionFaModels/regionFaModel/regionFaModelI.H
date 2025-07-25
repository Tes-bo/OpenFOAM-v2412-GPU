/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

#include "regionFaModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const Foam::faMesh& Foam::regionModels::regionFaModel::regionMesh() const
{
    if (!regionMeshPtr_)
    {
        FatalErrorInFunction
            << "Region mesh not available"
            << abort(FatalError);
    }

    return *regionMeshPtr_;
}


inline Foam::faMesh& Foam::regionModels::regionFaModel::regionMesh()
{
    if (!regionMeshPtr_)
    {
        FatalErrorInFunction
            << "Region mesh not available"
            << abort(FatalError);
    }

    return *regionMeshPtr_;
}


inline const Foam::IOdictionary&
Foam::regionModels::regionFaModel::outputProperties() const
{
    if (!outputPropertiesPtr_)
    {
        FatalErrorInFunction
            << "outputProperties dictionary not available"
            << abort(FatalError);
    }
    return *outputPropertiesPtr_;
}


inline Foam::IOdictionary&
Foam::regionModels::regionFaModel::outputProperties()
{
    if (!outputPropertiesPtr_)
    {
        FatalErrorInFunction
            << "outputProperties dictionary not available"
            << abort(FatalError);
    }

    return *outputPropertiesPtr_;
}


inline const Foam::dictionary&
Foam::regionModels::regionFaModel::solution() const
{
    return regionMesh().solutionDict();
}


inline const Foam::labelList&
Foam::regionModels::regionFaModel::primaryPatchIDs() const
{
    return regionMesh().whichPolyPatches();
}


inline bool Foam::regionModels::regionFaModel::isRegionPatch
(
    const label patchi
) const
{
    return primaryPatchIDs().found(patchi);
}


// ************************************************************************* //
