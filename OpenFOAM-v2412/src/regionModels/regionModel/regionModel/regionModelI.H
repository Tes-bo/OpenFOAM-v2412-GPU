/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "regionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const Foam::fvMesh& Foam::regionModels::regionModel::regionMesh() const
{
    const fvMesh* regionPtr = time_.getObjectPtr<fvMesh>(regionName_);

    if (!regionPtr)
    {
        FatalErrorInFunction
            << "Region mesh not available" << abort(FatalError);
    }

    return *regionPtr;
}


inline Foam::fvMesh& Foam::regionModels::regionModel::regionMesh()
{
    fvMesh* regionPtr = time_.getObjectPtr<fvMesh>(regionName_);

    if (!regionPtr)
    {
        FatalErrorInFunction
            << "Region mesh not available" << abort(FatalError);
    }

    return *regionPtr;
}


inline const Foam::dictionary&
Foam::regionModels::regionModel::solution() const
{
    return regionMesh().solutionDict();
}


inline const Foam::IOdictionary&
Foam::regionModels::regionModel::outputProperties() const
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
Foam::regionModels::regionModel::outputProperties()
{
    if (!outputPropertiesPtr_)
    {
        FatalErrorInFunction
            << "outputProperties dictionary not available"
            << abort(FatalError);
    }

    return *outputPropertiesPtr_;
}


inline const Foam::labelList&
Foam::regionModels::regionModel::primaryPatchIDs() const noexcept
{
    return primaryPatchIDs_;
}


inline const Foam::labelList&
Foam::regionModels::regionModel::intCoupledPatchIDs() const noexcept
{
    return intCoupledPatchIDs_;
}


inline bool Foam::regionModels::regionModel::isCoupledPatch
(
    const label regionPatchi
) const
{
    return intCoupledPatchIDs_.found(regionPatchi);
}


inline bool Foam::regionModels::regionModel::isRegionPatch
(
    const label primaryPatchi
) const
{
    return primaryPatchIDs_.found(primaryPatchi);
}


inline Foam::label Foam::regionModels::regionModel::regionPatchID
(
    const label primaryPatchi
) const
{
    const label i = primaryPatchIDs_.find(primaryPatchi);

    return (i >= 0 ? intCoupledPatchIDs_[i] : -1);
}


// ************************************************************************* //
