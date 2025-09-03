/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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

#include "Time.H"
#include "sampledSurface.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

inline bool Foam::functionObjects::fieldValues::surfaceFieldValue::
withSurfaceFields() const
{
    return
    (
        stFaceZone == regionType_
     || stPatch == regionType_
     || (sampledPtr_ && sampledPtr_->withSurfaceFields())
    );
}


inline bool Foam::functionObjects::fieldValues::surfaceFieldValue::
withTopologicalMerge() const noexcept
{
    return (stFaceZone == regionType_ || stPatch == regionType_);
}


inline const Foam::labelList&
Foam::functionObjects::fieldValues::surfaceFieldValue::faceId()
const noexcept
{
    return faceId_;
}


inline const Foam::labelList&
Foam::functionObjects::fieldValues::surfaceFieldValue::facePatch()
const noexcept
{
    return facePatchId_;
}


inline const Foam::boolList&
Foam::functionObjects::fieldValues::surfaceFieldValue::faceFlip()
const noexcept
{
    return faceFlip_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::functionObjects::fieldValues::surfaceFieldValue::
is_magOp()
const noexcept
{
    // Operation specifically tagged to use mag
    return (operation_ & typeAbsolute);
}


inline bool Foam::functionObjects::fieldValues::surfaceFieldValue::
is_weightedOp() const noexcept
{
    // Operation specifically tagged to require a weight field
    return (operation_ & typeWeighted);
}


inline Foam::functionObjects::fieldValues::surfaceFieldValue::regionTypes
Foam::functionObjects::fieldValues::surfaceFieldValue::regionType()
const noexcept
{
    return regionType_;
}


inline Foam::fileName
Foam::functionObjects::fieldValues::surfaceFieldValue::outputDir() const
{
    return baseFileDir()/name()/"surface"/time_.timeName();
}


// ************************************************************************* //
