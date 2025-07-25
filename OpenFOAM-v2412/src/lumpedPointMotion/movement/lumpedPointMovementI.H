/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2020 OpenCFD Ltd.
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


inline bool Foam::lumpedPointMovement::empty() const
{
    return state0_.empty();
}


inline Foam::label Foam::lumpedPointMovement::size() const
{
    return state0_.size();
}


inline Foam::label Foam::lumpedPointMovement::ownerId() const
{
    return ownerId_;
}


inline void Foam::lumpedPointMovement::ownerId(label id)
{
    ownerId_ = id;
}


inline bool
Foam::lumpedPointMovement::hasPatchControl(const label patchIndex) const
{
    return patchControls_.found(patchIndex);
}


inline bool
Foam::lumpedPointMovement::hasInterpolator(const label patchIndex) const
{
    const auto iter = patchControls_.cfind(patchIndex);

    return (iter.good() && iter().interp_.size());
}


inline const Foam::externalFileCoupler&
Foam::lumpedPointMovement::coupler() const
{
    return coupler_;
}


inline Foam::externalFileCoupler& Foam::lumpedPointMovement::coupler()
{
    return coupler_;
}


inline const Foam::lumpedPointState& Foam::lumpedPointMovement::state0() const
{
    return state0_;
}


inline const Foam::lumpedPointState& Foam::lumpedPointMovement::state() const
{
    return state_;
}


inline const Foam::point& Foam::lumpedPointMovement::origin() const
{
    return origin_;
}


inline void Foam::lumpedPointMovement::scalePoints
(
    lumpedPointState& state
) const
{
    state.scalePoints(scaleInput_[scalingType::LENGTH]);
}


inline Foam::scalar Foam::lumpedPointMovement::relax() const
{
    return relax_;
}


inline Foam::scalar& Foam::lumpedPointMovement::relax()
{
    return relax_;
}


inline const Foam::word& Foam::lumpedPointMovement::inputName() const
{
    return inputName_;
}


inline const Foam::word& Foam::lumpedPointMovement::outputName() const
{
    return outputName_;
}


inline const Foam::word& Foam::lumpedPointMovement::logName() const
{
    return logName_;
}


inline Foam::lumpedPointState::inputFormatType
Foam::lumpedPointMovement::inputFormat() const
{
    return inputFormat_;
}


Foam::lumpedPointMovement::outputFormatType
Foam::lumpedPointMovement::outputFormat() const
{
    return outputFormat_;
}


inline Foam::quaternion::eulerOrder
Foam::lumpedPointMovement::rotationOrder() const
{
    return state0().rotationOrder();
}


inline bool Foam::lumpedPointMovement::degrees() const
{
    return state0().degrees();
}


inline bool Foam::lumpedPointMovement::hasMapping() const
{
    return !patchControls_.empty();
}


// ************************************************************************* //
