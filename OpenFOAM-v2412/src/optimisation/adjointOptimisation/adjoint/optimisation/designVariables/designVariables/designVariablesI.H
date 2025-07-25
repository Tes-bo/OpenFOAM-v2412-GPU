/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 PCOpt/NTUA
    Copyright (C) 2021 FOSS GP
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::labelList&
Foam::designVariables::activeDesignVariables() const
{
    return activeDesignVariables_;
}


bool Foam::designVariables::isMaxInitChangeSet() const
{
    return bool(maxInitChange_);
}


const Foam::autoPtr<Foam::scalar>& Foam::designVariables::getMaxInitChange()
{
    return maxInitChange_;
}


void Foam::designVariables::setMaxInitChange(const scalar maxInitChange)
{
    maxInitChange_.reset(new scalar(maxInitChange));
}


bool Foam::designVariables::resetEta() const
{
    return false;
}


const Foam::autoPtr<Foam::scalarField>&
Foam::designVariables::lowerBounds() const
{
    return lowerBounds_;
}


const Foam::autoPtr<Foam::scalarField>&
Foam::designVariables::upperBounds() const
{
    return upperBounds_;
}


Foam::scalarField& Foam::designVariables::lowerBoundsRef()
{
    return lowerBounds_.ref();
}


Foam::scalarField& Foam::designVariables::upperBoundsRef()
{
    return upperBounds_.ref();
}


// ************************************************************************* //
