/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
    Copyright (C) 2019 OpenCFD Ltd.
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

inline const Foam::word& Foam::NURBS3DVolume::name() const
{
    return name_;
}


inline const Foam::boolList& Foam::NURBS3DVolume::getActiveCPs() const
{
    return activeControlPoints_;
}


inline const Foam::boolList&
Foam::NURBS3DVolume::getActiveDesignVariables() const
{
    return activeDesignVariables_;
}


inline const Foam::vectorField& Foam::NURBS3DVolume::getControlPoints() const
{
    return cps_;
}


inline Foam::vectorField& Foam::NURBS3DVolume::getControlPoints()
{
    return cps_;
}


inline bool Foam::NURBS3DVolume::confineUMovement() const
{
    return confineUMovement_;
}


inline bool Foam::NURBS3DVolume::confineVMovement() const
{
    return confineVMovement_;
}


inline bool Foam::NURBS3DVolume::confineWMovement() const
{
    return confineWMovement_;
}


const Foam::NURBSbasis& Foam::NURBS3DVolume::basisU() const
{
    return basisU_;
}


const Foam::NURBSbasis& Foam::NURBS3DVolume::basisV() const
{
    return basisV_;
}


const Foam::NURBSbasis& Foam::NURBS3DVolume::basisW() const
{
    return basisW_;
}


inline Foam::Vector<Foam::label>
Foam::NURBS3DVolume::nCPsPerDirection() const
{
    return Vector<label>(basisU_.nCPs(), basisV_.nCPs(), basisW_.nCPs());
}


const Foam::fvMesh& Foam::NURBS3DVolume::mesh() const
{
    return mesh_;
}


const Foam::dictionary& Foam::NURBS3DVolume::dict() const
{
    return dict_;
}


// ************************************************************************* //
