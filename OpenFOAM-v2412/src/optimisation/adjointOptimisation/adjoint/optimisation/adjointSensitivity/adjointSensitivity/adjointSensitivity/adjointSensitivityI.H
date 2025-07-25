/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 PCOpt/NTUA
    Copyright (C) 2023 FOSS GP
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
    along with OpenFOAM.  If not, see <Foam::http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::autoPtr<Foam::volTensorField>&
Foam::adjointSensitivity::gradDxDbMult() const
{
    return gradDxDbMult_;
}


inline Foam::autoPtr<Foam::volTensorField>&
Foam::adjointSensitivity::gradDxDbMult()
{
    return gradDxDbMult_;
}

inline const Foam::autoPtr<Foam::scalarField>&
Foam::adjointSensitivity::divDxDbMult() const
{
    return divDxDbMult_;
}


inline const Foam::autoPtr<Foam::boundaryVectorField>&
Foam::adjointSensitivity::dxdbMult() const
{
    return dxdbMult_;
}


inline const Foam::autoPtr<Foam::boundaryVectorField>&
Foam::adjointSensitivity::dSfdbMult() const
{
    return dSfdbMult_;
}

inline const Foam::autoPtr<Foam::boundaryVectorField>&
Foam::adjointSensitivity::dnfdbMult() const
{
    return dnfdbMult_;
}

inline const Foam::autoPtr<Foam::boundaryVectorField>&
Foam::adjointSensitivity::dxdbDirectMult() const
{
    return dxdbDirectMult_;
}

inline const Foam::autoPtr<Foam::pointBoundaryVectorField>&
Foam::adjointSensitivity::pointDxDbDirectMult() const
{
    return pointDxDbDirectMult_;
}

inline const Foam::autoPtr<Foam::boundaryVectorField>&
Foam::adjointSensitivity::bcDxDbMult() const
{
    return bcDxDbMult_;
}

inline const Foam::autoPtr<Foam::vectorField>&
Foam::adjointSensitivity::optionsDxDbMult() const
{
    return optionsDxDbMult_;
}



// ************************************************************************* //
