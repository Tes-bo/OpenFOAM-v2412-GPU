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

inline const Foam::volScalarField& Foam::incompressibleVars::pInst() const
{
    return pPtr_();
}


inline Foam::volScalarField& Foam::incompressibleVars::pInst()
{
    return pPtr_();
}


inline const Foam::volVectorField& Foam::incompressibleVars::UInst() const
{
    return UPtr_();
}


inline Foam::volVectorField& Foam::incompressibleVars::UInst()
{
    return UPtr_();
}


inline const Foam::surfaceScalarField&
Foam::incompressibleVars::phiInst() const
{
    return phiPtr_();
}


inline Foam::surfaceScalarField& Foam::incompressibleVars::phiInst()
{
    return phiPtr_();
}


inline const Foam::singlePhaseTransportModel&
Foam::incompressibleVars::laminarTransport() const
{
    return laminarTransportPtr_();
}


inline Foam::singlePhaseTransportModel&
Foam::incompressibleVars::laminarTransport()
{
    return laminarTransportPtr_();
}


inline const Foam::autoPtr<Foam::incompressible::turbulenceModel>&
Foam::incompressibleVars::turbulence() const
{
    return turbulence_;
}


inline Foam::autoPtr<Foam::incompressible::turbulenceModel>&
Foam::incompressibleVars::turbulence()
{
    return turbulence_;
}


inline const Foam::autoPtr<Foam::incompressible::RASModelVariables>&
Foam::incompressibleVars::RASModelVariables() const
{
    return RASModelVariables_;
}


inline Foam::autoPtr<Foam::incompressible::RASModelVariables>&
Foam::incompressibleVars::RASModelVariables()
{
    return RASModelVariables_;
}


// ************************************************************************* //
