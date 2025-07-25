/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

Typedef
    Foam::incompressible::turbulenceModel

Typedef
    Foam::incompressible::RASModel

Typedef
    Foam::incompressible::LESModel

Description
    Typedefs for turbulence, RAS and LES models for incompressible flow
    based on the standard laminar transport package. Density is explicity
    referenced to consider the isothermal variable-density effect, as
    described in:
    \verbatim
        Fan, W. & Anglart, H. (2020).
        varRhoTurbVOF: A new set of volume of fluid solvers for turbulent
        isothermal multiphase flows in OpenFOAM.
        Computer Physics Communications, 247, 106876
    \endverbatim
VoFphaseTurbulentTransportModel.H
SourceFiles
    VoFphaseTurbulentTransportModel.C
    VoFphaseTurbulentTransportModels.C


\*---------------------------------------------------------------------------*/

#ifndef VoFphaseTurbulentTransportModel_H
#define VoFphaseTurbulentTransportModel_H

#include "PhaseIncompressibleTurbulenceModel.H"
#include "laminarModel.H"
#include "RASModel.H"
#include "LESModel.H"
#include "incompressible/transportModel/transportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace incompressible
    {
        typedef PhaseIncompressibleTurbulenceModel<transportModel>
            phaseIncompressibleTurbulenceModel;

        template<class BasicCompressibleTurbulenceModel>
        autoPtr<BasicCompressibleTurbulenceModel> New
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const typename BasicCompressibleTurbulenceModel::transportModel&
                transport,
            const word& propertiesName = turbulenceModel::propertiesName
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "VoFphaseTurbulentTransportModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
