/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "LaunderSharmaKE.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASVariables
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(LaunderSharmaKE, 0);
addToRunTimeSelectionTable(RASModelVariables, LaunderSharmaKE, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

LaunderSharmaKE::LaunderSharmaKE
(
    const fvMesh& mesh,
    const solverControl& SolverControl
)
:
    RASModelVariables(mesh, SolverControl)
{
    TMVar1BaseName_ = "k";
    TMVar2BaseName_ = "epsilon";

    TMVar1Ptr_.ref(mesh_.lookupObjectRef<volScalarField>(TMVar1BaseName_));
    TMVar2Ptr_.ref(mesh_.lookupObjectRef<volScalarField>(TMVar2BaseName_));
    nutPtr_.ref(mesh_.lookupObjectRef<volScalarField>(nutBaseName_));

    allocateInitValues();
    allocateMeanFields();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASVariables
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
