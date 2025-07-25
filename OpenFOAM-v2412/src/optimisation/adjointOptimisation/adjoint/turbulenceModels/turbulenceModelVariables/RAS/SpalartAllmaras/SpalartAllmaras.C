/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
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

#include "SpalartAllmaras.H"
#include "wallDist.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASVariables
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(SpalartAllmaras, 0);
addToRunTimeSelectionTable(RASModelVariables, SpalartAllmaras, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

SpalartAllmaras::SpalartAllmaras
(
    const fvMesh& mesh,
    const solverControl& SolverControl
)
:
    RASModelVariables(mesh, SolverControl)
{
    TMVar1BaseName_ = "nuTilda";

    TMVar1Ptr_.ref(mesh_.lookupObjectRef<volScalarField>(TMVar1BaseName_));
    nutPtr_.ref(mesh_.lookupObjectRef<volScalarField>(nutBaseName_));

    // The wall dist name can vary depending on how wallDist was
    // constructed. Grab the field directly from wallDist

    distPtr_.ref
    (
        const_cast<volScalarField&>(wallDist::New(mesh_).y())
    );

    allocateInitValues();
    allocateMeanFields();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volScalarField> SpalartAllmaras::nutJacobianVar1
(
    const singlePhaseTransportModel& laminarTransport
) const
{
    auto tnutJacobian = volScalarField::New
    (
        "nutJacobianVar1",
        IOobject::NO_REGISTER,
        mesh_,
        dimensionedScalar(dimless, Zero)
    );
    auto& nutJacobian = tnutJacobian.ref();

    const tmp<volScalarField> tnu(laminarTransport.nu());
    const volScalarField& nu = tnu();
    const volScalarField& nuTilda = TMVar1();

    volScalarField chi(nuTilda/nu);
    volScalarField chi3(pow3(chi));

    const scalar Cv13 = pow3(7.1);
    volScalarField fv1(chi3/(chi3 + Cv13));
    volScalarField fv1ByChi2Sqr(sqr(chi/(chi3 + Cv13)));
    volScalarField Cdfv1(3.0*Cv13*fv1ByChi2Sqr);

    nutJacobian = Cdfv1*chi + fv1;

    return tnutJacobian;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASVariables
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
