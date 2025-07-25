/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
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

#include "fvMatrices.H"
#include "volFields.H"
#include "DimensionedField.H"

Foam::tmp<Foam::DimensionedField<Foam::vector, Foam::volMesh>>
Foam::coalCloudList::UTrans() const
{
    auto tfld = volVectorField::Internal::New
    (
        "UTransEff",
        IOobject::NO_REGISTER,
        mesh_,
        dimensionedVector(dimMass*dimVelocity, Zero)
    );
    auto& fld = tfld.ref();

    forAll(*this, i)
    {
        fld += operator[](i).UTrans();
    }

    return tfld;
}


Foam::tmp<Foam::fvVectorMatrix> Foam::coalCloudList::SU
(
    volVectorField& U
) const
{
    auto tfvm = tmp<fvVectorMatrix>::New(U, dimForce);
    auto& fvm = tfvm.ref();

    forAll(*this, i)
    {
        fvm += operator[](i).SU(U);
    }

    return tfvm;
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::coalCloudList::hsTrans() const
{
    auto tfld = volScalarField::Internal::New
    (
        "hsTransEff",
        IOobject::NO_REGISTER,
        mesh_,
        dimensionedScalar(dimEnergy, Zero)
    );
    auto& fld = tfld.ref();

    forAll(*this, i)
    {
        fld += operator[](i).hsTrans();
    }

    return tfld;
}


Foam::tmp<Foam::fvScalarMatrix> Foam::coalCloudList::Sh
(
    volScalarField& hs
) const
{
    auto tfvm = tmp<fvScalarMatrix>::New(hs, dimEnergy/dimTime);
    auto& fvm = tfvm.ref();

    forAll(*this, i)
    {
        fvm += operator[](i).Sh(hs);
    }

    return tfvm;
}


Foam::tmp<Foam::fvScalarMatrix> Foam::coalCloudList::SYi
(
    const label ii,
    volScalarField& Yi
) const
{
    auto tfvm = tmp<fvScalarMatrix>::New(Yi, dimMass/dimTime);
    auto& fvm = tfvm.ref();

    forAll(*this, i)
    {
        fvm += operator[](i).SYi(ii, Yi);
    }

    return tfvm;
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::coalCloudList::rhoTrans() const
{
    auto tfld = volScalarField::Internal::New
    (
        "rhoTransEff",
        IOobject::NO_REGISTER,
        mesh_,
        dimensionedScalar(dimMass, Zero)
    );
    auto& fld = tfld.ref();

    forAll(*this, i)
    {
        forAll(operator[](i).rhoTrans(), j)
        {
            fld += operator[](i).rhoTrans()[j];
        }
    }

    return tfld;
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::coalCloudList::Srho() const
{
    auto tfld = volScalarField::Internal::New
    (
        "rhoTransEff",
        IOobject::NO_REGISTER,
        mesh_,
        dimensionedScalar(dimDensity/dimTime, Zero)
    );
    auto& fld = tfld.ref();

    forAll(*this, i)
    {
        fld += operator[](i).Srho();
    }

    return tfld;
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::coalCloudList::Srho
(
    const label i
) const
{
    auto tfld = volScalarField::Internal::New
    (
        "rhoTransEff",
        IOobject::NO_REGISTER,
        mesh_,
        dimensionedScalar(dimDensity/dimTime, Zero)
    );
    auto& fld = tfld.ref();

    forAll(*this, j)
    {
        fld += operator[](j).Srho(i);
    }

    return tfld;
}


Foam::tmp<Foam::fvScalarMatrix> Foam::coalCloudList::Srho
(
    volScalarField& rho
) const
{
    auto tfvm = tmp<fvScalarMatrix>::New(rho, dimMass/dimTime);
    auto& fvm = tfvm.ref();

    forAll(*this, i)
    {
        fvm += operator[](i).Srho(rho);
    }

    return tfvm;
}


// ************************************************************************* //
