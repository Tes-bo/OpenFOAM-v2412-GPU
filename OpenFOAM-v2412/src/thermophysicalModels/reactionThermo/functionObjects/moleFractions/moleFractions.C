/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "moleFractions.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ThermoType>
void Foam::moleFractions<ThermoType>::calcMoleFractions()
{
    const auto& thermo =
        mesh_.lookupObject<ThermoType>
        (
            IOobject::groupName(basicThermo::dictName, phaseName_)
        );

    const PtrList<volScalarField>& Y = thermo.composition().Y();

    const volScalarField W(thermo.W());

    forAll(Y, i)
    {
        const dimensionedScalar Wi
        (
            dimMass/dimMoles,
            thermo.composition().W(i)
        );

        X_[i] = W*Y[i]/Wi;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::moleFractions<ThermoType>::moleFractions
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    phaseName_(dict.getOrDefault<word>("phase", word::null))
{
    const word dictName
    (
        IOobject::groupName(basicThermo::dictName, phaseName_)
    );

    if (mesh_.foundObject<ThermoType>(dictName))
    {
        const auto& thermo = mesh_.lookupObject<ThermoType>(dictName);

        const PtrList<volScalarField>& Y = thermo.composition().Y();

        X_.setSize(Y.size());

        forAll(Y, i)
        {
            X_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "X_" + Y[i].name(),
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh_,
                    dimensionedScalar(dimless, Zero)
                )
            );
        }

        calcMoleFractions();
    }
    else
    {
        FatalErrorInFunction
            << "Cannot find thermodynamics model of type "
            << ThermoType::typeName;

        if (!phaseName_.empty())
        {
            FatalError
                << " for phase " << phaseName_;
        }

        FatalError << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
bool Foam::moleFractions<ThermoType>::read
(
    const dictionary& dict
)
{
    if (functionObjects::fvMeshFunctionObject::read(dict))
    {
        phaseName_ = dict.getOrDefault<word>("phase", word::null);

        return true;
    }

    return false;
}



template<class ThermoType>
bool Foam::moleFractions<ThermoType>::execute()
{
    calcMoleFractions();

    return true;
}


// ************************************************************************* //
