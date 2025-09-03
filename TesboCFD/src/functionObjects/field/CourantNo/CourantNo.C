/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2020-2024 OpenCFD Ltd.
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

#include "CourantNo.H"
#include "surfaceFields.H"
#include "fvcSurfaceIntegrate.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(CourantNo, 0);
    addToRunTimeSelectionTable(functionObject, CourantNo, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField::Internal>
Foam::functionObjects::CourantNo::byRho
(
    const tmp<volScalarField::Internal>& Co
) const
{
    if (Co().dimensions() == dimDensity)
    {
        return Co/obr_.lookupObject<volScalarField>(rhoName_);
    }

    return Co;
}


bool Foam::functionObjects::CourantNo::calc()
{
    if (foundObject<surfaceScalarField>(fieldName_))
    {
        const surfaceScalarField& phi =
            lookupObject<surfaceScalarField>(fieldName_);

        tmp<volScalarField::Internal> Coi
        (
            byRho
            (
                (0.5*mesh_.time().deltaT())
               *fvc::surfaceSum(mag(phi))()()
               /mesh_.V()
            )
        );

        auto* resultPtr = getObjectPtr<volScalarField>(resultName_);

        if (!resultPtr)
        {
            resultPtr = new volScalarField
            (
                IOobject
                (
                    resultName_,
                    mesh_.time().timeName(),
                    mesh_.thisDb(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::REGISTER
                ),
                mesh_,
                dimensionedScalar(dimless, Zero),
                fvPatchFieldBase::zeroGradientType()
            );
            regIOobject::store(resultPtr);
        }
        auto& Co = *resultPtr;

        Co.internalFieldRef() = Coi;
        Co.correctBoundaryConditions();

        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::CourantNo::CourantNo
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, "phi"),
    rhoName_("rho")
{
    setResultName("Co", "phi");
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::CourantNo::read(const dictionary& dict)
{
    fieldExpression::read(dict);

    rhoName_ = dict.getOrDefault<word>("rho", "rho");

    return true;
}


// ************************************************************************* //
