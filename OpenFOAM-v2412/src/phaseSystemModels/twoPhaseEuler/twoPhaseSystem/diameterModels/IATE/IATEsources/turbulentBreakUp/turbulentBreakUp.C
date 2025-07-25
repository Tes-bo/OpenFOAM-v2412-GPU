/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2015 OpenFOAM Foundation
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

#include "turbulentBreakUp.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace IATEsources
{
    defineTypeNameAndDebug(turbulentBreakUp, 0);
    addToRunTimeSelectionTable(IATEsource, turbulentBreakUp, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::IATEsources::turbulentBreakUp::
turbulentBreakUp
(
    const IATE& iate,
    const dictionary& dict
)
:
    IATEsource(iate),
    Cti_("Cti", dimless, dict),
    WeCr_("WeCr", dimless, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::diameterModels::IATEsources::turbulentBreakUp::R() const
{
    auto tR = volScalarField::New
    (
        "R",
        IOobject::NO_REGISTER,
        iate_.phase().U().mesh(),
        dimensionedScalar(dimless/dimTime, Zero)
    );
    auto R = tR();

    const scalar Cti = Cti_.value();
    const scalar WeCr = WeCr_.value();
    volScalarField Ut(this->Ut());
    volScalarField We(this->We());
    const tmp<volScalarField> td(iate_.d());
    const volScalarField& d = td();

    forAll(R, celli)
    {
        if (We[celli] > WeCr)
        {
            R[celli] =
                (1.0/3.0)
               *Cti/d[celli]
               *Ut[celli]
               *sqrt(1 - WeCr/We[celli])
               *exp(-WeCr/We[celli]);
        }
    }

    return tR;
}


// ************************************************************************* //
