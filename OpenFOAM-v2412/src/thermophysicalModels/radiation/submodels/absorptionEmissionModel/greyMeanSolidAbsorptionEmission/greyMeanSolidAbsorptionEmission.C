/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015 OpenCFD Ltd.
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

#include "greyMeanSolidAbsorptionEmission.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"
#include "extrapolatedCalculatedFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(greyMeanSolidAbsorptionEmission, 0);

        addToRunTimeSelectionTable
        (
            absorptionEmissionModel,
            greyMeanSolidAbsorptionEmission,
            dictionary
        );
    }
}

// * * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::radiation::
greyMeanSolidAbsorptionEmission::X(const word specie) const
{
    const volScalarField& T = thermo_.T();
    const volScalarField& p = thermo_.p();

    auto tXj = tmp<scalarField>::New(T.primitiveField().size(), Zero);
    auto& Xj = tXj.ref();

    auto tRhoInv = tmp<scalarField>::New(T.primitiveField().size(), Zero);
    auto& rhoInv = tRhoInv.ref();

    forAll(mixture_.Y(), specieI)
    {
        const scalarField& Yi = mixture_.Y()[specieI];

        forAll(rhoInv, iCell)
        {
            rhoInv[iCell] +=
                Yi[iCell]/mixture_.rho(specieI, p[iCell], T[iCell]);
        }
    }
    const scalarField& Yj = mixture_.Y(specie);
    const label mySpecieI = mixture_.species().find(specie);
    forAll(Xj, iCell)
    {
        Xj[iCell] = Yj[iCell]/mixture_.rho(mySpecieI, p[iCell], T[iCell]);
    }

    return (Xj/rhoInv);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::greyMeanSolidAbsorptionEmission::
greyMeanSolidAbsorptionEmission
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    absorptionEmissionModel(dict, mesh),
    coeffsDict_((dict.optionalSubDict(typeName + "Coeffs"))),
    thermo_(mesh.lookupObject<solidThermo>(basicThermo::dictName)),
    speciesNames_(0),
    mixture_(dynamic_cast<const basicSpecieMixture&>(thermo_)),
    solidData_(mixture_.Y().size())
{
    if (!isA<basicSpecieMixture>(thermo_))
    {
        FatalErrorInFunction
            << "Model requires a multi-component thermo package"
            << abort(FatalError);
    }

    label nFunc = 0;
    const dictionary& functionDicts = dict.optionalSubDict(typeName + "Coeffs");

    for (const entry& dEntry : functionDicts)
    {
        if (!dEntry.isDict())  // safety
        {
            continue;
        }

        const word& key = dEntry.keyword();
        const dictionary& dict = dEntry.dict();

        if (!mixture_.contains(key))
        {
            WarningInFunction
                << " specie: " << key << " is not found in the solid mixture"
                << nl
                << " specie is the mixture are:" << mixture_.species() << nl
                << nl << endl;
        }
        speciesNames_.insert(key, nFunc);

        dict.readEntry("absorptivity", solidData_[nFunc][absorptivity]);
        dict.readEntry("emissivity", solidData_[nFunc][emissivity]);

        nFunc++;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiation::greyMeanSolidAbsorptionEmission::
calc(const label propertyId) const
{
    auto ta = volScalarField::New
    (
        "a",
        IOobject::NO_REGISTER,
        mesh(),
        dimensionedScalar(dimless/dimLength, Zero),
        fvPatchFieldBase::extrapolatedCalculatedType()
    );
    auto& a = ta.ref().primitiveFieldRef();

    forAllConstIters(speciesNames_, iter)
    {
        if (mixture_.contains(iter.key()))
        {
            a += solidData_[iter.val()][propertyId]*X(iter.key());
        }
    }

    ta.ref().correctBoundaryConditions();
    return ta;
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::greyMeanSolidAbsorptionEmission::eCont
(
    const label bandI
) const
{
   return calc(emissivity);
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::greyMeanSolidAbsorptionEmission::aCont
(
    const label bandI
) const
{
   return calc(absorptivity);
}


// ************************************************************************* //
