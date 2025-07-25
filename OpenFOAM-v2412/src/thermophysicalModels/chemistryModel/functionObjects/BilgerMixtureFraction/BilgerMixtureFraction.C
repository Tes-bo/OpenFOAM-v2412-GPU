/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 Thorsten Zirwes
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

#include "BilgerMixtureFraction.H"
#include "basicThermo.H"
#include "reactingMixture.H"
#include "thermoPhysicsTypes.H"
#include "scalarRange.H"
#include "basicChemistryModel.H"
#include "psiReactionThermo.H"
#include "rhoReactionThermo.H"
#include "BasicChemistryModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(BilgerMixtureFraction, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        BilgerMixtureFraction,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::BilgerMixtureFraction::calcBilgerMixtureFraction()
{
    auto* resultPtr = mesh_.getObjectPtr<volScalarField>(resultName_);

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
            dimensionedScalar(dimless, Zero)
        );
        regIOobject::store(resultPtr);
    }
    auto& f_Bilger = *resultPtr;

    auto& Y = thermo_.Y();

    f_Bilger = -o2RequiredOx_;
    forAll(Y, i)
    {
        f_Bilger +=
            Y[i]
           *(nAtomsC_[i] + nAtomsS_[i] + 0.25*nAtomsH_[i] - 0.5*nAtomsO_[i])
           /thermo_.W(i);
    }
    f_Bilger /= o2RequiredFuelOx_;
    f_Bilger.clamp_range(zero_one{});
}


bool Foam::functionObjects::BilgerMixtureFraction::readComposition
(
    const dictionary& subDict,
    scalarField& comp
)
{
    auto& Y = thermo_.Y();
    const speciesTable& speciesTab = thermo_.species();

    // Read mass fractions of all species for the oxidiser or fuel
    forAll(Y, i)
    {
        comp[i] =
            subDict.getCheckOrDefault<scalar>
            (
                speciesTab[i],
                0,
                scalarRange::ge0()
            );
    }

    if (sum(comp) < SMALL)
    {
        FatalIOErrorInFunction(subDict)
            << "No composition is given" << nl
            << "Valid species are:" << nl
            << speciesTab
            << exit(FatalIOError);

        return false;
    }

    const word fractionBasisType
    (
        subDict.getOrDefault<word>("fractionBasis", "mass")
    );

    if (fractionBasisType == "mass")
    {
        // Normalize fractionBasis to the unity
        comp /= sum(comp);
    }
    else if (fractionBasisType == "mole")
    {
        // In case the fractionBasis is given in mole fractions,
        // convert from mole fractions to normalized mass fractions
        scalar W(0);
        forAll(comp, i)
        {
            comp[i] *= thermo_.W(i);
            W += comp[i];
        }
        comp /= W;
    }
    else
    {
        FatalIOErrorInFunction(subDict)
            << "The given fractionBasis type is invalid" << nl
            << "Valid fractionBasis types are" << nl
            << "  \"mass\" (default)" << nl
            << "  \"mole\""
            << exit(FatalIOError);

        return false;
    }

    return true;
}


Foam::scalar Foam::functionObjects::BilgerMixtureFraction::o2Required
(
    const scalarField& comp
) const
{
    scalar o2req(0);
    forAll(thermo_.Y(), i)
    {
        o2req +=
            comp[i]/thermo_.W(i)*(nAtomsC_[i] + nAtomsS_[i] + 0.25*nAtomsH_[i]);
    }

    return o2req;
}


Foam::scalar Foam::functionObjects::BilgerMixtureFraction::o2Present
(
    const scalarField& comp
) const
{
    scalar o2pres(0);
    forAll(thermo_.Y(), i)
    {
        o2pres += comp[i]/thermo_.W(i)*nAtomsO_[i];
    }

    return 0.5*o2pres;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::BilgerMixtureFraction::BilgerMixtureFraction
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    phaseName_(dict.getOrDefault<word>("phase", word::null)),
    resultName_
    (
        dict.getOrDefault<word>
        (
            "result",
            IOobject::groupName("f_Bilger", phaseName_)
        )
    ),
    thermo_
    (
        mesh_.lookupObject<basicSpecieMixture>
        (
            IOobject::groupName(basicThermo::dictName, phaseName_)
        )
    ),
    nSpecies_(thermo_.Y().size()),
    o2RequiredOx_(0),
    o2RequiredFuelOx_(0),
    nAtomsC_(nSpecies_, 0),
    nAtomsS_(nSpecies_, 0),
    nAtomsH_(nSpecies_, 0),
    nAtomsO_(nSpecies_, 0),
    Yoxidiser_(nSpecies_, 0),
    Yfuel_(nSpecies_, 0)
{
    read(dict);

    calcBilgerMixtureFraction();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::BilgerMixtureFraction::read(const dictionary& dict)
{
    if (!fvMeshFunctionObject::read(dict))
    {
        return false;
    }

    Info<< nl << type() << " " << name() << ":" << nl;

    phaseName_ = dict.getOrDefault<word>("phase", word::null);
    resultName_ =
        dict.getOrDefault<word>
        (
            "result",
            IOobject::groupName("f_Bilger", phaseName_)
        );

    nSpecies_ = thermo_.Y().size();

    if (nSpecies_ == 0)
    {
        FatalErrorInFunction
            << "Number of input species is zero"
            << exit(FatalError);
    }

    nAtomsC_.resize(nSpecies_, 0);
    nAtomsS_.resize(nSpecies_, 0);
    nAtomsH_.resize(nSpecies_, 0);
    nAtomsO_.resize(nSpecies_, 0);

    auto& Y = thermo_.Y();
    const speciesTable& speciesTab = thermo_.species();

    typedef BasicChemistryModel<psiReactionThermo> psiChemistryModelType;
    typedef BasicChemistryModel<rhoReactionThermo> rhoChemistryModelType;

    const auto* psiChemPtr =
        mesh_.cfindObject<psiChemistryModelType>("chemistryProperties");

    const auto* rhoChemPtr =
        mesh_.cfindObject<rhoChemistryModelType>("chemistryProperties");

    autoPtr<speciesCompositionTable> speciesCompPtr;

    if (psiChemPtr)
    {
        speciesCompPtr.reset((*psiChemPtr).thermo().specieComposition());
    }
    else if (rhoChemPtr)
    {
        speciesCompPtr.reset((*rhoChemPtr).thermo().specieComposition());
    }
    else
    {
        FatalErrorInFunction
            << "BasicChemistryModel not found"
            << exit(FatalError);
    }

    forAll(Y, i)
    {
        const List<specieElement>& curSpecieComposition =
            (speciesCompPtr.ref())[speciesTab[i]];

        forAll(curSpecieComposition, j)
        {
            const word& e = curSpecieComposition[j].name();
            const label nAtoms  = curSpecieComposition[j].nAtoms();

            if (e == "C")
            {
                nAtomsC_[i] = nAtoms;
            }
            else if (e == "S")
            {
                nAtomsS_[i] = nAtoms;
            }
            else if (e == "H")
            {
                nAtomsH_[i] = nAtoms;
            }
            else if (e == "O")
            {
                nAtomsO_[i] = nAtoms;
            }
        }
    }

    if (sum(nAtomsO_) == 0)
    {
        FatalErrorInFunction
            << "No specie contains oxygen"
            << exit(FatalError);
    }

    Yoxidiser_.resize(nSpecies_, 0);
    Yfuel_.resize(nSpecies_, 0);

    if
    (
        !readComposition(dict.subDict("oxidiser"), Yoxidiser_)
     || !readComposition(dict.subDict("fuel"), Yfuel_)
    )
    {
        return false;
    }

    o2RequiredOx_ = o2Required(Yoxidiser_) - o2Present(Yoxidiser_);

    if (o2RequiredOx_ > 0)
    {
        FatalErrorInFunction
            << "Oxidiser composition contains not enough oxygen" << endl
            << "Mixed up fuel and oxidiser compositions?"
            << exit(FatalError);
    }

    const scalar o2RequiredFuel = o2Required(Yfuel_) - o2Present(Yfuel_);

    if (o2RequiredFuel < 0)
    {
        FatalErrorInFunction
            << "Fuel composition contains too much oxygen" << endl
            << "Mixed up fuel and oxidiser compositions?"
            << exit(FatalError);
    }

    o2RequiredFuelOx_ = o2RequiredFuel - o2RequiredOx_;

    if (mag(o2RequiredFuelOx_) < SMALL)
    {
        FatalErrorInFunction
            << "Fuel and oxidiser have the same composition"
            << exit(FatalError);
    }

    return true;
}


bool Foam::functionObjects::BilgerMixtureFraction::execute()
{
    calcBilgerMixtureFraction();

    return true;
}


bool Foam::functionObjects::BilgerMixtureFraction::clear()
{
    return clearObject(resultName_);
}


bool Foam::functionObjects::BilgerMixtureFraction::write()
{
    Log << type() << " " << name() << " write:" << nl
        << "    writing field " << resultName_ << endl;

    return writeObject(resultName_);
}


// ************************************************************************* //
