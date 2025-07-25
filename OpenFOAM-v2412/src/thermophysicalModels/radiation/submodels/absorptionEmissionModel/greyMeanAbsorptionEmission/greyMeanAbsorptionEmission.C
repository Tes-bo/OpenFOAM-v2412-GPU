/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "greyMeanAbsorptionEmission.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"
#include "extrapolatedCalculatedFvPatchFields.H"
#include "basicSpecieMixture.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(greyMeanAbsorptionEmission, 0);

        addToRunTimeSelectionTable
        (
            absorptionEmissionModel,
            greyMeanAbsorptionEmission,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::greyMeanAbsorptionEmission::greyMeanAbsorptionEmission
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    absorptionEmissionModel(dict, mesh),
    coeffsDict_((dict.optionalSubDict(typeName + "Coeffs"))),
    speciesNames_(0),
    specieIndex_(Zero),
    lookUpTablePtr_(),
    thermo_(mesh.lookupObject<fluidThermo>(basicThermo::dictName)),
    EhrrCoeff_(coeffsDict_.get<scalar>("EhrrCoeff")),
    Yj_(nSpecies_)
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

        speciesNames_.insert(key, nFunc);

        coeffs_[nFunc].initialise(dict);
        nFunc++;
    }

    if
    (
        coeffsDict_.found("lookUpTableFileName")
     && "none" != coeffsDict_.get<word>("lookUpTableFileName")
    )
    {
        lookUpTablePtr_.reset
        (
            new interpolationLookUpTable<scalar>
            (
                coeffsDict_.get<fileName>("lookUpTableFileName"),
                mesh.time().constant(),
                mesh
            )
        );

        if (!mesh.foundObject<volScalarField>("ft"))
        {
            FatalErrorInFunction
                << "specie ft is not present to use with "
                << "lookUpTableFileName " << nl
                << exit(FatalError);
        }
    }

    // Check that all the species on the dictionary are present in the
    // look-up table and save the corresponding indices of the look-up table

    label j = 0;
    forAllConstIters(speciesNames_, iter)
    {
        const word& specieName = iter.key();
        const label index = iter.val();

        volScalarField* fldPtr = mesh.getObjectPtr<volScalarField>(specieName);

        if (lookUpTablePtr_)
        {
            if (lookUpTablePtr_().found(specieName))
            {
                const label fieldIndex =
                    lookUpTablePtr_().findFieldIndex(specieName);

                Info<< "specie: " << specieName << " found on look-up table "
                    << " with index: " << fieldIndex << endl;

                specieIndex_[index] = fieldIndex;
            }
            else if (fldPtr)
            {
                Yj_.set(j, fldPtr);
                specieIndex_[index] = 0;
                j++;
                Info<< "specie: " << iter.key() << " is being solved" << endl;
            }
            else
            {
                FatalErrorInFunction
                    << "specie: " << specieName
                    << " is neither in look-up table: "
                    << lookUpTablePtr_().tableName()
                    << " nor is being solved" << nl
                    << exit(FatalError);
            }
        }
        else if (fldPtr)
        {
            Yj_.set(j, fldPtr);
            specieIndex_[index] = 0;
            j++;
        }
        else
        {
            FatalErrorInFunction
                << "There is no lookup table and the specie" << nl
                << specieName << nl
                << " is not found " << nl
                << exit(FatalError);
        }
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::greyMeanAbsorptionEmission::~greyMeanAbsorptionEmission()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiation::greyMeanAbsorptionEmission::aCont(const label bandI) const
{
    const basicSpecieMixture& mixture =
        dynamic_cast<const basicSpecieMixture&>(thermo_);

    const volScalarField& T = thermo_.T();
    const volScalarField& p = thermo_.p();


    auto ta = volScalarField::New
    (
        "aCont" + name(bandI),
        IOobject::NO_REGISTER,
        mesh(),
        dimensionedScalar(dimless/dimLength, Zero),
        fvPatchFieldBase::extrapolatedCalculatedType()
    );
    auto& a = ta.ref().primitiveFieldRef();

    forAll(a, celli)
    {
        forAllConstIters(speciesNames_, iter)
        {
            label n = iter();
            scalar Xipi = 0.0;
            if (specieIndex_[n] != 0)
            {
                //Specie found in the lookUpTable.
                const volScalarField& ft =
                    mesh_.lookupObject<volScalarField>("ft");

                const List<scalar>& Ynft = lookUpTablePtr_().lookUp(ft[celli]);
                //moles x pressure [atm]
                Xipi = Ynft[specieIndex_[n]]*paToAtm(p[celli]);
            }
            else
            {
                scalar invWt = 0.0;
                forAll(mixture.Y(), s)
                {
                    invWt += mixture.Y(s)[celli]/mixture.W(s);
                }

                const label index = mixture.species().find(iter.key());
                scalar Xk = mixture.Y(index)[celli]/(mixture.W(index)*invWt);

                Xipi = Xk*paToAtm(p[celli]);
            }

            const absorptionCoeffs::coeffArray& b = coeffs_[n].coeffs(T[celli]);

            scalar Ti = T[celli];
            // negative temperature exponents
            if (coeffs_[n].invTemp())
            {
                Ti = 1.0/T[celli];
            }
            a[celli] +=
                Xipi
               *(
                    ((((b[5]*Ti + b[4])*Ti + b[3])*Ti + b[2])*Ti + b[1])*Ti
                  + b[0]
                );
        }
    }
    ta.ref().correctBoundaryConditions();
    return ta;
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::greyMeanAbsorptionEmission::eCont(const label bandI) const
{
   return aCont(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::greyMeanAbsorptionEmission::ECont(const label bandI) const
{
    auto E = volScalarField::New
    (
        "ECont" + name(bandI),
        IOobject::NO_REGISTER,
        mesh_,
        dimensionedScalar(dimMass/dimLength/pow3(dimTime), Zero)
    );

    const volScalarField* QdotPtr = mesh_.findObject<volScalarField>("Qdot");

    if (QdotPtr)
    {
        const volScalarField& Qdot = *QdotPtr;

        if (Qdot.dimensions() == dimEnergy/dimTime)
        {
            E.ref().primitiveFieldRef() = EhrrCoeff_*Qdot/mesh_.V();
        }
        else if (Qdot.dimensions() == dimEnergy/dimTime/dimVolume)
        {
            E.ref().primitiveFieldRef() = EhrrCoeff_*Qdot;
        }
        else
        {
            if (debug)
            {
                WarningInFunction
                    << "Incompatible dimensions for Qdot field" << endl;
            }
        }
    }
    else
    {
        WarningInFunction
          << "Qdot field not found in mesh" << endl;
    }

    return E;
}


// ************************************************************************* //
