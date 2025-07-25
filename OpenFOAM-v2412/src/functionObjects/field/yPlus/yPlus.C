/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2016-2024 OpenCFD Ltd.
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

#include "yPlus.H"
#include "volFields.H"
#include "turbulenceModel.H"
#include "nutWallFunctionFvPatchScalarField.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(yPlus, 0);
    addToRunTimeSelectionTable(functionObject, yPlus, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::yPlus::writeFileHeader(Ostream& os) const
{
    writeHeader(os, "y+ ()");

    writeCommented(os, "Time");
    writeTabbed(os, "patch");
    writeTabbed(os, "min");
    writeTabbed(os, "max");
    writeTabbed(os, "average");
    os  << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::yPlus::yPlus
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(obr_, name, typeName, dict),
    useWallFunction_(true),
    writeFields_(true)  // May change in the future
{
    read(dict);

    writeFileHeader(file());

    volScalarField* yPlusPtr
    (
        new volScalarField
        (
            IOobject
            (
                scopedName(typeName),
                mesh_.time().timeName(),
                mesh_.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::REGISTER
            ),
            mesh_,
            dimensionedScalar(dimless, Zero)
        )
    );

    regIOobject::store(yPlusPtr);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::yPlus::read(const dictionary& dict)
{
    if (fvMeshFunctionObject::read(dict) && writeFile::read(dict))
    {
        useWallFunction_ = true;
        writeFields_ = true;   // May change in the future

        dict.readIfPresent("useWallFunction", useWallFunction_);
        dict.readIfPresent("writeFields", writeFields_);

        return true;
    }

    return false;
}


bool Foam::functionObjects::yPlus::execute()
{
    auto& yPlus = lookupObjectRef<volScalarField>(scopedName(typeName));

    if (foundObject<turbulenceModel>(turbulenceModel::propertiesName))
    {
        volScalarField::Boundary& yPlusBf = yPlus.boundaryFieldRef();

        const turbulenceModel& model =
            lookupObject<turbulenceModel>
            (
                turbulenceModel::propertiesName
            );

        const nearWallDist nwd(mesh_);
        const volScalarField::Boundary& d = nwd.y();

        // nut needed for wall function patches
        tmp<volScalarField> tnut = model.nut();
        const volScalarField::Boundary& nutBf = tnut().boundaryField();

        // U needed for plain wall patches
        const volVectorField::Boundary& UBf = model.U().boundaryField();

        const fvPatchList& patches = mesh_.boundary();

        forAll(patches, patchi)
        {
            const fvPatch& patch = patches[patchi];

            const auto* nutWallPatch =
                isA<nutWallFunctionFvPatchScalarField>(nutBf[patchi]);

            if (useWallFunction_ && nutWallPatch)
            {
                yPlusBf[patchi] = nutWallPatch->yPlus();
            }
            else if (isA<wallFvPatch>(patch))
            {
                yPlusBf[patchi] =
                    d[patchi]
                   *sqrt
                    (
                        model.nuEff(patchi)
                       *mag(UBf[patchi].snGrad())
                    )/model.nu(patchi);
            }
        }
    }
    else
    {
        WarningInFunction
            << "Unable to find turbulence model in the "
            << "database: yPlus will not be calculated" << endl;

        if (postProcess)
        {
            WarningInFunction
                << "Please try to use the solver option -postProcess, e.g.:"
                << " <solver> -postProcess -func yPlus" << endl;
        }

        return false;
    }

    return true;
}


bool Foam::functionObjects::yPlus::write()
{
    const auto& yPlus = obr_.lookupObject<volScalarField>(scopedName(typeName));

    Log << type() << ' ' << name() << " write:" << nl;

    if (writeFields_)
    {
        Log << "    writing field " << yPlus.name() << endl;
        yPlus.write();
    }

    const volScalarField::Boundary& yPlusBf = yPlus.boundaryField();
    const fvPatchList& patches = mesh_.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& patch = patches[patchi];

        if (isA<wallFvPatch>(patch))
        {
            const scalarField& yPlusp = yPlusBf[patchi];

            const scalar minYplus = gMin(yPlusp);
            const scalar maxYplus = gMax(yPlusp);
            const scalar avgYplus = gAverage(yPlusp);

            if (UPstream::master())
            {
                writeCurrentTime(file());
                file()
                    << token::TAB << patch.name()
                    << token::TAB << minYplus
                    << token::TAB << maxYplus
                    << token::TAB << avgYplus
                    << endl;
            }

            Log << "    patch " << patch.name()
                << " y+ : min = " << minYplus
                << ", max = " << maxYplus
                << ", average = " << avgYplus << endl;
        }
    }

    return true;
}


// ************************************************************************* //
