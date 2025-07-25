/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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

#include "limitTurbulenceViscosity.H"
#include "fvMesh.H"
#include "transportModel.H"
#include "fluidThermo.H"
#include "turbulenceModel.H"
#include "processorFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(limitTurbulenceViscosity, 0);
    addToRunTimeSelectionTable(option, limitTurbulenceViscosity, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fv::limitTurbulenceViscosity::writeFileHeader(Ostream& os)
{
    writeHeaderValue(os, "Nut", nutName_);
    writeHeaderValue(os, "c", c_);
    writeCommented(os, "Time");
    writeTabbed(os, "nLimitedCells_[count]");
    writeTabbed(os, "nLimitedCells_[%]");

    os  << endl;

    writtenHeader_ = true;
}


Foam::tmp<Foam::volScalarField> Foam::fv::limitTurbulenceViscosity::nu() const
{
    const auto* turbPtr =
        mesh_.cfindObject<turbulenceModel>(turbulenceModel::propertiesName);
    if (turbPtr)
    {
        return turbPtr->nu();
    }

    const auto* thermoPtr =
        mesh_.cfindObject<fluidThermo>(fluidThermo::dictName);
    if (thermoPtr)
    {
        return thermoPtr->nu();
    }

    const auto* laminarPtr =
        mesh_.cfindObject<transportModel>("transportProperties");
    if (laminarPtr)
    {
        return laminarPtr->nu();
    }

    const auto* dictPtr = mesh_.cfindObject<dictionary>("transportProperties");
    if (dictPtr)
    {
        return volScalarField::New
        (
            "nu",
            IOobject::NO_REGISTER,
            mesh_,
            dimensionedScalar("nu", dimViscosity, *dictPtr)
        );
    }

    FatalErrorInFunction
        << "No valid model for laminar viscosity"
        << exit(FatalError);

    return nullptr;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::limitTurbulenceViscosity::limitTurbulenceViscosity
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::cellSetOption(name, modelType, dict, mesh),
    writeFile(mesh, name, typeName, dict, false),
    nutName_("nut"),
    c_(1e5)
{
    if (isActive())
    {
        read(dict);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::limitTurbulenceViscosity::read(const dictionary& dict)
{
    if (!(fv::cellSetOption::read(dict) && writeFile::read(dict)))
    {
        return false;
    }

    coeffs_.readIfPresent("nut", nutName_);
    coeffs_.readIfPresent("c", c_);

    fieldNames_.resize(1, nutName_);

    fv::option::resetApplied();

    if (canResetFile())
    {
        resetFile(typeName);
    }

    if (canWriteHeader())
    {
        writeFileHeader(file());
    }


    return true;
}


void Foam::fv::limitTurbulenceViscosity::correct(volScalarField& nut)
{
    const auto& tnu = this->nu();
    const auto& nu = tnu();

    const label nTotCells(returnReduce(cells_.size(), sumOp<label>()));

    label nAboveMax = 0;
    for (const label celli : cells_)
    {
        const scalar nutLim = c_*nu[celli];

        if (nut[celli] > nutLim)
        {
            nut[celli] = nutLim;
            ++nAboveMax;
        }
    }

    reduce(nAboveMax, sumOp<label>());

    // Percent, max 2 decimal places
    const auto percent = [](scalar num, label denom) -> scalar
    {
        return (denom ? 1e-2*round(1e4*num/denom) : 0);
    };

    const scalar nAboveMaxPercent = percent(nAboveMax, nTotCells);

    Info<< type() << ' ' << name_ << " limited " << nAboveMax << " ("
        << nAboveMaxPercent << "%) of cells" << endl;

    if (canWriteToFile())
    {
        file()
            << mesh_.time().timeOutputValue() << token::TAB
            << nAboveMax << token::TAB
            << nAboveMaxPercent
            << endl;
    }


    // Handle boundaries in the case of 'all'
    if (!cellSetOption::useSubMesh())
    {
        const volScalarField::Boundary& nubf = nu.boundaryField();
        volScalarField::Boundary& nutbf = nut.boundaryFieldRef();

        for (auto& nutp : nutbf)
        {
            // if (!nutp.fixesValue())
            {
                const scalarField& nup = nubf[nutp.patch().index()];

                forAll(nutp, facei)
                {
                    scalar nutLim = c_*nup[facei];
                    if (nutp[facei] > nutLim)
                    {
                        nutp[facei] = nutLim;
                    }
                }
            }
        }
    }


    if (nAboveMax)
    {
        // We've changed internal values so give boundary conditions
        // opportunity to correct

        // Note: calling nut.correctBoundaryConditions() will re-evaluate wall
        // functions and potentially (re-)introduce out-of-bounds values
        //nut.correctBoundaryConditions();
        if (UPstream::parRun())
        {
            nut.boundaryFieldRef().evaluateCoupled<processorFvPatch>();
        }
    }
}


// ************************************************************************* //
