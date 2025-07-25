/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenFOAM Foundation
    Copyright (C) 2018-2023 OpenCFD Ltd.
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

#include "limitVelocity.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(limitVelocity, 0);
    addToRunTimeSelectionTable(option, limitVelocity, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fv::limitVelocity::writeFileHeader(Ostream& os)
{
    writeHeaderValue(os, "UMax", Foam::name(max_));
    writeCommented(os, "Time");
    writeTabbed(os, "nDampedCells_[count]");
    writeTabbed(os, "nDampedCells_[%]");
    writeTabbed(os, "nDampedFaces_[count]");
    writeTabbed(os, "nDampedFaces_[%]");

    os  << endl;

    writtenHeader_ = true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::limitVelocity::limitVelocity
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::cellSetOption(name, modelType, dict, mesh),
    writeFile(mesh, name, typeName, dict, false),
    UName_("U"),
    max_(0)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::limitVelocity::read(const dictionary& dict)
{
    if (!(fv::cellSetOption::read(dict) && writeFile::read(dict)))
    {
        return false;
    }

    coeffs_.readEntry("max", max_);
    coeffs_.readIfPresent("U", UName_);


    fieldNames_.resize(1, UName_);

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


void Foam::fv::limitVelocity::correct(volVectorField& U)
{
    const scalar maxSqrU = sqr(max_);

    // Count nTotCells ourselves
    // (maybe only applying on a subset)
    label nCellsAbove(0);
    const label nTotCells(returnReduce(cells_.size(), sumOp<label>()));

    vectorField& Uif = U.primitiveFieldRef();

    for (const label celli : cells_)
    {
        auto& Uval = Uif[celli];

        const scalar magSqrUi = magSqr(Uval);

        if (magSqrUi > maxSqrU)
        {
            Uval *= sqrt(maxSqrU/magSqrUi);
            ++nCellsAbove;
        }
    }

    // Handle boundaries in the case of 'all'

    label nFacesAbove(0);
    label nTotFaces(0);

    if (!cellSetOption::useSubMesh())
    {
        for (fvPatchVectorField& Up : U.boundaryFieldRef())
        {
            if (!Up.fixesValue())
            {
                // Do not count patches that fix velocity themselves
                nTotFaces += Up.size();

                for (auto& Uval : Up)
                {
                    const scalar magSqrUi = magSqr(Uval);

                    if (magSqrUi > maxSqrU)
                    {
                        Uval *= sqrt(maxSqrU/magSqrUi);
                        ++nFacesAbove;
                    }
                }
            }
        }
    }

    // Percent, max 2 decimal places
    const auto percent = [](scalar num, label denom) -> scalar
    {
        return (denom ? 1e-2*round(1e4*num/denom) : 0);
    };


    reduce(nCellsAbove, sumOp<label>());

    const scalar nCellsAbovePercent = percent(nCellsAbove, nTotCells);

    // Report total numbers and percent
    Info<< type() << ' ' << name_ << " Limited ";

    Info<< nCellsAbove << " ("
        << nCellsAbovePercent
        << "%) of cells";

    reduce(nTotFaces, sumOp<label>());
    reduce(nFacesAbove, sumOp<label>());
    scalar nFacesAbovePercent(0);
    if (nTotFaces)
    {
        nFacesAbovePercent = percent(nFacesAbove, nTotFaces);

        Info<< ", " << nFacesAbove << " ("
            << nFacesAbovePercent
            << "%) of faces";
    }
    Info<< ", with max limit " << max_ << endl;

    if (nCellsAbove || nFacesAbove)
    {
        // We've changed internal values so give
        // boundary conditions opportunity to correct
        U.correctBoundaryConditions();
    }


    if (canWriteToFile())
    {
        file()
            << mesh_.time().timeOutputValue() << token::TAB
            << nCellsAbove << token::TAB
            << nCellsAbovePercent << token::TAB
            << nFacesAbove << token::TAB
            << nFacesAbovePercent
            << endl;
    }
}


// ************************************************************************* //
