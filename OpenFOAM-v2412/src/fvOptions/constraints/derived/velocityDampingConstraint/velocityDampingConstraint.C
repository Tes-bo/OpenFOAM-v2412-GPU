/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

#include "velocityDampingConstraint.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMatrix.H"
#include "volFields.H"
#include "cellCellStencilObject.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(velocityDampingConstraint, 0);
    addToRunTimeSelectionTable
    (
        option,
        velocityDampingConstraint,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fv::velocityDampingConstraint::addDamping(fvMatrix<vector>& eqn)
{
    // Note: we want to add
    //      deltaU/deltaT
    // where deltaT is a local time scale:
    //      U/(cbrt of volume)
    // Since directly manipulating the diagonal we multiply by volume.

    const scalarField& vol = mesh_.V();
    const volVectorField& U = eqn.psi();
    scalarField& diag = eqn.diag();

    // Count nTotCells ourselves
    // (maybe only applying on a subset)
    label nDamped(0);
    const label nTotCells(returnReduce(cells_.size(), sumOp<label>()));

    for (const label celli : cells_)
    {
        const scalar magU = mag(U[celli]);
        if (magU > UMax_)
        {
            const scalar scale = sqr(Foam::cbrt(vol[celli]));

            diag[celli] += C_*scale*(magU-UMax_);

            ++nDamped;
        }
    }

    reduce(nDamped, sumOp<label>());

    // Percent, max 2 decimal places
    const auto percent = [](scalar num, label denom) -> scalar
    {
        return (denom ? 1e-2*round(1e4*num/denom) : 0);
    };

    const scalar nDampedPercent = percent(nDamped, nTotCells);

    Info<< type() << ' ' << name_ << " damped "
        << nDamped << " ("
        << nDampedPercent
        << "%) of cells, with max limit " << UMax_ << endl;


    if (canWriteToFile())
    {
        file()
            << mesh_.time().timeOutputValue() << token::TAB
            << nDamped << token::TAB
            << nDampedPercent
            << endl;
    }
}


void Foam::fv::velocityDampingConstraint::writeFileHeader(Ostream& os)
{
    writeHeaderValue(os, "UMax", Foam::name(UMax_));
    writeCommented(os, "Time");
    writeTabbed(os, "nDamped_[count]");
    writeTabbed(os, "nDamped_[%]");

    os  << endl;

    writtenHeader_ = true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::velocityDampingConstraint::velocityDampingConstraint
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::cellSetOption(name, modelType, dict, mesh),
    writeFile(mesh, name, typeName, dict, false),
    UMax_(GREAT),  // overwritten later
    C_(1)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::velocityDampingConstraint::constrain
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    addDamping(eqn);
}


void Foam::fv::velocityDampingConstraint::writeData(Ostream& os) const
{
    dict_.writeEntry(name_, os);
}


bool Foam::fv::velocityDampingConstraint::read(const dictionary& dict)
{
    if (!(fv::cellSetOption::read(dict) && writeFile::read(dict)))
    {
        return false;
    }

    coeffs_.readEntry("UMax", UMax_);
    coeffs_.readIfPresent("C", C_);

    if (!coeffs_.readIfPresent("UNames", fieldNames_))
    {
        fieldNames_.resize(1);
        fieldNames_.first() = coeffs_.getOrDefault<word>("U", "U");
    }

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


// ************************************************************************* //
