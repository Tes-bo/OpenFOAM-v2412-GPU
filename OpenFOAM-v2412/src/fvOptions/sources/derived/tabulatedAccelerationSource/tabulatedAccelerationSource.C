/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2016 OpenFOAM Foundation
    Copyright (C) 2020-2021 OpenCFD Ltd.
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

#include "tabulatedAccelerationSource.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "geometricOneField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(tabulatedAccelerationSource, 0);
    addToRunTimeSelectionTable(option, tabulatedAccelerationSource, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::tabulatedAccelerationSource::tabulatedAccelerationSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::option(name, modelType, dict, mesh),
    motion_(coeffs_, mesh.time()),
    UName_(coeffs_.getOrDefault<word>("U", "U")),
    g0_("g0", dimAcceleration, Zero)
{
    fieldNames_.resize(1, UName_);
    fv::option::resetApplied();

    if (mesh.time().foundObject<uniformDimensionedVectorField>("g"))
    {
        g0_ = mesh.time().lookupObject<uniformDimensionedVectorField>("g");
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::tabulatedAccelerationSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    addSup<geometricOneField>(geometricOneField(), eqn, fieldi);
}


void Foam::fv::tabulatedAccelerationSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    addSup<volScalarField>(rho, eqn, fieldi);
}


bool Foam::fv::tabulatedAccelerationSource::read(const dictionary& dict)
{
    if (fv::option::read(dict))
    {
        return motion_.read(coeffs_);
    }

    return false;
}


// ************************************************************************* //
