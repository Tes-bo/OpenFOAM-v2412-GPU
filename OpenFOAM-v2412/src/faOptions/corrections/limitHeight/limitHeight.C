/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "limitHeight.H"
#include "areaFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fa
{
    defineTypeNameAndDebug(limitHeight, 0);
    addToRunTimeSelectionTable
    (
        option,
        limitHeight,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fa::limitHeight::limitHeight
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    faceSetOption(name, modelType, dict, mesh),
    hName_("h"),
    max_(0)  // overwritten later
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fa::limitHeight::read(const dictionary& dict)
{
    if (!faceSetOption::read(dict))
    {
        return false;
    }

    coeffs_.readIfPresent("h", hName_);
    coeffs_.readEntry("max", max_);

    fieldNames_.resize(1, hName_);

    applied_.resize(1, false);

    return true;
}


void Foam::fa::limitHeight::correct(areaScalarField& h)
{
    // Count nTotFaces ourselves
    // (maybe only applying on a subset)
    label nFacesAbove = 0;
    const label nTotFaces = returnReduce(faces_.size(), sumOp<label>());

    scalarField& hif = h.primitiveFieldRef();

    for (const label facei : faces_)
    {
        auto& hval = hif[facei];

        if (hval > max_)
        {
            hval *= max_/max(hval, SMALL);
            ++nFacesAbove;
        }
    }

    // Handle boundaries in the case of 'all'
    label nEdgesAbove = 0;

    if (!faceSetOption::useSubMesh())
    {
        for (faPatchScalarField& hp : h.boundaryFieldRef())
        {
            if (!hp.fixesValue())
            {
                for (auto& hval : hp)
                {
                    if (hval > max_)
                    {
                        hval *= max_/max(hval, SMALL);
                        ++nEdgesAbove;
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


    reduce(nFacesAbove, sumOp<label>());
    reduce(nEdgesAbove, sumOp<label>());

    Info<< type() << ' ' << name_ << " Limited "
        << nFacesAbove << " ("
        << percent(nFacesAbove, nTotFaces)
        << "%) of faces, with max limit " << max_ << endl;

    if (nFacesAbove || nEdgesAbove)
    {
        // We've changed internal values so give
        // boundary conditions opportunity to correct
        h.correctBoundaryConditions();
    }
}


// ************************************************************************* //
