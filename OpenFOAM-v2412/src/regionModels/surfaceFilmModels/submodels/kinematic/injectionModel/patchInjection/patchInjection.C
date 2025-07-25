/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2017 OpenFOAM Foundation
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

#include "patchInjection.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(patchInjection, 0);
addToRunTimeSelectionTable(injectionModel, patchInjection, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

patchInjection::patchInjection
(
    surfaceFilmRegionModel& film,
    const dictionary& dict
)
:
    injectionModel(type(), film, dict),
    deltaStable_(coeffDict_.getOrDefault<scalar>("deltaStable", 0))
{
    const polyBoundaryMesh& pbm = film.regionMesh().boundaryMesh();

    const label nPatches
    (
        pbm.size()
      - film.regionMesh().globalData().processorPatches().size()
    );

    wordRes patchNames;
    if (coeffDict_.readIfPresent("patches", patchNames))
    {
        // Can also use pbm.indices(), but no warnings...
        patchIDs_ = pbm.patchSet(patchNames).sortedToc();

        Info<< "        applying to patches:" << nl;

        for (const label patchi : patchIDs_)
        {
            Info<< "            " << pbm[patchi].name() << endl;
        }
    }
    else
    {
        Info<< "            applying to all patches" << endl;

        patchIDs_ = identity(nPatches);
    }

    patchInjectedMasses_.resize(patchIDs_.size(), Zero);

    if (patchIDs_.empty())
    {
        FatalErrorInFunction
            << "No patches selected"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void patchInjection::correct
(
    scalarField& availableMass,
    scalarField& massToInject,
    scalarField& diameterToInject
)
{
    // Do not correct if no patches selected
    if (patchIDs_.empty()) return;

    const scalarField& delta = film().delta();
    const scalarField& rho = film().rho();
    const scalarField& magSf = film().magSf();

    const polyBoundaryMesh& pbm = film().regionMesh().boundaryMesh();

    forAll(patchIDs_, pidi)
    {
        label patchi = patchIDs_[pidi];
        const polyPatch& pp = pbm[patchi];
        const labelList& faceCells = pp.faceCells();

        // Accumulate the total mass removed from patch
        scalar dMassPatch = 0;

        for (const label celli : faceCells)
        {
            scalar ddelta = max(0.0, delta[celli] - deltaStable_);
            scalar dMass = ddelta*rho[celli]*magSf[celli];
            massToInject[celli] += dMass;
            availableMass[celli] -= dMass;
            dMassPatch += dMass;
        }

        patchInjectedMasses_[pidi] += dMassPatch;
        addToInjectedMass(dMassPatch);
    }

    injectionModel::correct();

    if (writeTime())
    {
        scalarField patchInjectedMasses0
        (
            getModelProperty<scalarField>
            (
                "patchInjectedMasses",
                scalarField(patchInjectedMasses_.size(), Zero)
            )
        );

        scalarField patchInjectedMassTotals(patchInjectedMasses_);
        Pstream::listCombineGather(patchInjectedMassTotals, plusEqOp<scalar>());
        patchInjectedMasses0 += patchInjectedMassTotals;

        setModelProperty<scalarField>
        (
            "patchInjectedMasses",
            patchInjectedMasses0
        );

        patchInjectedMasses_ = 0;
    }
}


void patchInjection::patchInjectedMassTotals(scalarField& patchMasses) const
{
    // Do not correct if no patches selected
    if (!patchIDs_.size()) return;

    scalarField patchInjectedMasses
    (
        getModelProperty<scalarField>
        (
            "patchInjectedMasses",
            scalarField(patchInjectedMasses_.size(), Zero)
        )
    );

    scalarField patchInjectedMassTotals(patchInjectedMasses_);
    Pstream::listCombineGather(patchInjectedMassTotals, plusEqOp<scalar>());

    forAll(patchIDs_, pidi)
    {
        label patchi = patchIDs_[pidi];
        patchMasses[patchi] +=
            patchInjectedMasses[pidi] + patchInjectedMassTotals[pidi];
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
