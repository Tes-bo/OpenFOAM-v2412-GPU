/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2023 OpenCFD Ltd.
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

#include "singleDirectionUniformBin.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace binModels
{
    defineTypeNameAndDebug(singleDirectionUniformBin, 0);
    addToRunTimeSelectionTable(binModel, singleDirectionUniformBin, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::binModels::singleDirectionUniformBin::initialise()
{
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    // Use geometry limits if not specified by the user
    const bool useGeomLimits
    (
        binLimits_.min() == GREAT
     || binLimits_.max() == GREAT
    );

    if (useGeomLimits)
    {
        // Determine extents of patches/cells in a given direction
        scalarMinMax geomLimits;

        for (const label patchi : patchIDs_)
        {
            for (const vector& p : pbm[patchi].faceCentres())
            {
                geomLimits.add(p & binDir_);
            }
        }

        for (const label zonei : cellZoneIDs_)
        {
            for (const label celli : mesh_.cellZones()[zonei])
            {
                geomLimits.add(mesh_.C()[celli] & binDir_);
            }
        }

        // Globally consistent
        reduce(geomLimits, minMaxOp<scalar>());

        if (!geomLimits.good())
        {
            FatalErrorInFunction
                << "No patches/cellZones provided"
                << exit(FatalError);
        }

        // Slightly boost max so that region of interest is fully within bounds
        // TBD: also adjust min?
        const scalar adjust(1e-4*geomLimits.span());
        geomLimits.max() += adjust;

        // Use geometry limits if not specified by the user
        if (binLimits_.min() == GREAT)
        {
            binLimits_.min() = geomLimits.min();
        }
        if (binLimits_.max() == GREAT)
        {
            binLimits_.max() = geomLimits.max();
        }
    }

    binWidth_ = binLimits_.span()/scalar(nBin_);

    if (binWidth_ <= 0)
    {
        FatalErrorInFunction
            << "Max bound must be greater than min bound" << nl
            << "    d           = " << binWidth_ << nl
            << "    min         = " << binLimits_.min() << nl
            << "    max         = " << binLimits_.max() << nl
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::binModels::singleDirectionUniformBin::singleDirectionUniformBin
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& outputPrefix
)
:
    binModel(dict, mesh, outputPrefix),
    binWidth_(0),
    binLimits_(GREAT),
    binDir_(Zero)
{
    read(dict);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::binModels::singleDirectionUniformBin::read(const dictionary& dict)
{
    if (!binModel::read(dict))
    {
        return false;
    }

    Info<< "    Activating a set of single-direction bins" << endl;

    const dictionary& binDict = dict.subDict("binData");

    nBin_ = binDict.getCheck<label>("nBin", labelMinMax::ge(1));

    Info<< "    Employing " << nBin_ << " bins" << nl;

    if (binDict.readIfPresent("min", binLimits_.min()))
    {
        Info<< "    - min        : " << binLimits_.min() << nl;
    }
    if (binDict.readIfPresent("max", binLimits_.max()))
    {
        Info<< "    - max        : " << binLimits_.max() << nl;
    }

    cumulative_ = binDict.getOrDefault<bool>("cumulative", false);
    Info<< "    - cumulative    : " << cumulative_ << nl
        << "    - decomposePatchValues    : " << decomposePatchValues_ << nl;

    binDir_ = binDict.get<vector>("direction");
    if (binDir_.mag() < SMALL)
    {
        FatalIOErrorInFunction(dict)
            << "Input direction should not be zero valued" << nl
            << "    direction = " << binDir_ << nl
            << exit(FatalIOError);
    }
    binDir_.normalise();

    Info<< "    - direction     : " << binDir_ << nl << endl;

    initialise();

    return true;
}


void Foam::binModels::singleDirectionUniformBin::apply()
{
    forAll(fieldNames_, i)
    {
        const bool ok =
        (
            processField<scalar>(i)
         || processField<vector>(i)
         || processField<sphericalTensor>(i)
         || processField<symmTensor>(i)
         || processField<tensor>(i)
        );

        if (!ok)
        {
            WarningInFunction
                << "Unable to find field " << fieldNames_[i]
                << ". Avaliable objects are "
                << mesh_.objectRegistry::sortedToc()
                << endl;
        }
    }

    writtenHeader_ = true;
}


// ************************************************************************* //
