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

#include "unequalBinWidth.H"
#include "histogramModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace histogramModels
{
    defineTypeNameAndDebug(unequalBinWidth, 0);
    addToRunTimeSelectionTable(histogramModel, unequalBinWidth, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::histogramModels::unequalBinWidth::unequalBinWidth
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    histogramModel(name, mesh, dict),
    nBins_(-1),
    ranges_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::histogramModels::unequalBinWidth::read(const dictionary& dict)
{
    if (!histogramModel::read(dict))
    {
        return false;
    }

    ranges_ = dict.get<List<scalarMinMax>>("ranges");
    nBins_ = ranges_.size();

    forAll(ranges_, bini)
    {
        const auto& range = ranges_[bini];

        if (!range.good())
        {
            FatalIOErrorInFunction(dict)
                << "Histogram bin-" << bini
                << " has invalid range: " << range
                << abort(FatalIOError);
        }
    }

    if (nBins_ < 1)
    {
        FatalIOErrorInFunction(dict)
            << "Invalid number of histogram bins: " << nBins_
            << abort(FatalIOError);
    }

    return true;
}


bool Foam::histogramModels::unequalBinWidth::write(const bool log)
{
    // Retrieve operand field
    const volScalarField& field = histogramModel::getOrReadField(fieldName());

    // Calculate the mid-points of bins for the graph axis
    pointField midBin(nBins_, Zero);

    forAll(ranges_, bini)
    {
        midBin[bini].x() = ranges_[bini].centre();
    }


    // Calculate the histogram data
    scalarField dataNormalised(nBins_, Zero);
    labelField dataCount(nBins_, Zero);
    const scalarField& V = mesh().V();

    forAll(field, celli)
    {
        forAll(ranges_, bini)
        {
            const auto& range = ranges_[bini];

            // Like range.contains(field[celli]) but exclusive on max()

            if (field[celli] >= range.min() && field[celli] < range.max())
            {
                dataNormalised[bini] += V[celli];
                dataCount[bini]++;
                break;
            }
        }
    }
    Pstream::listCombineGather(dataNormalised, plusEqOp<scalar>());
    Pstream::listCombineGather(dataCount, plusEqOp<label>());


    // Write histogram data
    histogramModel::write(dataNormalised, dataCount, mag(midBin));

    return true;
}


// ************************************************************************* //
