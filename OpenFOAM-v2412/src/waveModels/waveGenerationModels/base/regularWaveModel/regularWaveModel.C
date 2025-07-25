/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 IH-Cantabria
    Copyright (C) 2016-2018 OpenCFD Ltd.
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

#include "regularWaveModel.H"
#include "unitConversion.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace waveModels
{
    defineTypeNameAndDebug(regularWaveModel, 0);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::word Foam::waveModels::regularWaveModel::waveType() const
{
    scalar waveK = 2.0*mathematical::pi/waveLength_;

    word generation = "Intermediate";
    if (waveK*waterDepthRef_ > mathematical::pi)
    {
        generation = "Deep";
    }
    else if (waveK*waterDepthRef_ < 0.1*mathematical::pi)
    {
        generation = "Shallow";
    }

    return generation;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveModels::regularWaveModel::regularWaveModel
(
    const dictionary& dict,
    const fvMesh& mesh,
    const polyPatch& patch,
    const bool readFields
)
:
    irregularWaveModel(dict, mesh, patch, false),
    waveHeight_(0),
    waveAngle_(0),
    wavePeriod_(0),
    waveLength_(0),
    wavePhase_(1.5*mathematical::pi)
{
    if (readFields)
    {
        readDict(dict);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::waveModels::regularWaveModel::readDict
(
    const dictionary& overrideDict
)
{
    if (irregularWaveModel::readDict(overrideDict))
    {
        waveHeight_ = readWaveHeight();
        waveAngle_ = readWaveAngle();

        readEntry("wavePeriod", wavePeriod_);

        if (wavePeriod_ < 0)
        {
            FatalIOErrorInFunction(*this)
                << "Wave period must be greater than zero.  Supplied"
                << " value wavePeriod = " << wavePeriod_
                << exit(FatalIOError);
        }

        readIfPresent("wavePhase", wavePhase_);

        // Note: waveLength to be set in derived classes

        return true;
    }

    return false;
}


void Foam::waveModels::regularWaveModel::info(Ostream& os) const
{
    irregularWaveModel::info(os);

    os  << "    Wave height : " << waveHeight_ << nl
        << "    Wave angle  : " << radToDeg(waveAngle_) << nl
        << "    Wave period : " << wavePeriod_ << nl
        << "    Wave length : " << waveLength_ << nl
        << "    Wave phase : " << wavePhase_ << nl;
}


// ************************************************************************* //
