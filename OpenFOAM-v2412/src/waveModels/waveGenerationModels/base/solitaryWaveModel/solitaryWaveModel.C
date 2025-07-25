/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 IH-Cantabria
    Copyright (C) 2016-2017 OpenCFD Ltd.
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

#include "solitaryWaveModel.H"
#include "polyPatch.H"
#include "SubField.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace waveModels
{
    defineTypeNameAndDebug(solitaryWaveModel, 0);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::waveModels::solitaryWaveModel::timeCoeff
(
    const scalar t
) const
{
    // Ramping not applicable to solitary waves
    return 1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveModels::solitaryWaveModel::solitaryWaveModel
(
    const dictionary& dict,
    const fvMesh& mesh,
    const polyPatch& patch,
    const bool readFields
)
:
    waveGenerationModel(dict, mesh, patch, false),
    waveHeight_(0),
    waveAngle_(0),
    x_
    (
        patch.faceCentres().component(0)*cos(waveAngle_)
      + patch.faceCentres().component(1)*sin(waveAngle_)
    ),
    x0_(gMin(x_))
{
    if (readFields)
    {
        readDict(dict);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::waveModels::solitaryWaveModel::readDict
(
    const dictionary& overrideDict
)
{
    if (waveGenerationModel::readDict(overrideDict))
    {
        waveHeight_ = readWaveHeight();
        waveAngle_ = readWaveAngle();

        return true;
    }

    return false;
}


void Foam::waveModels::solitaryWaveModel::info(Ostream& os) const
{
    waveGenerationModel::info(os);

    os  << "    Wave height : " << waveHeight_ << nl
        << "    Wave angle  : " << radToDeg(waveAngle_) << nl
        << "    x0: " << x0_ << nl;
}


// ************************************************************************* //
