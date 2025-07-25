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

#include "waveGenerationModel.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace waveModels
{
    defineTypeNameAndDebug(waveGenerationModel, 0);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::waveModels::waveGenerationModel::readWaveHeight() const
{
    const scalar h(get<scalar>("waveHeight"));
    if (h < 0)
    {
        FatalIOErrorInFunction(*this)
            << "Wave height must be greater than zero.  Supplied"
            << " value waveHeight = " << h
            << exit(FatalIOError);
    }

    return h;
}


Foam::scalar Foam::waveModels::waveGenerationModel::readWaveAngle() const
{
    return degToRad(get<scalar>("waveAngle"));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveModels::waveGenerationModel::waveGenerationModel
(
    const dictionary& dict,
    const fvMesh& mesh,
    const polyPatch& patch,
    const bool readFields
)
:
    waveModel(dict, mesh, patch, false)
{
    if (readFields)
    {
        readDict(dict);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::waveModels::waveGenerationModel::readDict
(
    const dictionary& overrideDict
)
{
    if (waveModel::readDict(overrideDict))
    {
        readEntry("activeAbsorption", activeAbsorption_);

        return true;
    }

    return false;
}


void Foam::waveModels::waveGenerationModel::info(Ostream& os) const
{
    waveModel::info(os);
}


// ************************************************************************* //
