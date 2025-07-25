/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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

#include "DMD.H"
#include "readFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace ROMmodels
{
    defineTypeNameAndDebug(DMD, 0);
    addToRunTimeSelectionTable(ROMmodel, DMD, dictionary);
}
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Implementation
#include "DMDImpl.C"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::wordList Foam::ROMmodels::DMD::modeNames(const word& modeType) const
{
    wordList modeNames(modes_.size());
    for (const label modei : modes_)
    {
        modeNames[modei] =
            word
            (
                "mode"+modeType+"_"+name(modei)+"_"+fieldName_+"_"+objectName_
            );
    }
    return modeNames;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ROMmodels::DMD::DMD
(
    Time& runTime,
    fvMesh& mesh,
    const dictionary& dict,
    const instantList& times
)
:
    ROMmodel(runTime, mesh, dict, times),
    fieldName_(),
    objectName_(),
    deltaT_(),
    time_(),
    startTime_(),
    modes_(),
    dims_(),
    amps_(),
    evals_()
{
    read(dict);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::ROMmodels::DMD::read(const dictionary& dict)
{
    dict.readEntry("field", fieldName_);
    dict.readEntry("object", objectName_);
    dict.readEntry("deltaT", deltaT_);
    dict.readEntry("time", time_);
    dict.readIfPresent("startTime", startTime_);

    if (deltaT_ < SMALL || time_ < SMALL || startTime_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "Out-of-range values for " << nl
            << tab << "deltaT: " << deltaT_
            << tab << "time: " << time_ << nl
            << tab << "startTime: " << startTime_ << nl
            << exit(FatalIOError);
    }

    dict.readEntry("modes", modes_);

    if (modes_.empty())
    {
        FatalIOErrorInFunction(dict)
            << "Empty list for the mode indices " << nl
            << tab << "modes: " << modes_ << nl
            << exit(FatalIOError);
    }

    dims_.reset(dict.get<dimensionSet>("dimensions"));

    // Load complex amplitudes and eigenvalues
    dict.readEntry("amplitudes", amps_);
    dict.readEntry("eigenvalues", evals_);

    if ((amps_.size() != modes_.size()) || (evals_.size() != modes_.size()))
    {
        FatalIOErrorInFunction(dict)
            << "Inconsistent input sizes for "
            << tab << "modes: " << modes_.size() << nl
            << tab << "amplitudes: " << amps_.size() << nl
            << tab << "eigenvalues: " << evals_.size() << nl
            << exit(FatalIOError);
    }

    // Create mode-field names
    const wordList modeReNames(modeNames(word("Re")));
    const wordList modeImNames(modeNames(word("Im")));

    // Load mode fields
    runTime_.setTime(time_, 0);

    readFieldsHandler(mesh_).execute(modeReNames);
    readFieldsHandler(mesh_).execute(modeImNames);

    return true;
}


bool Foam::ROMmodels::DMD::createAndWrite()
{
    do
    {
        #undef  doLocalCode
        #define doLocalCode(InputType)                                      \
        {                                                                   \
            createAndWriteImpl<VolumeField<InputType>>();                   \
            createAndWriteImpl<SurfaceField<InputType>>();                  \
        }

        doLocalCode(scalar);
        doLocalCode(vector);
        doLocalCode(sphericalTensor);
        doLocalCode(symmTensor);
        doLocalCode(tensor);

        #undef doInnerLocalCode
        #undef doLocalCode
    }
    while (false);

    return true;
}


// ************************************************************************* //
