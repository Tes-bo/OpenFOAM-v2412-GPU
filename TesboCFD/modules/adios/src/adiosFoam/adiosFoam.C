/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "adiosFoam.H"
#include "fileOperation.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::word Foam::adiosFoam::defaultDirectory("adiosData");
const Foam::string Foam::adiosFoam::foamAttribute("/openfoam");
const Foam::string Foam::adiosFoam::timeAttribute("/time");


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

// adios-2.4.0 and earlier:
// - <time>.bp files and <time>.bp.dir directories
//
// later:
// - <time> directories

namespace Foam
{

static Foam::instantList findTimes_bp(const fileName& directory)
{
    // Read directory entries into a list
    const fileNameList dirEntries(Foam::readDir(directory, fileName::FILE, false));

    label nTimes = 0;

    // Check for .bp files
    for (const fileName& dirEntry : dirEntries)
    {
        if (dirEntry.hasExt("bp"))
        {
            ++nTimes;
        }
    }

    if (!nTimes)
    {
        return instantList();
    }

    instantList times(nTimes);
    nTimes = 0;

    // Parse directory entry into time value/name
    for (const fileName& dirEntry : dirEntries)
    {
        if (dirEntry.hasExt("bp"))
        {
            const word dirValue(dirEntry.lessExt());

            if (readScalar(dirValue, times[nTimes].value()))
            {
                times[nTimes].name() = word(dirEntry);
                ++nTimes;
            }
        }
    }

    times.resize(nTimes);

    if (nTimes > 1)
    {
        std::sort(times.begin(), times.end(), instant::less());
    }

    return times;
}

} // End namespace Foam


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

Foam::instantList Foam::adiosFoam::findTimes
(
    const fileName& directory,
    const word& constantName
)
{
    instantList times(findTimes_bp(directory));

    if (times.size())
    {
        // Older format with .bp files
        return times;
    }

    // Newer format, uses directory names only

    fileNameList dirEntries(Foam::readDir(directory, fileName::DIRECTORY));

    return fileOperation::sortTimes(dirEntries, constantName);
}


// ************************************************************************* //
