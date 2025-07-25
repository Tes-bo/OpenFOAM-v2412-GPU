/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2017-2019 OpenCFD Ltd.
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

#include "TRIReader.H"
#include "surfaceFormatsCore.H"
#include "IFstream.H"
#include "IOmanip.H"
#include "StringStream.H"
#include "DynamicList.H"
#include "mergePoints.H"
#include "Map.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
static inline STLpoint getSTLpoint(Istream& is)
{
    scalar a = readScalar(is);
    scalar b = readScalar(is);
    scalar c = readScalar(is);

    return STLpoint(a, b, c);
}
} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::fileFormats::TRIReader::readFile(const fileName& filename)
{
    // Clear everything
    this->clear();
    sorted_ = true;

    IFstream is(filename);
    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot read file " << filename << nl
            << exit(FatalError);
    }

    // uses similar structure as STL, just some points
    // the rest of the reader resembles the STL binary reader
    DynamicList<STLpoint> dynPoints;
    DynamicList<label> dynZones;
    DynamicList<word>  dynNames;
    DynamicList<label> dynSizes;
    HashTable<label>   lookup;

    // place faces without a group in zone0
    label zoneI = 0;
    dynSizes.append(zoneI);
    lookup.insert("zoneI", zoneI);

    while (is.good())
    {
        string line = this->getLineNoComment(is);

        if (line.empty())
        {
            break;
        }

        // Do not handle continuations

        IStringStream lineStream(line);

        STLpoint p(getSTLpoint(lineStream));

        if (!lineStream) break;

        dynPoints.append(p);
        dynPoints.append(getSTLpoint(lineStream));
        dynPoints.append(getSTLpoint(lineStream));

        // zone/colour in .tri file starts with 0x. Skip.
        // ie, instead of having 0xFF, skip 0 and leave xFF to
        // get read as a word and name it "zoneFF"

        char zeroChar;
        lineStream >> zeroChar;

        const word rawName(lineStream);
        const word name("zone" + rawName.substr(1));

        const auto iter = lookup.cfind(name);
        if (iter.good())
        {
            if (zoneI != iter.val())
            {
                sorted_ = false; // Group appeared out of order
                zoneI = iter.val();
            }
        }
        else
        {
            zoneI = dynSizes.size();
            if (lookup.insert(name, zoneI))
            {
                dynNames.append(name);
                dynSizes.append(0);
            }
        }

        dynZones.append(zoneI);
        dynSizes[zoneI]++;
    }

    // skip empty groups
    label nZone = 0;
    forAll(dynSizes, zonei)
    {
        if (dynSizes[zonei])
        {
            if (nZone != zonei)
            {
                dynNames[nZone] = dynNames[zonei];
                dynSizes[nZone] = dynSizes[zonei];
            }
            ++nZone;
        }
    }

    // truncate addressed size
    dynNames.setCapacity(nZone);
    dynSizes.setCapacity(nZone);

    // transfer to normal lists
    points_.transfer(dynPoints);
    zoneIds_.transfer(dynZones);
    names_.transfer(dynNames);
    sizes_.transfer(dynSizes);

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileFormats::TRIReader::TRIReader
(
    const fileName& filename
)
:
    sorted_(true),
    points_(),
    zoneIds_(),
    names_(),
    sizes_()
{
    readFile(filename);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fileFormats::TRIReader::clear()
{
    sorted_ = true;
    points_.clear();
    zoneIds_.clear();
    names_.clear();
    sizes_.clear();
}


Foam::label Foam::fileFormats::TRIReader::mergePointsMap
(
    labelList& pointMap
) const
{
    // Use merge tolerance as per STL ASCII
    return mergePointsMap
    (
        100 * doubleScalarSMALL,
        pointMap
    );
}


Foam::label Foam::fileFormats::TRIReader::mergePointsMap
(
    const scalar mergeTol,
    labelList& pointMap
) const
{
    return Foam::mergePoints
    (
        points_,
        mergeTol,
        false, // verbose
        pointMap
    );
}


// ************************************************************************* //
