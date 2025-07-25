/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

#include "volRegion.H"
#include "volMesh.H"
#include "cellSet.H"
#include "globalMeshData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(volRegion, 0);
}
}


const Foam::Enum
<
    Foam::functionObjects::volRegion::regionTypes
>
Foam::functionObjects::volRegion::regionTypeNames_
({
    { regionTypes::vrtAll, "all" },
    { regionTypes::vrtCellSet, "cellSet" },
    { regionTypes::vrtCellZone, "cellZone" },
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::volRegion::calculateCache()
{
    cellIds_.clear();
    regionIDs_.clear();

    // Update now. Need a valid state for the cellIDs() call
    requireUpdate_ = false;

    switch (regionType_)
    {
        case vrtAll:
        {
            nCells_ = volMesh_.globalData().nTotalCells();
            V_ = gSum(volMesh_.V());
            return;
            break;
        }

        case vrtCellSet:
        {
            cellIds_ = cellSet(volMesh_, regionName_).sortedToc();
            break;
        }

        case vrtCellZone:
        {
            regionIDs_ = volMesh_.cellZones().indices(regionName_);

            if (regionIDs_.empty())
            {
                FatalErrorInFunction
                    << "Unknown cell zone: " << regionName_ << nl
                    << "    Valid zones : "
                    << flatOutput(volMesh_.cellZones().names()) << nl
                    << "    Valid groups: "
                    << flatOutput(volMesh_.cellZones().groupNames()) << nl
                    << exit(FatalError);
            }

            if (regionIDs_.size() > 1)
            {
                cellIds_ =
                    volMesh_.cellZones().selection(regionIDs_).sortedToc();
            }
            break;
        }
    }


    // Calculate cache value for nCells() and V()
    const labelList& selected = this->cellIDs();

    V_ = 0;
    for (const label celli : selected)
    {
        V_ += volMesh_.V()[celli];
    }

    nCells_ = returnReduce(selected.size(), sumOp<label>());
    reduce(V_, sumOp<scalar>());

    if (!nCells_)
    {
        FatalErrorInFunction
            << regionTypeNames_[regionType_]
            << '(' << regionName_ << "):" << nl
            << "    Region has no cells" << nl
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::volRegion::writeFileHeader
(
    const writeFile& wf,
    Ostream& file
) const
{
    wf.writeCommented(file, "Region");
    file<< setw(1) << ':' << setw(1) << ' '
        << regionTypeNames_[regionType_] << ' ' << regionName_ << endl;
    wf.writeHeaderValue(file, "Cells", nCells_);
    wf.writeHeaderValue(file, "Volume", V_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::volRegion::volRegion
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    volMesh_(mesh),
    cellIds_(),
    regionIDs_(),
    nCells_(0),
    V_(Zero),
    requireUpdate_(true),
    regionType_
    (
        regionTypeNames_.getOrDefault
        (
            "regionType",
            dict,
            regionTypes::vrtAll
        )
    ),
    regionName_(volMesh_.name())
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::volRegion::read(const dictionary& dict)
{
    switch (regionType_)
    {
        case vrtAll:
        {
            regionName_ = volMesh_.name();
            break;
        }

        case vrtCellSet:
        case vrtCellZone:
        {
            dict.readEntry("name", regionName_);
            break;
        }

        default:
        {
            FatalIOErrorInFunction(dict)
                << "Unknown region type. Valid region types: "
                << flatOutput(regionTypeNames_.names()) << nl
                << exit(FatalIOError);
            break;
        }
    }

    calculateCache();
    return true;
}


const Foam::labelList& Foam::functionObjects::volRegion::cellIDs() const
{
    #ifdef FULLDEBUG
    if (requireUpdate_)
    {
        FatalErrorInFunction
            << "Retrieving cached values that are not up-to-date" << nl
            << exit(FatalError);
    }
    #endif

    switch (regionType_)
    {
        case vrtCellSet:
        {
            return cellIds_;
            break;
        }

        case vrtCellZone:
        {
            if (regionIDs_.size() == 1)
            {
                return volMesh_.cellZones()[regionIDs_.first()];
            }
            else
            {
                return cellIds_;
            }
            break;
        }

        default:
            break;
    }

    return labelList::null();
}


bool Foam::functionObjects::volRegion::update()
{
    if (requireUpdate_)
    {
        calculateCache();
        return true;
    }

    return false;
}


void Foam::functionObjects::volRegion::updateMesh(const mapPolyMesh&)
{
    requireUpdate_ = true;
}


void Foam::functionObjects::volRegion::movePoints(const polyMesh&)
{
    requireUpdate_ = true;
}


// ************************************************************************* //
