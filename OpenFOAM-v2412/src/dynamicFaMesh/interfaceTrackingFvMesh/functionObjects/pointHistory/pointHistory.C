/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 Zeljko Tukovic, FSB Zagreb.
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

\*----------------------------------------------------------------------------*/

#include "pointHistory.H"
#include "addToRunTimeSelectionTable.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pointHistory, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        pointHistory,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::pointHistory::writeData()
{
    const fvMesh& mesh = time_.lookupObject<fvMesh>(polyMesh::defaultRegion);

    vector position(Zero);

    if (processor_ == Pstream::myProcNo())
    {
        position = mesh.points()[historyPointID_];
    }

    reduce(position, sumOp<vector>());

    if (Pstream::master())
    {
        historyFilePtr_() << setprecision(12);

        historyFilePtr_()
            << time_.time().value() << tab
            << position.x() << tab
            << position.y() << tab
            << position.z() << endl;
    }

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointHistory::pointHistory
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(runTime),
    regionName_(polyMesh::defaultRegion),
    historyPointID_(-1),
    refHistoryPoint_(dict.lookup("refHistoryPoint")),
    processor_(-1),
    fileName_(dict.get<word>("fileName")),
    historyFilePtr_(nullptr)
{
    Info<< "Creating " << this->name() << " function object." << endl;

    dict.readIfPresent("region", regionName_);
    dict.readIfPresent("historyPointID", historyPointID_);

    const fvMesh& mesh = time_.lookupObject<fvMesh>(regionName_);

    const vectorField& points = mesh.points();

    List<scalar> minDist(Pstream::nProcs(), GREAT);

    if (historyPointID_ == -1)
    {
        forAll(points, pointI)
        {
            scalar dist = mag(refHistoryPoint_ - points[pointI]);

            if (dist < minDist[Pstream::myProcNo()])
            {
                minDist[Pstream::myProcNo()] = dist;
                historyPointID_ = pointI;
            }
        }
    }

    Pstream::allGatherList(minDist);

    processor_ = -1;
    scalar min = GREAT;

    forAll(minDist, procI)
    {
        if (minDist[procI] < min)
        {
            min = minDist[procI];
            processor_ = procI;
        }
    }

    if (processor_ == Pstream::myProcNo())
    {
        Pout<< "History point ID: " << historyPointID_ << nl
            << "History point coordinates: "
            << points[historyPointID_] << nl
            << "Reference point coordinates: " << refHistoryPoint_
            << endl;
    }

    // Create history file if not already created
    if (!historyFilePtr_)
    {
        // File update
        if (Pstream::master())
        {
            const word startTimeName =
                time_.timeName(mesh.time().startTime().value());

            const fileName historyDir =
                time_.globalPath()/"history"/startTimeName;

            // Create directory if does not exist.
            mkDir(historyDir);


            // Open new file at start up

            // OStringStream FileName;
            // FileName() << "point_" << historyPointID_ << ".dat";

            historyFilePtr_.reset
            (
                new OFstream(historyDir/fileName_)
            );

            // Add headers to output data
            if (historyFilePtr_)
            {
                historyFilePtr_()
                    << "# Time" << tab << "X" << tab << "Y" << tab << "Z"
                    << endl;
            }
        }
    }

    // Write start time data
    writeData();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::pointHistory::execute()
{
    return writeData();
}


bool Foam::pointHistory::read(const dictionary& dict)
{
    dict.readIfPresent("region", regionName_);

    return true;
}


// ************************************************************************* //
