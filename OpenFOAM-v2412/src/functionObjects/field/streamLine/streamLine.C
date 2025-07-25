/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

#include "streamLine.H"
#include "streamLineParticleCloud.H"
#include "sampledSet.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(streamLine, 0);
    addToRunTimeSelectionTable(functionObject, streamLine, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::streamLine::track()
{
    // Start with empty cloud
    streamLineParticleCloud particles(mesh_, Foam::zero{}, cloudName_);

    const sampledSet& seedPoints = sampledSetPoints();

    forAll(seedPoints, seedi)
    {
        particles.addParticle
        (
            new streamLineParticle
            (
                mesh_,
                seedPoints[seedi],
                seedPoints.cells()[seedi],
                (trackDir_ == trackDirType::FORWARD),
                lifeTime_
            )
        );

        if (trackDir_ == trackDirType::BIDIRECTIONAL)
        {
            // Add additional particle for the forward bit of the track
            particles.addParticle
            (
                new streamLineParticle
                (
                    mesh_,
                    seedPoints[seedi],
                    seedPoints.cells()[seedi],
                    true,
                    lifeTime_
                )
            );
        }
    }

    label nSeeds = returnReduce(particles.size(), sumOp<label>());

    Log << "    seeded " << nSeeds << " particles" << endl;

    // Field interpolators
    // Velocity interpolator
    PtrList<interpolation<scalar>> vsInterp;
    PtrList<interpolation<vector>> vvInterp;

    refPtr<interpolation<vector>> UInterp
    (
        initInterpolations(nSeeds, vsInterp, vvInterp)
    );

    // Additional particle info
    streamLineParticle::trackingData td
    (
        particles,
        vsInterp,
        vvInterp,
        UInterp.cref(), // velocity interpolator (possibly within vvInterp)
        nSubCycle_,     // automatic track control:step through cells in steps?
        trackLength_,   // fixed track length

        allTracks_,
        allScalars_,
        allVectors_
    );


    // Set very large dt. Note: cannot use GREAT since 1/GREAT is SMALL
    // which is a trigger value for the tracking...
    const scalar trackTime = Foam::sqrt(GREAT);

    // Track
    particles.move(particles, td, trackTime);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::streamLine::streamLine
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    streamLineBase(name, runTime, dict)
{
    read(dict_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::streamLine::read(const dictionary& dict)
{
    if (streamLineBase::read(dict))
    {
        bool subCycling = dict.found("nSubCycle");
        bool fixedLength = dict.found("trackLength");

        if (subCycling && fixedLength)
        {
            FatalIOErrorInFunction(dict)
                << "Cannot both specify automatic time stepping (through '"
                << "nSubCycle' specification) and fixed track length (through '"
                << "trackLength')"
                << exit(FatalIOError);
        }

        nSubCycle_ = 1;
        if (dict.readIfPresent("nSubCycle", nSubCycle_))
        {
            trackLength_ = VGREAT;
            nSubCycle_ = max(nSubCycle_, 1);

            Info<< "    automatic track length specified through"
                << " number of sub cycles : " << nSubCycle_ << nl
                << endl;
        }
    }
    return true;
}


// ************************************************************************* //
