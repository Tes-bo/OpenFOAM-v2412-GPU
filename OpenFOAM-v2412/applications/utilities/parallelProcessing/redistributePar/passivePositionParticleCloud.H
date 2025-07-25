/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

Class
    Foam::passivePositionParticleCloud

Description
    A Cloud of passive position particles

SourceFiles
    passivePositionParticleCloud.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_passivePositionParticleCloud_H
#define Foam_passivePositionParticleCloud_H

#include "Cloud.H"
#include "passivePositionParticle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class passivePositionParticleCloud Declaration
\*---------------------------------------------------------------------------*/

class passivePositionParticleCloud
:
    public Cloud<passivePositionParticle>
{
public:

    //- Type of parcel within the cloud
    typedef passivePositionParticle parcelType;

    //- No copy construct
    passivePositionParticleCloud(const passivePositionParticleCloud&) = delete;

    //- No copy assignment
    void operator=(const passivePositionParticleCloud&) = delete;


    // Constructors

        //- Read construct
        explicit passivePositionParticleCloud
        (
            const polyMesh& mesh,
            const word& cloudName = cloud::defaultName,
            bool readFields = true
        );

        //- Construct without particles
        passivePositionParticleCloud
        (
            const polyMesh& mesh,
            const Foam::zero,
            const word& cloudName = cloud::defaultName
        )
        :
            Cloud<passivePositionParticle>(mesh, Foam::zero{}, cloudName)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
