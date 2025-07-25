/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
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
    Foam::streamLineParticleCloud

Description
    A Cloud of streamLine particles

SourceFiles
    streamLineParticleCloud.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_streamLineParticleCloud_H
#define Foam_streamLineParticleCloud_H

#include "Cloud.H"
#include "streamLineParticle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class streamLineParticleCloud Declaration
\*---------------------------------------------------------------------------*/

class streamLineParticleCloud
:
    public Cloud<streamLineParticle>
{
public:

    //- Type of parcel within the cloud
    typedef streamLineParticle parcelType;

    //- No copy construct
    streamLineParticleCloud(const streamLineParticleCloud&) = delete;

    //- No copy assignment
    void operator=(const streamLineParticleCloud&) = delete;


    // Constructors

        //- Read construct
        explicit streamLineParticleCloud
        (
            const polyMesh& mesh,
            const word& cloudName = cloud::defaultName,
            bool readFields = true
        );

        //- Construct without particles
        streamLineParticleCloud
        (
            const polyMesh& mesh,
            const Foam::zero,
            const word& cloudName = cloud::defaultName
        )
        :
            Cloud<streamLineParticle>(mesh, Foam::zero{}, cloudName)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
