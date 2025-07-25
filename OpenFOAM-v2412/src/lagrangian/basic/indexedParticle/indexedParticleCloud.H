/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::indexedParticleCloud

Description
    A Cloud of particles carrying an additional index.

SourceFiles
    indexedParticleCloud.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_indexedParticleCloud_H
#define Foam_indexedParticleCloud_H

#include "Cloud.H"
#include "indexedParticle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class indexedParticleCloud Declaration
\*---------------------------------------------------------------------------*/

class indexedParticleCloud
:
    public Cloud<indexedParticle>
{
public:

    //- Type of parcel within the cloud
    typedef indexedParticle parcelType;

    //- No copy construct
    indexedParticleCloud(const indexedParticleCloud&) = delete;

    //- No copy assignment
    void operator=(const indexedParticleCloud&) = delete;


    // Constructors

        //- Read construct
        explicit indexedParticleCloud
        (
            const polyMesh& mesh,
            const word& cloudName = cloud::defaultName,
            bool readFields = true
        );

        //- Construct without particles
        explicit indexedParticleCloud
        (
            const polyMesh& mesh,
            const Foam::zero,
            const word& cloudName = cloud::defaultName
        )
        :
            Cloud<indexedParticle>(mesh, Foam::zero{}, cloudName)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
