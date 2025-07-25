/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2019 OpenCFD Ltd.
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
    Foam::injectedParticleCloud

Group

Description

SourceFiles
    injectedParticleCloud.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_injectedParticleCloud_H
#define Foam_injectedParticleCloud_H

#include "Cloud.H"
#include "injectedParticle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class injectedParticleCloud Declaration
\*---------------------------------------------------------------------------*/

class injectedParticleCloud
:
    public Cloud<injectedParticle>
{
public:

    //- Type of parcel within the cloud
    typedef injectedParticle parcelType;

    //- No copy construct
    injectedParticleCloud(const injectedParticleCloud&) = delete;

    //- No copy assignment
    void operator=(const injectedParticleCloud&) = delete;


    // Constructors

        //- Read construct
        injectedParticleCloud
        (
            const polyMesh& mesh,
            const word& cloudName,
            const bool readFields = true
        );

        //- Construct without particles
        injectedParticleCloud
        (
            const polyMesh& mesh,
            const Foam::zero,
            const word& cloudName = cloud::defaultName
        )
        :
            Cloud<injectedParticle>(mesh, Foam::zero{}, cloudName)
        {}

        //- Copy constructor with new name
        injectedParticleCloud(const injectedParticleCloud& c, const word& name);

        //- Construct and return clone based on (this) with new name
        virtual autoPtr<injectedParticleCloud> clone(const word& name) const
        {
            return autoPtr<injectedParticleCloud>
            (
                new injectedParticleCloud(*this, name)
            );
        }


    //- Destructor
    virtual ~injectedParticleCloud() = default;


    // Member Functions

        //- Read particle fields as objects from the obr registry
        virtual void readObjects(const objectRegistry& obr);

        //- Write particle fields as objects into the obr registry
        virtual void writeObjects(objectRegistry& obr) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
