/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2024 OpenCFD Ltd.
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
    Foam::passiveParticle

Description
    Copy of base particle

SourceFiles
    passiveParticle.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_passiveParticle_H
#define Foam_passiveParticle_H

#include "particle.H"
#include "IOstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class passiveParticle Declaration
\*---------------------------------------------------------------------------*/

class passiveParticle
:
    public particle
{
public:

    // Constructors

        //- Construct from components
        passiveParticle
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti
        )
        :
            particle(mesh, coordinates, celli, tetFacei, tetPti)
        {}


        //- Construct from a position and a cell.
        //  Searches for the rest of the required topology.
        passiveParticle
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli = -1
        )
        :
            particle(mesh, position, celli)
        {}


        //- Construct from Istream
        passiveParticle
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true,
            bool newFormat = true
        )
        :
            particle(mesh, is, readFields, newFormat)
        {}


        //- Construct as copy
        passiveParticle(const passiveParticle& p)
        :
            particle(p)
        {}


        //- Return a clone
        virtual autoPtr<particle> clone() const
        {
            return particle::Clone(*this);
        }


        //- Factory class to read-construct particles (for parallel transfer)
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<passiveParticle> operator()(Istream& is) const
            {
                return autoPtr<passiveParticle>::New(mesh_, is, true);
            }
        };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
