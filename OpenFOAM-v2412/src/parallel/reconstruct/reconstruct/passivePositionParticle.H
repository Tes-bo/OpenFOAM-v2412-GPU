/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2024 OpenCFD Ltd.
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
    Foam::passivePositionParticle

Description
    Copy of base particle but without particle locating and preserving
    read location.

    Used in reconstructing lagrangian positions generated outside the
    mesh domain (can happen in extractEulerianParticles functionObject)

SourceFiles
    passivePositionParticle.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_passivePositionParticle_H
#define Foam_passivePositionParticle_H

#include "particle.H"
#include "IOstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class passivePositionParticle Declaration
\*---------------------------------------------------------------------------*/

class passivePositionParticle
:
    public particle
{
    // Private Data

        //- Raw location
        point location_;


public:

    // Constructors

        //- Construct from components (known location)
        passivePositionParticle
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti
        )
        :
            particle(mesh, coordinates, celli, tetFacei, tetPti),
            location_(position())
        {}

        //- Construct from components (supplied location)
        passivePositionParticle
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti,
            const vector& position
        )
        :
            particle(mesh, coordinates, celli, tetFacei, tetPti),
            location_(position)
        {}

        //- Construct from Istream
        passivePositionParticle
        (
            const polyMesh& mesh,
            Istream& is,
            const bool readFields = true,
            const bool newFormat = true
        )
        :
            //particle(mesh, is, readFields, newFormat)
            particle
            (
                mesh,
                Zero,   // position
                -1,     // celli
                -1,     // tetFacei
                -1,     // tetPti
                false   // doLocate
            )
        {
            readData
            (
                is,
                location_,
                readFields,
                newFormat,
                false       //doLocate
            );
        }


        //- Construct as copy
        passivePositionParticle(const passivePositionParticle& p)
        :
            particle(p),
            location_(p.location_)
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

            autoPtr<passivePositionParticle> operator()(Istream& is) const
            {
                return autoPtr<passivePositionParticle>::New(mesh_, is, true);
            }
        };


    // Member Functions

        //- The cached particle position
        const point& location() const noexcept
        {
            return location_;
        }

        //- Write the particle position and cell id
        //  Uses cached location() instead of calculated position()
        virtual void writePosition(Ostream& os) const
        {
            if (os.format() == IOstreamOption::ASCII)
            {
                os  << location() << token::SPACE << cell();
            }
            else
            {
                positionsCompat1706 p;

                const size_t s =
                (
                    offsetof(positionsCompat1706, facei)
                  - offsetof(positionsCompat1706, position)
                );

                p.position = location();
                p.celli = cell();

                os.write(reinterpret_cast<const char*>(&p.position), s);
            }

            // Check state of Ostream
            os.check(FUNCTION_NAME);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
