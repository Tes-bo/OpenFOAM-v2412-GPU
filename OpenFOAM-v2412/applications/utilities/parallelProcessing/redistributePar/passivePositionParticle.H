/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2024 OpenCFD Ltd.
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
    Passive particle, transferring in old format (i.e. position instead of
    coordinates). Used for e.g. redistributePar.

SourceFiles
    passivePositionParticle.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_passivePositionParticle_H
#define Foam_passivePositionParticle_H

#include "passiveParticle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class passivePositionParticle;
Ostream& operator<<(Ostream& os, const passivePositionParticle& ppi);

/*---------------------------------------------------------------------------*\
                   Class passivePositionParticle Declaration
\*---------------------------------------------------------------------------*/

class passivePositionParticle
:
    public passiveParticle
{
    // Private Member Data

        //- Cached position
        point cachedPosition_;


public:

    // Constructors

        //- Construct from Istream in old format
        passivePositionParticle
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields,
            bool newFormat
        )
        :
            passiveParticle(mesh, is, readFields, newFormat),
            cachedPosition_(position())
        {}

        //- Construct from a position and a cell.
        //  Searches for the rest of the required topology.
        passivePositionParticle
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli = -1
        )
        :
            passiveParticle(mesh, position, celli),
            cachedPosition_(position)
        {}

        //- Construct as copy
        passivePositionParticle(const passivePositionParticle& p)
        :
            passiveParticle(p),
            cachedPosition_(p.cachedPosition_)
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
                return autoPtr<passivePositionParticle>
                (
                    // Read in old format
                    new passivePositionParticle(mesh_, is, true, false)
                );
            }
        };


    // Member Functions

        const point& cachedPosition() const noexcept
        {
            return cachedPosition_;
        }


    // Friend Operators

        friend Ostream& operator<<
        (
            Ostream& os,
            const passivePositionParticle& ppi
        )
        {
            // Copy data into old format structure. Exact opposite of
            // particleIO.C reading old format.
            particle::positionsCompat1706 p;

            p.position = ppi.cachedPosition_;
            p.celli = ppi.cell();
            p.facei = ppi.face();
            p.stepFraction = ppi.stepFraction();
            p.tetFacei = ppi.tetFace();
            p.tetPti = ppi.tetPt();
            p.origProc = ppi.origProc();
            p.origId = ppi.origId();

            if (os.format() == IOstreamOption::ASCII)
            {
                os  << p.position
                    << token::SPACE << p.celli
                    << token::SPACE << p.facei
                    << token::SPACE << p.stepFraction
                    << token::SPACE << p.tetFacei
                    << token::SPACE << p.tetPti
                    << token::SPACE << p.origProc
                    << token::SPACE << p.origId;
            }
            else
            {
                const std::size_t sizeofFields
                (
                    sizeof(particle::positionsCompat1706)
                  - offsetof(particle::positionsCompat1706, position)
                );

                os.write
                (
                    reinterpret_cast<const char*>(&p.position),
                    sizeofFields
                );
            }
            return os;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
