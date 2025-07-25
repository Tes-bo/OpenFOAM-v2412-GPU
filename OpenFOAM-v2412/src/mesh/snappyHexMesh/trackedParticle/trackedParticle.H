/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::trackedParticle

Description
    Particle class that marks cells it passes through. Used to mark cells
    visited by feature edges.

SourceFiles
    trackedParticle.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_trackedParticle_H
#define Foam_trackedParticle_H

#include "bitSet.H"
#include "particle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations

class trackedParticleCloud;
class trackedParticle;

Ostream& operator<<(Ostream&, const trackedParticle&);


/*---------------------------------------------------------------------------*\
                     Class trackedParticle Declaration
\*---------------------------------------------------------------------------*/

class trackedParticle
:
    public particle
{
    // Private Data

        //- Start point to track from
        point start_;

        //- End point to track to
        point end_;

        //- Level of this particle
        label level_;

        //- Passive label (used to store feature edge mesh)
        label i_;

        //- Passive label (used to store feature edge point)
        label j_;

        //- Passive label (used to store feature edge label)
        label k_;


public:

    friend class Cloud<trackedParticle>;

    //- Class used to pass tracking data to the trackToFace function
    class trackingData
    :
        public particle::trackingData
    {
    public:

        labelList& maxLevel_;

        List<bitSet>& featureEdgeVisited_;


        // Constructors

            trackingData
            (
                Cloud<trackedParticle>& cloud,
                labelList& maxLevel,
                List<bitSet>& featureEdgeVisited
            )
            :
                particle::trackingData(cloud),
                maxLevel_(maxLevel),
                featureEdgeVisited_(featureEdgeVisited)
            {}

    };


    // Static Data Members

        //- Size in bytes of the fields
        static const std::size_t sizeofFields_;


    // Constructors

        //- Construct from components
        trackedParticle
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPtI,
            const point& end,
            const label level,
            const label i,
            const label j,
            const label k
        );

        //- Construct from a position and a cell,
        //- searching for the rest of the required topology
        trackedParticle
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli,
            const point& end,
            const label level,
            const label i,
            const label j,
            const label k
        );

        //- Construct from Istream
        trackedParticle
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true,
            bool newFormat = true
        );

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

            autoPtr<trackedParticle> operator()(Istream& is) const
            {
                return autoPtr<trackedParticle>::New(mesh_, is, true);
            }
        };


    // Member Functions

        //- Point to track from
        point& start() noexcept { return start_; }

        //- Point to track to
        point& end() noexcept { return end_; }

        //- Transported label
        label i() const noexcept { return i_; }

        //- Transported label
        label& i() noexcept { return i_; }

        //- Transported label
        label j() const noexcept { return j_; }

        //- Transported label
        label& j() noexcept { return j_; }

        //- Transported label
        label k() const noexcept { return k_; }

        //- Transported label
        label& k() noexcept { return k_; }



        // Tracking

            //- Track all particles to their end point
            bool move(Cloud<trackedParticle>&, trackingData&, const scalar);

            //- Overridable function to handle the particle hitting a patch
            //  Executed before other patch-hitting functions
            bool hitPatch(Cloud<trackedParticle>&, trackingData&);

            //- Overridable function to handle the particle hitting a wedge
            void hitWedgePatch(Cloud<trackedParticle>&, trackingData&);

            //- Overridable function to handle the particle hitting a
            //  symmetry plane
            void hitSymmetryPlanePatch(Cloud<trackedParticle>&, trackingData&);

            //- Overridable function to handle the particle hitting a
            //  symmetry patch
            void hitSymmetryPatch(Cloud<trackedParticle>&, trackingData&);

            //- Overridable function to handle the particle hitting a cyclic
            void hitCyclicPatch(Cloud<trackedParticle>&, trackingData&);

            //- Overridable function to handle the particle hitting a cyclicAMI
            void hitCyclicAMIPatch
            (
                Cloud<trackedParticle>&,
                trackingData&,
                const vector&
            );

            //- Overridable function to handle the particle hitting a cyclicACMI
            void hitCyclicACMIPatch
            (
                Cloud<trackedParticle>&,
                trackingData&,
                const vector&
            );

            //- Overridable function to handle the particle hitting a
            //- processorPatch
            void hitProcessorPatch(Cloud<trackedParticle>&, trackingData&);

            //- Overridable function to handle the particle hitting a wallPatch
            void hitWallPatch(Cloud<trackedParticle>&, trackingData&);

            //- Convert processor patch addressing to the global equivalents
            //  and set the celli to the face-neighbour
            void correctAfterParallelTransfer(const label, trackingData&);


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const trackedParticle&);
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Contiguous data for trackedParticle
template<> struct is_contiguous<trackedParticle> : std::true_type {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
