/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
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
    Foam::findCellParticle

Description
    Particle class that finds cells by tracking

SourceFiles
    findCellParticle.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_findCellParticle_H
#define Foam_findCellParticle_H

#include "particle.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class findCellParticleCloud;
class findCellParticle;

Ostream& operator<<(Ostream&, const findCellParticle&);


/*---------------------------------------------------------------------------*\
                     Class findCellParticle Declaration
\*---------------------------------------------------------------------------*/

class findCellParticle
:
    public particle
{
    // Private data

        //- Start point to track from
        point start_;

        //- End point to track to
        point end_;

        //- Passive data
        label data_;


public:

    friend class Cloud<findCellParticle>;

    //- Class used to pass tracking data to the trackToFace function
    class trackingData
    :
        public particle::trackingData
    {
        labelListList& cellToData_;
        List<List<point>>& cellToEnd_;

    public:

        // Constructors

            trackingData
            (
                Cloud<findCellParticle>& cloud,
                labelListList& cellToData,
                List<List<point>>& cellToEnd
            )
            :
                particle::trackingData(cloud),
                cellToData_(cellToData),
                cellToEnd_(cellToEnd)
            {}


        // Member Functions

            labelListList& cellToData() noexcept { return cellToData_; }

            List<List<point>>& cellToEnd() noexcept { return cellToEnd_; }
    };


    // Constructors

        //- Construct from components
        findCellParticle
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti,
            const point& end,
            const label data
        );

        //- Construct from a position and a cell,
        //- searching for the rest of the required topology
        findCellParticle
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli,
            const point& end,
            const label data
        );

        //- Construct from Istream
        findCellParticle
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

            autoPtr<findCellParticle> operator()(Istream& is) const
            {
                return autoPtr<findCellParticle>
                (
                    new findCellParticle(mesh_, is, true)
                );
            }
        };


    // Member Functions

        //- Point to track from
        const point& start() const noexcept { return start_; }

        //- Point to track from
        point& start() noexcept { return start_; }

        //- Point to track to
        const point& end() const noexcept { return end_; }

        //- Point to track to
        point& end() noexcept { return end_; }

        //- Transported label
        label data() const noexcept { return data_; }

        //- Transported label
        label& data() noexcept { return data_; }


        // Tracking

            //- Track all particles to their end point
            bool move(Cloud<findCellParticle>&, trackingData&, const scalar);

            //- Overridable function to handle the particle hitting a patch
            //  Executed before other patch-hitting functions
            bool hitPatch(Cloud<findCellParticle>&, trackingData&);

            //- Overridable function to handle the particle hitting a wedge
            void hitWedgePatch(Cloud<findCellParticle>&, trackingData&);

            //- Overridable function to handle the particle hitting a
            //  symmetry plane
            void hitSymmetryPlanePatch(Cloud<findCellParticle>&, trackingData&);

            //- Overridable function to handle the particle hitting a
            //  symmetry patch
            void hitSymmetryPatch(Cloud<findCellParticle>&, trackingData&);

            //- Overridable function to handle the particle hitting a cyclic
            void hitCyclicPatch(Cloud<findCellParticle>&, trackingData&);

            //- Overridable function to handle the particle hitting a cyclicAMI
            void hitCyclicAMIPatch
            (
                Cloud<findCellParticle>&,
                trackingData&,
                const vector&
            );

            //- Overridable function to handle the particle hitting a cyclicACMI
            void hitCyclicACMIPatch
            (
                Cloud<findCellParticle>&,
                trackingData&,
                const vector&
            );

            //- Overridable function to handle the particle hitting a
            //- processorPatch
            void hitProcessorPatch(Cloud<findCellParticle>&, trackingData&);

            //- Overridable function to handle the particle hitting a wallPatch
            void hitWallPatch(Cloud<findCellParticle>&, trackingData&);


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const findCellParticle&);
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Contiguous data for findCellParticle
template<> struct is_contiguous<findCellParticle> : std::true_type {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
