/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2024 OpenCFD Ltd.
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
    Foam::solidParticle

Description
    Simple solid spherical particle class with one-way coupling with the
    continuous phase.

SourceFiles
    solidParticleI.H
    solidParticle.C
    solidParticleIO.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_solidParticle_H
#define Foam_solidParticle_H

#include "particle.H"
#include "IOstream.H"
#include "interpolationCellPoint.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class solidParticle;
class solidParticleCloud;

Ostream& operator<<(Ostream&, const solidParticle&);

/*---------------------------------------------------------------------------*\
                           Class solidParticle Declaration
\*---------------------------------------------------------------------------*/

class solidParticle
:
    public particle
{
    // Private Data

        //- Diameter
        scalar d_;

        //- Velocity of parcel
        vector U_;


public:

    friend class Cloud<solidParticle>;

    //- Class used to pass tracking data to the trackToFace function
    class trackingData
    :
        public particle::trackingData
    {
        // Interpolators for continuous phase fields

            const interpolationCellPoint<scalar>& rhoInterp_;
            const interpolationCellPoint<vector>& UInterp_;
            const interpolationCellPoint<scalar>& nuInterp_;

        //- Local gravitational or other body-force acceleration
        const vector& g_;


    public:

        // Constructors

            trackingData
            (
                const solidParticleCloud& cld,
                const interpolationCellPoint<scalar>& rhoInterp,
                const interpolationCellPoint<vector>& UInterp,
                const interpolationCellPoint<scalar>& nuInterp,
                const vector& g
            )
            :
                particle::trackingData(cld),
                rhoInterp_(rhoInterp),
                UInterp_(UInterp),
                nuInterp_(nuInterp),
                g_(g)
            {}


        // Member Functions

            const interpolationCellPoint<scalar>& nuInterp() const noexcept
            {
                return nuInterp_;
            }

            const interpolationCellPoint<vector>& UInterp() const noexcept
            {
                return UInterp_;
            }

            const interpolationCellPoint<scalar>& rhoInterp() const noexcept
            {
                return rhoInterp_;
            }

            const vector& g() const noexcept { return g_; }
    };


    // Static data members

        //- Size in bytes of the fields
        static const std::size_t sizeofFields;


    // Constructors

        //- Construct from a position and a cell
        //  Searches for the rest of the required topology.
        //  Other properties are zero initialised.
        solidParticle
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli = -1
        )
        :
            particle(mesh, position, celli),
            d_(0),
            U_(Zero)
        {}

        //- Construct from components
        solidParticle
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti,
            const scalar d,
            const vector& U
        )
        :
            particle(mesh, coordinates, celli, tetFacei, tetPti),
            d_(d),
            U_(U)
        {}

        //- Construct from Istream
        solidParticle
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

            autoPtr<solidParticle> operator()(Istream& is) const
            {
                return autoPtr<solidParticle>::New(mesh_, is, true);
            }
        };


    // Member Functions

    // Access

        //- Return diameter
        scalar d() const noexcept { return d_; }

        //- Return velocity
        const vector& U() const noexcept { return U_; }


    // Tracking

        //- Move
        bool move(solidParticleCloud&, trackingData&, const scalar);


        // Patch interactions

            //- Overridable function to handle the particle hitting a patch
            //  Executed before other patch-hitting functions
            bool hitPatch(solidParticleCloud& cloud, trackingData& td);

            //- Overridable function to handle the particle hitting a
            //  processorPatch
            void hitProcessorPatch(solidParticleCloud& cloud, trackingData& td);

            //- Overridable function to handle the particle hitting a wallPatch
            void hitWallPatch(solidParticleCloud& cloud, trackingData& td);

            //- Transform the physical properties of the particle
            //  according to the given transformation tensor
            virtual void transformProperties(const tensor& T);

            //- Transform the physical properties of the particle
            //  according to the given separation vector
            virtual void transformProperties(const vector& separation);


    // I-O

        static void readFields(Cloud<solidParticle>& c);

        static void writeFields(const Cloud<solidParticle>& c);


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const solidParticle&);
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Contiguous data for solidParticle
template<> struct is_contiguous<solidParticle> : std::true_type {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
