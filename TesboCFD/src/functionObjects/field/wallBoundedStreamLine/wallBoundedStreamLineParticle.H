/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::wallBoundedStreamLineParticle

Description
    Particle class that samples fields as it passes through. Used in streamline
    calculation.

SourceFiles
    wallBoundedStreamLineParticle.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_wallBoundedStreamLineParticle_H
#define Foam_wallBoundedStreamLineParticle_H

#include "wallBoundedParticle.H"
#include "interpolation.H"
#include "vectorList.H"
#include "InfoProxy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class wallBoundedStreamLineParticle;
class wallBoundedStreamLineParticleCloud;

Ostream& operator<<(Ostream&, const wallBoundedStreamLineParticle&);


/*---------------------------------------------------------------------------*\
                Class wallBoundedStreamLineParticle Declaration
\*---------------------------------------------------------------------------*/

class wallBoundedStreamLineParticle
:
    public wallBoundedParticle
{
public:

    //- Class used to pass tracking data to the trackToEdge function
    class trackingData
    :
        public wallBoundedParticle::trackingData
    {
    public:

        // Public Data

        const PtrList<interpolation<scalar>>& vsInterp_;
        const PtrList<interpolation<vector>>& vvInterp_;
        const interpolation<vector>& UInterp_;
        const scalar trackLength_;

        DynamicList<vectorList>& allPositions_;
        List<DynamicList<scalarList>>& allScalars_;
        List<DynamicList<vectorList>>& allVectors_;


        // Constructors

            template<class TrackCloudType>
            trackingData
            (
                TrackCloudType& cloud,
                const PtrList<interpolation<scalar>>& vsInterp,
                const PtrList<interpolation<vector>>& vvInterp,
                const interpolation<vector>& UInterp,
                const scalar trackLength,
                const bitSet& isWallPatch,

                DynamicList<List<point>>& allPositions,
                List<DynamicList<scalarList>>& allScalars,
                List<DynamicList<vectorList>>& allVectors
            )
            :
                wallBoundedParticle::trackingData(cloud, isWallPatch),
                vsInterp_(vsInterp),
                vvInterp_(vvInterp),
                UInterp_(UInterp),
                trackLength_(trackLength),

                allPositions_(allPositions),
                allScalars_(allScalars),
                allVectors_(allVectors)
            {}

            virtual ~trackingData() = default;
    };


protected:

    // Protected Data

        //- Track with +U or -U
        bool trackForward_;

        //- Lifetime of particle. Particle dies when reaches 0.
        label lifeTime_;

        //- Sampled positions
        DynamicList<point> sampledPositions_;

        //- Sampled scalars
        List<DynamicList<scalar>> sampledScalars_;

        //- Sampled vectors
        List<DynamicList<vector>> sampledVectors_;


    // Protected Member Functions

        vector interpolateFields
        (
            const trackingData& td,
            const point& position,
            const label celli,
            const label facei
        );

        vector sample(trackingData& td);


public:

    // Constructors

        //- Construct from components
        wallBoundedStreamLineParticle
        (
            const polyMesh& c,
            const point& position,
            const label celli,
            const label tetFacei,
            const label tetPti,
            const label meshEdgeStart,
            const label diagEdge,
            const bool trackForward,
            const label lifeTime
        );

        //- Construct from Istream
        wallBoundedStreamLineParticle
        (
            const polyMesh& c,
            Istream& is,
            bool readFields = true,
            bool newFormat = true
        );

        //- Construct copy
        wallBoundedStreamLineParticle(const wallBoundedStreamLineParticle& p);

        //- Return a clone
        autoPtr<particle> clone() const
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

            autoPtr<wallBoundedStreamLineParticle> operator()
            (
                Istream& is
            ) const
            {
                return autoPtr<wallBoundedStreamLineParticle>
                (
                    new wallBoundedStreamLineParticle(mesh_, is, true)
                );
            }
        };


    // Member Functions

        // Tracking

            //- Track all particles to their end point
            template<class TrackCloudType>
            bool move
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar trackTime
            );


        // I-O

            //- Read
            static void readFields(Cloud<wallBoundedStreamLineParticle>&);

            //- Write
            static void writeFields
            (
                const Cloud<wallBoundedStreamLineParticle>&
            );


    // Ostream Operator

        friend Ostream& operator<<
        (
            Ostream&,
            const wallBoundedStreamLineParticle&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "wallBoundedStreamLineParticleTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
