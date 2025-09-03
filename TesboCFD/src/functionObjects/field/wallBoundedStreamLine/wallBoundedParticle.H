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
    Foam::wallBoundedParticle

Description
    Particle class that tracks on triangles of boundary faces. Use
    trackToEdge similar to trackToFace on particle.

SourceFiles
    wallBoundedParticle.C
    wallBoundedParticleTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_wallBoundedParticle_H
#define Foam_wallBoundedParticle_H

#include "particle.H"
#include "InfoProxy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class wallBoundedParticle;

Ostream& operator<<(Ostream&, const wallBoundedParticle&);
Ostream& operator<<(Ostream&, const InfoProxy<wallBoundedParticle>&);


/*---------------------------------------------------------------------------*\
                     Class wallBoundedParticle Declaration
\*---------------------------------------------------------------------------*/

class wallBoundedParticle
:
    public particle
{
public:

    //- Class used to pass tracking data to the trackToFace function
    class trackingData
    :
        public particle::trackingData
    {
    public:

        const bitSet& isWallPatch_;

        // Constructors

            template<class TrackCloudType>
            trackingData
            (
                const TrackCloudType& cloud,
                const bitSet& isWallPatch
            )
            :
                particle::trackingData(cloud),
                isWallPatch_(isWallPatch)
            {}
    };


protected:

    // Protected Data

        //- Particle position is updated locally as opposed to via track
        //  functions of the base Foam::particle class
        point localPosition_;

        //- Particle is on mesh edge:
        //      const face& f = mesh.faces()[tetFace()]
        //      const edge e(f[meshEdgeStart_], f.nextLabel(meshEdgeStart_));
        //  Note that this real edge
        //  is also one of the edges of the face-triangle (from
        //  tetFace()+tetPt()).
        label meshEdgeStart_;

        //- Particle is on diagonal edge:
        //      const face& f = mesh.faces()[tetFace()]
        //      label faceBasePti = mesh.tetBasePtIs()[facei];
        //      label diagPti = (faceBasePti+diagEdge_)%f.size();
        //      const edge e(f[faceBasePti], f[diagPti]);
        label diagEdge_;


    // Protected Member Functions

        //- Construct current edge
        edge currentEdge() const;

        //- Replacement for particle::crossEdgeConnectedFace that avoids bombing
        //  out on invalid tet decomposition (tetBasePtIs = -1)
        void crossEdgeConnectedFace
        (
            const label& celli,
            label& tetFacei,
            label& tetPti,
            const edge& e
        );

        //- Cross mesh edge into different face on same cell
        void crossEdgeConnectedFace(const edge& meshEdge);

        //- Cross diagonal edge into different triangle on same face,cell
        void crossDiagonalEdge();

        //- Track through single triangle
        scalar trackFaceTri(const vector& n, const vector& endPosition, label&);

        //- Is current triangle in the track direction
        bool isTriAlongTrack(const vector& n, const point& endPosition) const;


public:

    // Static Data Members

        //- Size in bytes of the fields
        static const std::size_t sizeofFields_;


    // Constructors

        //- Construct from components
        wallBoundedParticle
        (
            const polyMesh& c,
            const point& position,
            const label celli,
            const label tetFacei,
            const label tetPti,
            const label meshEdgeStart,
            const label diagEdge
        );

        //- Construct from Istream
        wallBoundedParticle
        (
            const polyMesh& c,
            Istream& is,
            bool readFields = true,
            bool newFormat = true
        );

        //- Construct copy
        wallBoundedParticle(const wallBoundedParticle& p);

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

            autoPtr<wallBoundedParticle> operator()
            (
                Istream& is
            ) const
            {
                return autoPtr<wallBoundedParticle>
                (
                    new wallBoundedParticle(mesh_, is, true)
                );
            }
        };


    // Member Functions

    // Access

        //- The mesh edge label or -1
        label meshEdgeStart() const noexcept { return meshEdgeStart_; }

        //- The diagonal edge label or -1
        label diagEdge() const noexcept { return diagEdge_; }


        // Track

            //- Equivalent of trackToFace
            template<class TrackCloudType>
            scalar trackToEdge
            (
                TrackCloudType& cloud,
                trackingData& td,
                const vector& endPosition
            );


            // Patch interactions

                //- Do all patch interaction
                template<class TrackCloudType>
                void patchInteraction
                (
                    TrackCloudType& cloud,
                    trackingData& td,
                    const scalar trackFraction
                );

                //- Overridable function to handle the particle hitting a
                //- processorPatch
                template<class TrackCloudType>
                void hitProcessorPatch
                (
                    TrackCloudType& cloud,
                    trackingData& td
                );

                //- Overridable function to handle the particle hitting a wallPatch
                template<class TrackCloudType>
                void hitWallPatch
                (
                    TrackCloudType& cloud,
                    trackingData& td
                );


    // Info

        //- Return info proxy,
        //- used to print particle information to a stream
        InfoProxy<wallBoundedParticle> info() const noexcept
        {
            return *this;
        }


    // I-O

        //- Read
        template<class CloudType>
        static void readFields(CloudType&);

        //- Write
        template<class CloudType>
        static void writeFields(const CloudType&);


    // Ostream Operator

        friend Ostream& operator<<
        (
            Ostream&,
            const wallBoundedParticle&
        );

        friend Ostream& operator<<
        (
            Ostream&,
            const InfoProxy<wallBoundedParticle>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "wallBoundedParticleTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
