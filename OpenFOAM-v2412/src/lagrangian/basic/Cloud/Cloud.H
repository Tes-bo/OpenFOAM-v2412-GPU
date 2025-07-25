/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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
    Foam::Cloud

Description
    Base cloud calls templated on particle type

SourceFiles
    Cloud.C
    CloudIO.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_Cloud_H
#define Foam_Cloud_H

#include "cloud.H"
#include "IDLList.H"
#include "IOField.H"
#include "CompactIOField.H"
#include "polyMesh.H"
#include "bitSet.H"
#include "wordRes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class ParticleType> class Cloud;
template<class ParticleType> class IOPosition;


/*---------------------------------------------------------------------------*\
                           Class Cloud Declaration
\*---------------------------------------------------------------------------*/

template<class ParticleType>
class Cloud
:
    public cloud,
    public IDLList<ParticleType>
{
    // Private Data

        //- Reference to the mesh database
        const polyMesh& polyMesh_;

        //- Temporary storage for addressing. Used in findTris.
        mutable DynamicList<label> labels_;

        //- Does the cell have wall faces
        mutable autoPtr<bitSet> cellWallFacesPtr_;

        //- Temporary storage for the global particle positions
        mutable autoPtr<vectorField> globalPositionsPtr_;


    // Private Member Functions

        //- Check patches
        void checkPatches() const;

        //- Initialise cloud on IO constructor
        void initCloud(const bool checkClass);

        //- Find all cells which have wall faces
        void calcCellWallFaces() const;

        //- Read cloud properties dictionary
        void readCloudUniformProperties();

        //- Write cloud properties dictionary
        void writeCloudUniformProperties() const;


protected:

        //- Geometry type
        cloud::geometryType geometryType_;


public:

    friend class particle;
    template<class ParticleT>
    friend class IOPosition;

    typedef ParticleType particleType;

    //- Parcels are just particles
    typedef ParticleType parcelType;


    //- Runtime type information
    TypeName("Cloud");


    // Static Data

        //- Name of cloud properties dictionary
        static word cloudPropertiesName;


    // Constructors

        //- Construct without particles
        Cloud
        (
            const polyMesh& mesh,
            const Foam::zero,
            const word& cloudName
        );

        //- Construct from mesh and a list of particles
        Cloud
        (
            const polyMesh& mesh,
            const word& cloudName,
            const IDLList<ParticleType>& particles
        );

        //- Read construct from file(s) with given cloud instance.
        Cloud
        (
            const polyMesh& pMesh,
            const word& cloudName,
            //! Disable checking of class name (eg, for post-processing)
            const bool checkClass = true
        );


    // Member Functions

        // Access

            //- Return the polyMesh reference
            const polyMesh& pMesh() const noexcept
            {
                return polyMesh_;
            }

            //- Return the number of particles in the cloud
            using IDLList<ParticleType>::size;

            //- Return the number of particles in the cloud
            virtual label nParcels() const
            {
                return IDLList<ParticleType>::size();
            };

            //- Return temporary addressing
            DynamicList<label>& labels() const
            {
                return labels_;
            }


    // Iterators

        using typename IDLList<ParticleType>::iterator;
        using typename IDLList<ParticleType>::const_iterator;

        using IDLList<ParticleType>::begin;
        using IDLList<ParticleType>::cbegin;
        using IDLList<ParticleType>::end;
        using IDLList<ParticleType>::cend;


        // Edit

            //- Clear the particle list
            using IDLList<ParticleType>::clear;

            //- Transfer particle to cloud
            void addParticle(ParticleType* pPtr);

            //- Remove particle from cloud and delete
            void deleteParticle(ParticleType& p);

            //- Remove lost particles from cloud and delete
            void deleteLostParticles();

            //- Reset the particles
            void cloudReset(const Cloud<ParticleType>& c);

            //- Move the particles
            template<class TrackCloudType>
            void move
            (
                TrackCloudType& cloud,
                typename ParticleType::trackingData& td,
                const scalar trackTime
            );

            //- Remap the cells of particles corresponding to the
            //  mesh topology change
            void autoMap(const mapPolyMesh&);


        // Read

            //- Check lagrangian data field
            template<class DataType>
            void checkFieldIOobject
            (
                const Cloud<ParticleType>& c,
                const IOField<DataType>& data
            ) const;

            //- Check lagrangian data fieldfield
            template<class DataType>
            void checkFieldFieldIOobject
            (
                const Cloud<ParticleType>& c,
                const CompactIOField<Field<DataType>, DataType>& data
            ) const;

            //- Helper function to store a cloud field on its registry
            template<class Type>
            bool readStoreFile
            (
                const IOobject& io,
                const IOobject& ioNew
            ) const;

            //- Read from files into objectRegistry
            void readFromFiles
            (
                objectRegistry& obr,
                const wordRes& selectFields,
                const wordRes& excludeFields = wordRes::null()
            ) const;


        // Write

            //- Write the field data for the cloud of particles Dummy at
            //  this level.
            virtual void writeFields() const;

            //- Write using stream options.
            //  Only writes the cloud file if the Cloud isn't empty
            virtual bool writeObject
            (
                IOstreamOption streamOpt,
                const bool writeOnProc
            ) const;

            //- Write positions to \<cloudName\>_positions.obj file
            void writePositions() const;

            //- Call this before a topology change.
            //  Stores the particles global positions in the database
            //  for use during mapping.
            void storeGlobalPositions() const;


    // Housekeeping

        //- Helper to construct IOobject for field and current time.
        FOAM_DEPRECATED_FOR(2024-01, "newIOobject()")
        IOobject fieldIOobject
        (
            const word& fieldName,
            IOobjectOption::readOption rOpt = IOobjectOption::NO_READ
        ) const
        {
            return this->newIOobject(fieldName, rOpt);
        }
};


// Ostream Operator

template<class ParticleType>
Ostream& operator<<(Ostream& os, const Cloud<ParticleType>& c);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Cloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
