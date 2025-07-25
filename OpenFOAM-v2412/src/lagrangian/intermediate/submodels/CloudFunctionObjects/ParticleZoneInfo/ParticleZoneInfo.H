/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::ParticleZoneInfo

Group
    grpLagrangianIntermediateFunctionObjects

Description
    Reports cloud information for particles passing through a specified cell
    zone.

Usage
    Example usage:
    \verbatim
    cloudFunctions
    {
        particleZoneInfo1
        {
            // Mandatory entries
            type            particleZoneInfo;
            cellZone        leftFluid;

            // Optional entries
            writer          vtk;
        }
    }
    \endverbatim

    Results are written to file:
    - \<case\>/postProcessing/lagrangian/\<cloudName\>/\<functionName\>/\<time\>

    \verbatim
    # cellZone        : leftFluid
    # time            : 1.0000000000e+00
    #
    # origID    origProc    (x y z)    time0    age    d0    d    mass0    mass
    \endverbatim

    Where
    - origID : particle ID
    - origProc : processor ID
    - (x y z) : Cartesian co-ordinates
    - time0 : time particle enters the cellZone
    - age : time spent in the cellZone
    - d0 : diameter on entry to the cellZone
    - d : current diameter
    - mass0 : mass on entry to the cellZone
    - mass : current mass

    If the optional \c writer entry is supplied, cloud data is written in the
    specified format.

    During the run, output statistics are reported after the cloud solution,
    e.g.:

    \verbatim
    particleZoneInfo:
        Cell zone                       = leftFluid
        Contributions                   = 257
    \endverbatim

    Here, 'Contributions' refers to the number of incremental particle-move
    contributions recorded during this time step. At write times, the output
    is extended, e.g.:

    \verbatim
    particleZoneInfo:
        Cell zone                       = leftFluid
        Contributions                   = 822
        Number of particles             = 199
        Written data to "postProcessing/lagrangian/reactingCloud1/
    \endverbatim

SourceFiles
    ParticleZoneInfo.C

\*---------------------------------------------------------------------------*/

#ifndef ParticleZoneInfo_H
#define ParticleZoneInfo_H

#include "CloudFunctionObject.H"
#include "writeFile.H"
#include "coordSetWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        struct particleInfo Declaration
\*---------------------------------------------------------------------------*/

struct particleInfo
{
    label origID = -1;
    label origProc = -1;
    vector position = Zero;
    scalar time0 = 0;
    scalar age = 0;
    scalar d0 = 0;
    scalar d = 0;
    scalar mass0 = 0;
    scalar mass = 0;

    void operator+=(const particleInfo& p)
    {
        // Increment age
        age += p.age;

        // Set current values
        position = p.position;
        d = p.d;
        mass = p.mass;
    }

    scalar isOlderThan(const particleInfo& p) const
    {
        // Cannot just use time0 - particle may leave/re-enter and
        // so age is decoupled
        return (p.time0 + p.age) < (time0 + age);
    }

    friend bool operator==(const particleInfo& a, const particleInfo& b)
    {
        return
            a.origID == b.origID
         && a.origProc == b.origProc
         && a.position == b.position
         && a.time0 == b.time0
         && a.age == b.age
         && a.d0 == b.d0
         && a.d == b.d
         && a.mass0 == b.mass0
         && a.mass == b.mass;
    }

    friend bool operator!=(const particleInfo& a, const particleInfo& b)
    {
        return !(a == b);
    }


// IOstream Operators

    friend Istream& operator>>(Istream& is, particleInfo& pi)
    {
        is  >> pi.origID
            >> pi.origProc
            >> pi.position
            >> pi.time0
            >> pi.age
            >> pi.d0
            >> pi.d
            >> pi.mass0
            >> pi.mass;

        return is;
    }

    friend Ostream& operator<<(Ostream& os, const particleInfo& pi)
    {
        os  << pi.origID
            << " " << pi.origProc
            << " " << pi.position
            << " " << pi.time0
            << " " << pi.age
            << " " << pi.d0
            << " " << pi.d
            << " " << pi.mass0
            << " " << pi.mass;

        return os;
    }
};

/*---------------------------------------------------------------------------*\
                       Class ParticleZoneInfo Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ParticleZoneInfo
:
    public CloudFunctionObject<CloudType>,
    public functionObjects::writeFile

{
    // Private Data

        // Typedefs

            //- Convenience typedef for parcel type
            typedef typename CloudType::parcelType parcelType;


        //- Cell zone name
        word cellZoneName_;

        //- Cell zone index
        label cellZoneId_;

        //- Stored data
        DynamicList<particleInfo> data_;

        //- Work storage
        DynamicList<particleInfo> movedParticles_;

        //- Maximum particle ID per processor
        labelList maxIDs_;

        //- Set writer
        autoPtr<coordSetWriter> writerPtr_;


    // Private Member Functions

        //- Write output file header
        void writeFileHeader(Ostream& os) const;

        //- Return true if celli is in the cellZone
        bool inZone(const label celli) const;

        //- Return the index of the particle in the storage (data_)
        //- Returns -1 if not found
        label getParticleID(const particleInfo& p) const;

        //- Write fields using writerPtr_
        void writeWriter(const UList<particleInfo>& data);


public:

    //- Runtime type information
    TypeName("particleZoneInfo");


    // Constructors

        //- Construct from dictionary
        ParticleZoneInfo
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        ParticleZoneInfo(const ParticleZoneInfo<CloudType>& pe);

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObject<CloudType>> clone() const
        {
            return autoPtr<CloudFunctionObject<CloudType>>
            (
                new ParticleZoneInfo<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ParticleZoneInfo() = default;


    // Member Functions

        //- Pre-evolve hook
        virtual void preEvolve
        (
            const typename parcelType::trackingData& td
        );

        //- Post-evolve hook
        virtual void postEvolve
        (
            const typename parcelType::trackingData& td
        );

        //- Post-move hook
        virtual bool postMove
        (
            parcelType& p,
            const scalar dt,
            const point& position0,
            const typename parcelType::trackingData& td
        );

        //- Write
        virtual void write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ParticleZoneInfo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
