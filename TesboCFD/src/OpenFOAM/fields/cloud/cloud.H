/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::cloud

Description
    A cloud is a registry collection of lagrangian particles.

SourceFiles
    cloud.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_cloud_H
#define Foam_cloud_H

#include "objectRegistry.H"
#include "Enum.H"
#include "point.H"
#include "IOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                            Class cloud Declaration
\*---------------------------------------------------------------------------*/

class cloud
:
    public objectRegistry
{
public:

    //- Cloud geometry type (internal or IO representations)
    enum class geometryType
    {
        COORDINATES,    //!< barycentric coordinates
        POSITIONS       //!< positions
    };

    //- Named enumerations "coordinates", "positions"
    static const Enum<geometryType> geometryTypeNames;


    // Static Data Members

        //- The prefix to local: %lagrangian
        static const word prefix;

        //- The default cloud name: %defaultCloud
        static word defaultName;


    //- Runtime type information
    TypeName("cloud");


    // Generated Methods

        //- No copy construct
        cloud(const cloud&) = delete;

        //- No copy assignment
        void operator=(const cloud&) = delete;


    // Constructors

        //- Construct for given objectRegistry. Uses default cloud name
        explicit cloud(const objectRegistry& obr);

        //- Construct for given objectRegistry and named cloud
        cloud(const objectRegistry& obr, const word& cloudName);


    //- Destructor
    virtual ~cloud() = default;


    // Member Functions

        // Sizes

            //- Number of parcels for the hosting cloud
            virtual label nParcels() const;


        // Edit

            //- Remap the cells of particles corresponding to the
            //- mesh topology change
            virtual void autoMap(const mapPolyMesh&);


        // I-O

            //- Read particle fields from objects in the obr registry
            virtual void readObjects(const objectRegistry& obr);

            //- Write particle fields as objects into the obr registry
            virtual void writeObjects(objectRegistry& obr) const;

            //- Helper to construct IOField on a supplied object registry
            template<class Type>
            static IOField<Type>& createIOField
            (
                const word& fieldName,
                const label nParticle,
                objectRegistry& obr
            );

            //- Locate an IOField within object registry
            //  \return nullptr if not found or wrong type
            template<class Type>
            static const IOField<Type>* findIOField
            (
                const word& fieldName,
                const objectRegistry& obr
            )
            {
                return obr.cfindObject<IOField<Type>>(fieldName);
            }

            //- Locate the "position" IOField within object registry
            //  \return nullptr if not found or wrong type
            static const IOField<point>* findIOPosition
            (
                const objectRegistry& obr
            )
            {
                return obr.cfindObject<IOField<point>>("position");
            }

            //- Lookup an IOField within object registry
            //  Fatal if not found or wrong type
            template<class Type>
            static const IOField<Type>& lookupIOField
            (
                const word& fieldName,
                const objectRegistry& obr
            )
            {
                return obr.lookupObject<IOField<Type>>(fieldName);
            }

            //- Lookup an IOField within object registry
            //  Fatal if not found or wrong type
            //
            //  Note: const char signature to avoid spurious
            //  -Wdangling-reference with gcc-13
            template<class Type>
            static const IOField<Type>& lookupIOField
            (
                const char* fieldName,
                const objectRegistry& obr
            )
            {
                return obr.lookupObject<IOField<Type>>(word(fieldName));
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "cloudTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
