/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2022-2024 OpenCFD Ltd.
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
    Foam::IOmapDistributePolyMesh

Description
    IOmapDistributePolyMesh is derived from mapDistributePolyMesh and
    IOobject to give the mapDistributePolyMesh
    automatic IO functionality via the objectRegistry.

SourceFiles
    IOmapDistributePolyMesh.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_IOmapDistributePolyMesh_H
#define Foam_IOmapDistributePolyMesh_H

#include "mapDistributePolyMesh.H"
#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class IOmapDistributePolyMesh Declaration
\*---------------------------------------------------------------------------*/

class IOmapDistributePolyMesh
:
    public regIOobject,
    public mapDistributePolyMesh
{
    // Private Member Functions

        //- Read if IOobject flags set. Return true if read.
        bool readIOcontents();

public:

    //- Runtime type information
    TypeName("mapDistributePolyMesh");


    // Constructors

        //- Construct from IOobject
        explicit IOmapDistributePolyMesh(const IOobject& io);

        //- Construct from IOobject, copying mapDistributePolyMesh contents
        IOmapDistributePolyMesh
        (
            const IOobject& io,
            const mapDistributePolyMesh& map
        );

        //- Construct from IOobject, moving mapDistributePolyMesh contents
        IOmapDistributePolyMesh
        (
            const IOobject& io,
            mapDistributePolyMesh&& map
        );


    //- Destructor
    virtual ~IOmapDistributePolyMesh() = default;


    // Member Functions

        //- The readData method for regIOobject read operation
        virtual bool readData(Istream& is);

        //- The writeData method for regIOobject write operation
        virtual bool writeData(Ostream& os) const;
};


/*---------------------------------------------------------------------------*\
                 Class IOmapDistributePolyMeshRef Declaration
\*---------------------------------------------------------------------------*/

//- A IOmapDistributePolyMesh wrapper for using referenced external data.
class IOmapDistributePolyMeshRef
:
    public regIOobject
{
    // Private Data

        //- Reference to the external content
        refPtr<mapDistributePolyMesh> contentRef_;

public:

    //- Type is identical to IOmapDistributePolyMesh
    virtual const word& type() const
    {
        return IOmapDistributePolyMesh::typeName;
    }


    // Generated Methods

        //- No default construct
        IOmapDistributePolyMeshRef() = delete;

        //- No copy construct
        IOmapDistributePolyMeshRef(const IOmapDistributePolyMeshRef&) = delete;

        //- No copy assignment
        void operator=(const IOmapDistributePolyMeshRef&) = delete;


    // Constructors

        //- Construct from IOobject and const data reference
        IOmapDistributePolyMeshRef
        (
            const IOobject& io,
            const mapDistributePolyMesh& map
        );


    //- Destructor
    virtual ~IOmapDistributePolyMeshRef() = default;


    // Member Functions

        //- Allow cast to const content
        //  Fatal if content is not set
        operator const mapDistributePolyMesh&() const
        {
            return contentRef_.cref();
        }

        //- The readData method for regIOobject read operation
        //  Fatal if content is constant (or not set)
        virtual bool readData(Istream& is);

        //- The writeData method for regIOobject write operation
        //  Fatal if content is not set
        virtual bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
