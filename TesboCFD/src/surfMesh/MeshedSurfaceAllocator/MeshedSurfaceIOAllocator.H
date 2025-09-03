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
    Foam::Detail::MeshedSurfaceIOAllocator

Description
    A helper class for storing points and faces with IO capabilities.

SourceFiles
    MeshedSurfaceIOAllocator.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_MeshedSurfaceIOAllocator_H
#define Foam_MeshedSurfaceIOAllocator_H

#include "pointIOField.H"
#include "faceIOList.H"  // faceIOList and faceCompactIOList

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Detail
{

/*---------------------------------------------------------------------------*\
              Class Detail::MeshedSurfaceIOAllocator Declaration
\*---------------------------------------------------------------------------*/

class MeshedSurfaceIOAllocator
{
    // Private Data

        //- Points
        pointIOField points_;

        //- Faces
        faceCompactIOList faces_;

public:

    // Generated Methods

        //- No copy construct
        MeshedSurfaceIOAllocator(const MeshedSurfaceIOAllocator&) = delete;

        //- No copy assignment
        void operator=(const MeshedSurfaceIOAllocator&) = delete;


    // Constructors

        //- Read construct from IOobjects
        MeshedSurfaceIOAllocator
        (
            const IOobject& ioPoints,
            const IOobject& ioFaces
        );

        //- Construct from IOobjects, copying components
        MeshedSurfaceIOAllocator
        (
            const IOobject& ioPoints, const pointField& points,
            const IOobject& ioFaces,  const faceList& faces
        );

        //- Construct from IOobjects, moving components
        MeshedSurfaceIOAllocator
        (
            const IOobject& ioPoints, pointField&& points,
            const IOobject& ioFaces,  faceList&& faces
        );


    //- Destructor
    virtual ~MeshedSurfaceIOAllocator() = default;


    // Member Functions

    // Database

        //- Set the instance for mesh files
        void setInstance(const fileName& inst);

        //- Adjust the write option for all components
        void setWriteOption(IOobjectOption::writeOption wOpt);


    // Access

        //- Non-const access to the points
        pointIOField& storedIOPoints() noexcept
        {
            return points_;
        }

        //- Non-const access to the faces
        faceCompactIOList& storedIOFaces() noexcept
        {
            return faces_;
        }

        //- Const access to the points
        const pointIOField& storedIOPoints() const noexcept
        {
            return points_;
        }

        //- Const access to the faces
        const faceCompactIOList& storedIOFaces() const noexcept
        {
            return faces_;
        }


    // Storage Management

        //- Clear primitive data (points, faces)
        void clear();


    // Writing

        //- Write using stream options
        virtual bool writeObject
        (
            IOstreamOption streamOpt,
            const bool writeOnProc
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Detail
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
