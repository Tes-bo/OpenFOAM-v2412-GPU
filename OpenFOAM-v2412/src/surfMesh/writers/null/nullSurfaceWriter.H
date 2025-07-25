/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    Foam::surfaceWriters::nullWriter

Description
    Suppresses output of geometry and fields

SourceFiles
    nullSurfaceWriter.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_surfaceWriters_nullWriter_H
#define Foam_surfaceWriters_nullWriter_H

#include "surfaceWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceWriters
{

/*---------------------------------------------------------------------------*\
                         Class nullWriter Declaration
\*---------------------------------------------------------------------------*/

class nullWriter
:
    public surfaceWriter
{
public:

    //- Declare type-name, virtual type (without debug switch)
    TypeNameNoDebug("none");


    // Constructors

        //- Default construct
        nullWriter();

        //- Construct with some output options
        explicit nullWriter(const dictionary& options);

        //- Construct from components
        nullWriter
        (
            const meshedSurf& surf,
            const fileName& outputPath,
            bool parallel = UPstream::parRun(),
            const dictionary& options = dictionary()
        );

        //- Construct from components
        nullWriter
        (
            const pointField& points,
            const faceList& faces,
            const fileName& outputPath,
            bool parallel = UPstream::parRun(),
            const dictionary& options = dictionary()
        );


    //- Destructor
    virtual ~nullWriter() = default;


    // Member Functions

    // Capability

        //- False: never needs an update.
        virtual bool needsUpdate() const;

        //- True: like a /dev/null device.
        virtual bool wroteData() const;

        //- False: The null writer is never enabled, which lets the caller
        //- skip various (possibly expensive) preparatory operations.
        virtual bool enabled() const;


    // Surface Association

        //- Change association with a surface (no-op).
        virtual void setSurface
        (
            const meshedSurf& s,
            bool parallel
        ); // override

        //- Change association with a surface (no-op).
        virtual void setSurface
        (
            const pointField& points,
            const faceList& faces,
            bool parallel
        ); // override


    // Output

        //- Open for output on specified path, using existing surface (no-op)
        virtual void open(const fileName& outputPath); // override


    // Write

        virtual fileName write(); // override;

        declareSurfaceWriterWriteMethod(label);
        declareSurfaceWriterWriteMethod(scalar);
        declareSurfaceWriterWriteMethod(vector);
        declareSurfaceWriterWriteMethod(sphericalTensor);
        declareSurfaceWriterWriteMethod(symmTensor);
        declareSurfaceWriterWriteMethod(tensor);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceWriters
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
