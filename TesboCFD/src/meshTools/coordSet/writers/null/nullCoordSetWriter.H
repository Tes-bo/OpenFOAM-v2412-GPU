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
    Foam::coordSetWriters::nullWriter

Description
    A coordSet(s) writer with suppressed output.

SourceFiles
    nullCoordSetWriter.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_coordSetWriters_nullWriter_H
#define Foam_coordSetWriters_nullWriter_H

#include "coordSetWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace coordSetWriters
{

/*---------------------------------------------------------------------------*\
                         Class nullWriter Declaration
\*---------------------------------------------------------------------------*/

class nullWriter
:
    public coordSetWriter
{
    // Private Member Functions

        //- Templated write operation
        template<class Type>
        fileName writeTemplate
        (
            const word& fieldName,      //!< Name of field
            const Field<Type>& values   //!< Local field values to write
        )
        {
            return fileName();
        }

        //- Templated write operation
        template<class Type>
        fileName writeTemplate
        (
            const word& fieldName,
            const List<Field<Type>>& fieldValues
        )
        {
            return fileName();
        }


public:

    //- Runtime type information (no debug)
    TypeNameNoDebug("none");


    // Constructors

        //- Default construct
        nullWriter();

        //- Construct with some output options
        explicit nullWriter(const dictionary& options);

        //- Construct from components
        nullWriter
        (
            const coordSet& coords,
            const fileName& outputPath,
            const dictionary& options = dictionary()
        );

        //- Construct from components
        nullWriter
        (
            const UPtrList<coordSet>& tracks,
            const fileName& outputPath,
            const dictionary& options = dictionary()
        );


    //- Destructor
    virtual ~nullWriter();


    // Member Functions

    // Capability

        //- False: The null writer is never enabled, which lets the caller
        //- skip various (possibly expensive) preparatory operations.
        virtual bool enabled() const;

        //- False: no internal buffering possible
        virtual bool buffering() const;  // override

        //- False: never needs an update.
        virtual bool needsUpdate() const;

        //- True: like a /dev/null device.
        virtual bool wroteData() const;


    // Content Association

        //- Set coordinates (no-op).
        virtual void setCoordinates(const coordSet*); // override

        //- Set coordinates (no-op).
        virtual void setCoordinates(const coordSet&); // override

        //- Set tracks (no-op).
        virtual void setTracks(const UPtrList<coordSet>&); // override


    // Output

        //- Characteristic output file name - information only.
        //- Always an empty fileName
        virtual fileName path() const;  // override

        //- Open for output on specified path, using existing content (no-op)
        virtual void open(const fileName& outputPath); // override


    // Member Functions

        declareCoordSetWriterWriteMethod(label);
        declareCoordSetWriterWriteMethod(scalar);
        declareCoordSetWriterWriteMethod(vector);
        declareCoordSetWriterWriteMethod(sphericalTensor);
        declareCoordSetWriterWriteMethod(symmTensor);
        declareCoordSetWriterWriteMethod(tensor);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace coordSetWriters
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
