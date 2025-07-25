/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2023 OpenCFD Ltd.
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
    Foam::surfaceWriters::mydebugWriter

Description
    A surfaceWriter for special purpose debugging.
    Its definition and behaviour are subject to change at any time.

    \verbatim
    formatOptions
    {
        debug
        {
            merge   true;
            write   true;
        }
    }
    \endverbatim

    Format options:
    \table
        Property    | Description                           | Required | Default
        commsType   | scheduled/nonBlocking/buffered        | no  | scheduled
        merge       | Enable geometry/field merging         | no  | true
        write       | Write file(s)                         | no  | false
        narrow      | Communicate with narrowed values      | no  | false
    \endtable

Note
    Disabling geometry/field merging (in parallel) implicitly deactivates
    writing as well.

    Output files (if any) are written as boundaryData (binary + header).

SourceFiles
    mydebugSurfaceWriter.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_mydebugSurfaceWriter_H
#define Foam_mydebugSurfaceWriter_H

#include "surfaceWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class regIOobject;

namespace surfaceWriters
{

/*---------------------------------------------------------------------------*\
                        Class debugWriter Declaration
\*---------------------------------------------------------------------------*/

class mydebugWriter
:
    public surfaceWriter
{
    // Private Data

        //- Enable/disable all merging in parallel
        bool enableMerge_;

        //- Output writing?
        bool enableWrite_;

        //- Output files with FoamFile header
        bool header_;

        //- Communicate with narrowed values
        bool narrowTransfer_;

        //- Output stream option
        IOstreamOption streamOpt_;


    // Private Member Functions

        //- Write serial surface geometry to "points" file.
        void serialWriteGeometry(const regIOobject&, const meshedSurf& surf);

        //- Gather (merge) fields with renumbering and shrinking for point data
        template<class Type>
        tmp<Field<Type>> mergeField(const Field<Type>& fld) const;

        //- Templated write operation
        template<class Type>
        fileName writeTemplate
        (
            const word& fieldName,          //!< Name of field
            const Field<Type>& localValues  //!< Local field values to write
        );


public:

    //- Declare type-name, virtual type (without debug switch)
    TypeNameNoDebug("mydebug");


    // Constructors

        //- Default construct
        mydebugWriter();

        //- Construct with some output options
        explicit mydebugWriter(const dictionary& options);

        //- Construct from components
        mydebugWriter
        (
            const meshedSurf& surf,
            const fileName& outputPath,
            bool parallel = UPstream::parRun(),
            const dictionary& options = dictionary()
        );

        //- Construct from components
        mydebugWriter
        (
            const pointField& points,
            const faceList& faces,
            const fileName& outputPath,
            bool parallel = UPstream::parRun(),
            const dictionary& options = dictionary()
        );


    //- Destructor
    virtual ~mydebugWriter() = default;


    // Member Functions

        //- Write surface geometry to file.
        virtual fileName write(); // override

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
