/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2015-2023 OpenCFD Ltd.
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
    Foam::surfaceWriters::ensightWriter

Description
    A surfaceWriter for Ensight format.

    \verbatim
    formatOptions
    {
        ensight
        {
            format          ascii;
            collateTimes    true;
        }
    }
    \endverbatim

    Format options for ensight:
    \table
        Property | Description                              | Required | Default
        format   | ascii/binary                             | no  | binary
        collateTimes | Use common geometry for times        | no  | true
        scale       | Output geometry scaling               | no  | 1
        transform   | Output coordinate transform           | no  |
        fieldLevel  | Subtract field level before scaling   | no  | empty dict
        fieldScale  | Output field scaling                  | no  | empty dict
        timeFormat  | Time format (ensight case)            | no  | scientific
        timePrecision | Time precision (ensight case)       | no  | 5
    \endtable

    The collated format maintains an internal list of the known times
    as well as a file-cached version with the field information.
    The information is used for restarts.

SourceFiles
    ensightSurfaceWriter.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_surfaceWriters_ensightWriter_H
#define Foam_surfaceWriters_ensightWriter_H

#include "surfaceWriter.H"
#include "ensightCase.H"
#include "ensightWriterCaching.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceWriters
{

/*---------------------------------------------------------------------------*\
                        Class ensightWriter Declaration
\*---------------------------------------------------------------------------*/

class ensightWriter
:
    public surfaceWriter
{
    // Private Data

        //- Ensight case options
        ensightCase::options caseOpts_;

        //- Collate times (default: true)
        bool collateTimes_;

        //- Cached information for times, geometry, fields (collated)
        ensightOutput::writerCaching caching_;


    // Private Member Functions

        //- Write geometry
        fileName writeCollated();

        //- Write geometry
        fileName writeUncollated();

        //- Templated write operation - one file per timestep
        template<class Type>
        fileName writeCollated
        (
            const word& fieldName,          //!< Name of field
            const Field<Type>& localValues  //!< Local field values to write
        );

        //- Templated write operation - all time steps in single file
        template<class Type>
        fileName writeUncollated
        (
            const word& fieldName,          //!< Name of field
            const Field<Type>& localValues  //!< Local field values to write
        );

        //- Templated write operation
        template<class Type>
        fileName writeTemplate
        (
            const word& fieldName,          //!< Name of field
            const Field<Type>& localValues  //!< Local field values to write
        );

public:

    //- Declare type-name, virtual type (without debug switch)
    TypeNameNoDebug("ensight");


    // Constructors

        //- Default construct
        ensightWriter();

        //- Construct with some output options
        explicit ensightWriter(const dictionary& options);

        //- Construct from components
        ensightWriter
        (
            const meshedSurf& surf,
            const fileName& outputPath,
            bool parallel = UPstream::parRun(),
            const dictionary& options = dictionary()
        );

        //- Construct from components
        ensightWriter
        (
            const pointField& points,
            const faceList& faces,
            const fileName& outputPath,
            bool parallel = UPstream::parRun(),
            const dictionary& options = dictionary()
        );


    //- Destructor
    virtual ~ensightWriter() = default;


    // Member Functions

        //- Finish output, clears output times.
        //  Later reuse will rebuild times from fieldsDict file cache.
        virtual void close(); // override

        //- True if the surface format supports geometry in a separate file.
        //  False if geometry and field must be in a single file
        virtual bool separateGeometry() const;

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
