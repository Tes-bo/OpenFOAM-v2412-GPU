/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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
    Foam::surfaceWriters::starcdWriter

Description
    A surfaceWriter for STARCD files.

    The formatOptions for starcd:
    \table
        Property    | Description                           | Required | Default
        compression | Use file compression                  | no  | false
        scale       | Output geometry scaling               | no  | 1
        transform   | Output coordinate transform           | no  |
        fieldLevel  | Subtract field level before scaling   | no  | empty dict
        fieldScale  | Output field scaling                  | no  | empty dict
    \endtable

    The geometry is written via the MeshedSurfaceProxy, the fields
    are written in a trivial ASCII format with ID and VALUE as
    so-called user data. These \c .usr files can be read into proSTAR
    with these types of commands. For element data:
    \verbatim
        getuser FILENAME.usr cell scalar free
        getuser FILENAME.usr cell vector free
    \endverbatim
    and for vertex data:
    \verbatim
        getuser FILENAME.usr vertex scalar free
        getuser FILENAME.usr vertex vector free
    \endverbatim

    \section Output file locations

    The \c rootdir normally corresponds to something like
    \c postProcessing/\<name\>

    \subsection Geometry
    \verbatim
    rootdir
    `-- timeName
        `-- surfaceName.{cel,vrt,inp}
    \endverbatim

    \subsection Fields
    \verbatim
    rootdir
    `-- timeName
        |-- <field0>_surfaceName.{usr}
        `-- <field1>_surfaceName.{usr}
    \endverbatim

SourceFiles
    starcdSurfaceWriter.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_surfaceWriters_starcdWriter_H
#define Foam_surfaceWriters_starcdWriter_H

#include "surfaceWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceWriters
{

/*---------------------------------------------------------------------------*\
                        Class starcdWriter Declaration
\*---------------------------------------------------------------------------*/

class starcdWriter
:
    public surfaceWriter
{
    // Private Data

        //- Output stream option
        IOstreamOption streamOpt_;


    // Private Member Functions

        //- Templated write operation
        template<class Type>
        fileName writeTemplate
        (
            const word& fieldName,          //!< Name of field
            const Field<Type>& localValues  //!< Local field values to write
        );


public:

    //- Declare type-name, virtual type (without debug switch)
    TypeNameNoDebug("starcd");


    // Constructors

        //- Default construct
        starcdWriter();

        //- Construct with some output options
        explicit starcdWriter(const dictionary& options);

        //- Construct from components
        starcdWriter
        (
            const meshedSurf& surf,
            const fileName& outputPath,
            bool parallel = UPstream::parRun(),
            const dictionary& options = dictionary()
        );

        //- Construct from components
        starcdWriter
        (
            const pointField& points,
            const faceList& faces,
            const fileName& outputPath,
            bool parallel = UPstream::parRun(),
            const dictionary& options = dictionary()
        );


    //- Destructor
    virtual ~starcdWriter() = default;


    // Member Functions

        //- True if the surface format supports geometry in a separate file.
        //  False if geometry and field must be in a single file
        virtual bool separateGeometry() const
        {
            return true;
        }

        //- Format uses faceIds as part of its output
        virtual bool usesFaceIds() const // override
        {
            return true;
        }

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
