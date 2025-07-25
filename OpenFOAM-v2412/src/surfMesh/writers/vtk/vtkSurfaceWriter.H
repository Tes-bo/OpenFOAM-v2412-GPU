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
    Foam::surfaceWriters::vtkWriter

Description
    A surfaceWriter for VTK legacy (.vtk) or XML (.vtp) format.

    The formatOptions for vtk:
    \table
        Property    | Description                           | Required | Default
        format      | ascii or binary format                | no  | binary
        legacy      | Legacy VTK output                     | no  | false
        precision   | Write precision in ascii         | no | same as IOstream
        scale       | Output geometry scaling               | no  | 1
        transform   | Output coordinate transform           | no  |
        fieldLevel  | Subtract field level before scaling   | no  | empty dict
        fieldScale  | Output field scaling                  | no  | empty dict
        normal      | Write face area-normal in output      | no  | false
    \endtable

    For example,
    \verbatim
    formatOptions
    {
        vtk
        {
            format      binary;
            legacy      false;
            precision   10;
            fieldLevel
            {
                p       1e5;    // Absolute -> gauge [Pa]
            }
            fieldScale
            {
                "p.*"   0.01;   // [Pa] -> [mbar]
            }
        }
    }
    \endverbatim

    \section Output file locations

    The \c rootdir normally corresponds to something like
    \c postProcessing/\<name\>

    \subsection Geometry and Fields
    \verbatim
    rootdir
    `-- timeName
        `-- surfaceName.{vtk,vtp}
    \endverbatim

SourceFiles
    vtkSurfaceWriter.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_surfaceWriters_vtkWriter_H
#define Foam_surfaceWriters_vtkWriter_H

#include "surfaceWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace vtk
{
// Forward Declarations
class outputOptions;
class surfaceWriter;
}

namespace surfaceWriters
{

/*---------------------------------------------------------------------------*\
                          Class vtkWriter Declaration
\*---------------------------------------------------------------------------*/

class vtkWriter
:
    public surfaceWriter
{
    // Private Data

        //- The VTK output format type.
        //  Stored as a raw value to avoid a header dependency on fileFormats
        unsigned fmtType_;

        //- ASCII write precision
        unsigned precision_;

        //- Output face area normal
        const bool writeNormal_;

        //- Backend writer - master only
        autoPtr<Foam::vtk::surfaceWriter> writer_;


    // Private Member Functions

        //- Templated write field operation
        template<class Type>
        fileName writeTemplate
        (
            const word& fieldName,          //!< Name of field
            const Field<Type>& localValues  //!< Local field values to write
        );


public:

    //- Declare type-name, virtual type (without debug switch)
    TypeNameNoDebug("vtk");


    // Constructors

        //- Default construct
        vtkWriter();

        //- Construct with some output options
        explicit vtkWriter(const vtk::outputOptions& opts);

        //- Construct with some output options
        explicit vtkWriter(const dictionary& options);

        //- Construct from components
        //  The file name is with/without an extension.
        vtkWriter
        (
            const meshedSurf& surf,
            const fileName& outputPath,
            bool parallel = UPstream::parRun(),
            const dictionary& options = dictionary()
        );

        //- Construct from components with specified output path.
        //  The file name is with/without an extension.
        vtkWriter
        (
            const pointField& points,
            const faceList& faces,
            const fileName& outputPath,
            bool parallel = UPstream::parRun(),
            const dictionary& options = dictionary()
        );


    //- Destructor. Calls close()
    virtual ~vtkWriter();


    // Member Functions

        //- Finish output, clears backend.
        virtual void close(); // override

        //- Begin time step. Clears existing backend.
        virtual void beginTime(const Time& t); // override

        //- Begin time step. Clears existing backend.
        virtual void beginTime(const instant& inst); // override

        //- End time step. Clears existing backend.
        virtual void endTime(); // override


    // Write

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
