/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::surfaceWriters::boundaryDataWriter

Description
    A surfaceWriter for outputting to a form usable for the
    timeVaryingMapped boundary condition. This reads the data from
    constant/boundaryData/\<patch\> directory.

    \verbatim
    formatOptions
    {
        boundaryData
        {
            header      true;
            format      ascii;
            compression false;
            normal      false;
        }
    }
    \endverbatim

    Format options:
    \table
        Property    | Description                           | Required | Default
        header      | Generate files with FoamFile header   | no  | true
        format      | ascii/binary                          | no  | ascii
        compression | Use file compression                  | no  | false
        scale       | Output geometry scaling               | no  | 1
        transform   | Output coordinate transform           | no  |
        fieldLevel  | Subtract field level before scaling   | no  | empty dict
        fieldScale  | Output field scaling                  | no  | empty dict
        normal      | Write face area normal in output      | no  | false
    \endtable

    Typical way of working:
    - use a sampledSurface of type 'patch' (to sample a patch):
    \verbatim
    surfaces
    {
        type            surfaces;
        fields          ( p );
        surfaceFormat   boundaryData;
        formatOptions
        {
            boundaryData
            {
                format  binary;
                normal  yes;
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
        surfaces
        {
            outlet
            {
                type            patch;
                patches         (outlet);
                interpolate     false;
            }
        }
    }
    \endverbatim

    - write using this writer.
    - move postProcessing/surfaces/outlet to constant/boundaryData/outlet
      in your destination case.
    - use a timeVaryingMappedFixedValue condition to read and interpolate
      the profile:
        type            timeVaryingMappedFixedValue;
        setAverage      false;  // do not use read average
        offset          0;      // do not apply offset to values

    Note:
    - with 'interpolate false' the data is on the face centres of the
      patch. Take care that a 2D geometry will only have a single row
      of face centres so might not provide a valid triangulation
      (this is what timeVaryingMappedFixedValue uses to do interpolation)
      (Alternatively use timeVaryingMappedFixedValue with mapMethod 'nearest')

    \heading Output file locations

    The \c rootdir normally corresponds to something like
    \c postProcessing/\<name\>

    where the geometry is written as:
    \verbatim
    rootdir
    `-- surfaceName
        `-- "points"
    \endverbatim

    and field data:
    \verbatim
    rootdir
    `-- surfaceName
        |-- "points"
        `-- timeName
            `-- field
    \endverbatim

SourceFiles
    boundaryDataSurfaceWriter.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_surfaceWriters_boundaryDataWriter_H
#define Foam_surfaceWriters_boundaryDataWriter_H

#include "surfaceWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class regIOobject;

namespace surfaceWriters
{

/*---------------------------------------------------------------------------*\
                     Class boundaryDataWriter Declaration
\*---------------------------------------------------------------------------*/

class boundaryDataWriter
:
    public surfaceWriter
{
    // Private Data

        //- Output stream option
        IOstreamOption streamOpt_;

        //- Output files with FoamFile header
        bool header_;

        //- Output face area normal
        const bool writeNormal_;


    // Private Member Functions

        //- Write serial surface geometry to "points" file.
        void serialWriteGeometry(const regIOobject&, const meshedSurf& surf);

        //- Templated write field operation
        template<class Type>
        fileName writeTemplate
        (
            const word& fieldName,          //!< Name of field
            const Field<Type>& localValues  //!< Local field values to write
        );


public:

    //- Declare type-name, virtual type (without debug switch)
    TypeNameNoDebug("boundaryData");


    // Constructors

        //- Default construct
        boundaryDataWriter();

        //- Construct with some output options
        explicit boundaryDataWriter(const dictionary& options);

        //- Construct from components
        boundaryDataWriter
        (
            const meshedSurf& surf,
            const fileName& outputPath,
            bool parallel = UPstream::parRun(),
            const dictionary& options = dictionary()
        );

        //- Construct from components
        boundaryDataWriter
        (
            const pointField& points,
            const faceList& faces,
            const fileName& outputPath,
            bool parallel = UPstream::parRun(),
            const dictionary& options = dictionary()
        );


    //- Destructor
    virtual ~boundaryDataWriter() = default;


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
