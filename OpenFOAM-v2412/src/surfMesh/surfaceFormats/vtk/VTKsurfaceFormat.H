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
    Foam::fileFormats::VTKsurfaceFormat

Description
    Read/write VTK legacy format (ASCII) for surfaces.

    The output is never sorted by zone.

    Output stream options: ignored (dictionary options only)

    \heading Format options
    \table
        Property    | Description                       | Required | Default
        format      | ascii or binary format            | no | ascii
        precision   | Write precision in ascii          | no | same as IOstream
    \endtable

SourceFiles
    VTKsurfaceFormat.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_VTKsurfaceFormat_H
#define Foam_VTKsurfaceFormat_H

#include "MeshedSurface.H"
#include "MeshedSurfaceProxy.H"
#include "UnsortedMeshedSurface.H"
#include "VTKsurfaceFormatCore.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                Class fileFormats::VTKsurfaceFormat Declaration
\*---------------------------------------------------------------------------*/

template<class Face>
class VTKsurfaceFormat
:
    public MeshedSurface<Face>,
    public VTKsurfaceFormatCore
{
    // Private Member Functions

        //- Write polygons
        static void writePolys
        (
            vtk::formatter& format,
            const UList<Face>& faces
        );


public:

    // Constructors

        //- Construct from file name
        explicit VTKsurfaceFormat(const fileName& filename);


    //- Destructor
    virtual ~VTKsurfaceFormat() = default;


    // Static Functions

        //- Write surface mesh components (by proxy) in legacy VTK format
        static void write
        (
            const fileName& filename,
            const MeshedSurfaceProxy<Face>& surf,
            IOstreamOption /*ignored*/ = IOstreamOption(),
            const dictionary& options = dictionary::null
        );

        //- Write UnsortedMeshedSurface, the output remains unsorted
        static void write
        (
            const fileName& fileName,
            const UnsortedMeshedSurface<Face>& surf,
            IOstreamOption /*ignored*/ = IOstreamOption(),
            const dictionary& options = dictionary::null
        );


    // Member Functions

        //- Read from file
        virtual bool read(const fileName& filename) override;

        //- Write surface mesh to file (by proxy)
        virtual void write
        (
            const fileName& name,
            IOstreamOption streamOpt = IOstreamOption(),
            const dictionary& options = dictionary::null
        ) const override
        {
            write(name, MeshedSurfaceProxy<Face>(*this), streamOpt, options);
        }

        //- Write surface mesh to file (by proxy)
        virtual void write
        (
            const fileName& name,
            const word& fileType,  /* ignored */
            IOstreamOption streamOpt = IOstreamOption(),
            const dictionary& options = dictionary::null
        ) const override
        {
            write(name, MeshedSurfaceProxy<Face>(*this), streamOpt, options);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "VTKsurfaceFormat.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
