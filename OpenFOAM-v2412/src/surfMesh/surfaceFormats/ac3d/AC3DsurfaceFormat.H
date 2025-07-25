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
    Foam::fileFormats::AC3DsurfaceFormat

Description
    Read/write AC3D format.

    http://www.inivis.com/ac3d/man/ac3dfileformat.html

    Output stream options:
       - ASCII only
       - compression on/off

    Output dictionary options: ignored

Note
    On input, the faces are already organized as zones.
    The output is always sorted by zones.
    In the absence of zones, a single zone will be assigned.

SourceFiles
    AC3DsurfaceFormat.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_AC3DsurfaceFormat_H
#define Foam_AC3DsurfaceFormat_H

#include "MeshedSurface.H"
#include "MeshedSurfaceProxy.H"
#include "UnsortedMeshedSurface.H"
#include "AC3DsurfaceFormatCore.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
               Class fileFormats::AC3DsurfaceFormat Declaration
\*---------------------------------------------------------------------------*/

template<class Face>
class AC3DsurfaceFormat
:
    public MeshedSurface<Face>,
    public AC3DsurfaceFormatCore
{
public:

    // Constructors

        //- Construct from file name
        explicit AC3DsurfaceFormat(const fileName& filename);


    //- Destructor
    virtual ~AC3DsurfaceFormat() = default;


    // Static Functions

        //- Write surface mesh components by proxy
        static void write
        (
            const fileName& filename,
            const MeshedSurfaceProxy<Face>& surf,
            IOstreamOption streamOpt = IOstreamOption(),
            const dictionary& /*unused*/ = dictionary::null
        );

        //- Write UnsortedMeshedSurface, the output is always sorted by zones.
        static void write
        (
            const fileName& filename,
            const UnsortedMeshedSurface<Face>& surf,
            IOstreamOption streamOpt = IOstreamOption(),
            const dictionary& /*unused*/ = dictionary::null
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
    #include "AC3DsurfaceFormat.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
