/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::fileFormats::FLMAsurfaceFormat

Description
    Provide a means of writing AVL/FIRE FLMA format.

    Output stream options:
       - ASCII only
       - compression : ignored

    Output dictionary options: ignored

Note
   The internal logic is still a bit odd (2020-02) since compressed files
   must have a \b .flmaz ending.
   To be revisited in the future.

SourceFiles
    FLMAsurfaceFormat.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_FLMAsurfaceFormat_H
#define Foam_FLMAsurfaceFormat_H

#include "MeshedSurface.H"
#include "MeshedSurfaceProxy.H"
#include "FIRECore.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
               Class fileFormats::FLMAsurfaceFormat Declaration
\*---------------------------------------------------------------------------*/

template<class Face>
class FLMAsurfaceFormat
:
    public MeshedSurface<Face>,
    public FIRECore
{
    // Private Member Functions

        static inline void writeShell(OSstream& os, const Face& f);
        static inline void writeType(OSstream& os, const Face& f);


protected:

    // Protected Member Functions

        //- Write surface mesh components (by proxy)
        static void write
        (
            OSstream& os,
            const MeshedSurfaceProxy<Face>& surf
        );

        //- Write surface mesh components (by proxy) with/without compression
        static void write
        (
            IOstreamOption::compressionType comp,
            const fileName& filename,
            const MeshedSurfaceProxy<Face>& surf
        );


public:

    // Constructors

        //- Default construct
        FLMAsurfaceFormat() = default;


    //- Destructor
    virtual ~FLMAsurfaceFormat() = default;


    // Static Functions

        //- Write surface mesh components (by proxy)
        static void write
        (
            const fileName& filename,
            const MeshedSurfaceProxy<Face>& surf,
            IOstreamOption /*ignored*/ = IOstreamOption(),
            const dictionary& /*unused*/ = dictionary::null
        );


    // Member Functions

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


/*---------------------------------------------------------------------------*\
                     Class FLMAZsurfaceFormat Declaration
\*---------------------------------------------------------------------------*/

template<class Face>
class FLMAZsurfaceFormat
:
    public FLMAsurfaceFormat<Face>
{
public:

    // Constructors

        //- Default construct
        FLMAZsurfaceFormat() = default;


    //- Destructor
    virtual ~FLMAZsurfaceFormat() = default;


    // Static Functions

        //- Write surface mesh components (by proxy)
        static void write
        (
            const fileName& filename,
            const MeshedSurfaceProxy<Face>& surf,
            IOstreamOption /*ignored*/ = IOstreamOption(),
            const dictionary& /*unused*/ = dictionary::null
        );


    // Member Functions

        //- Write surface mesh as flmaz file (by proxy)
        virtual void write
        (
            const fileName& name,
            IOstreamOption streamOpt = IOstreamOption(),
            const dictionary& options = dictionary::null
        ) const override
        {
            write(name, MeshedSurfaceProxy<Face>(*this), streamOpt, options);
        }

        //- Write surface mesh as flmaz file (by proxy)
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
#   include "FLMAsurfaceFormat.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
