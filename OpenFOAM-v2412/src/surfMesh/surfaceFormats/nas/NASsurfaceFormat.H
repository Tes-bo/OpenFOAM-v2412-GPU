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
    Foam::fileFormats::NASsurfaceFormat

Description
    Nastran surface reader/writer.

    - Uses the Ansa "$ANSA_NAME" or the Hypermesh "$HMNAME COMP" extensions
      to obtain zone names.
    - Handles Nastran short, long formats and comma-separated free format.
    - Properly handles the Nastran compact floating point notation: \n
    \verbatim
        GRID          28        10.20269-.030265-2.358-8
    \endverbatim

    The Nastran writer uses FREE format only.

    Output stream options:
       - ASCII only
       - compression on/off

    Output dictionary options: ignored

SourceFiles
    NASsurfaceFormat.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_NASsurfaceFormat_H
#define Foam_NASsurfaceFormat_H

#include "MeshedSurface.H"
#include "MeshedSurfaceProxy.H"
#include "UnsortedMeshedSurface.H"
#include "NASCore.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                Class fileFormats::NASsurfaceFormat Declaration
\*---------------------------------------------------------------------------*/

template<class Face>
class NASsurfaceFormat
:
    public MeshedSurface<Face>,
    public NASCore
{
    // Private Member Functions

        //- Output CTRIA3 or CQUAD4
        inline static label writeShell
        (
            Ostream& os,
            const Face& f,
            label elemId,           //!< 0-based element id
            const label groupId     //!< 0-based property/group id
        );


public:

    // Constructors

        //- Construct from file name
        explicit NASsurfaceFormat(const fileName& filename);


    //- Destructor
    virtual ~NASsurfaceFormat() = default;


    // Static Functions

        //- Write surface mesh components (by proxy) in NAS format
        static void write
        (
            const fileName& filename,
            const MeshedSurfaceProxy<Face>& surf,
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
    #include "NASsurfaceFormat.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
