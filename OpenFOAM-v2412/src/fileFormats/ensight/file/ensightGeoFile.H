/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::ensightGeoFile

Description
    A variant of ensightFile (Ensight writing) that includes
    the extra geometry file header information.

\*---------------------------------------------------------------------------*/

#ifndef Foam_ensightGeoFile_H
#define Foam_ensightGeoFile_H

#include "ensightFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class ensightGeoFile Declaration
\*---------------------------------------------------------------------------*/

class ensightGeoFile
:
    public ensightFile
{
public:

    // Static Functions

        //- Return a null ensightGeoFile
        static const ensightGeoFile& null() noexcept
        {
            return NullObjectRef<ensightGeoFile>();
        }


    // Generated Methods

        //- No copy construct
        ensightGeoFile(const ensightGeoFile&) = delete;

        //- No copy assignment
        void operator=(const ensightGeoFile&) = delete;


    // Constructors

        //- Construct from path-name and initialise the binary header.
        //  The path-name is adjusted for valid ensight file naming.
        //  Created as an atomic or in append mode (single-file format).
        ensightGeoFile
        (
            //! Append mode: NO_APPEND or APPEND_ATE
            IOstreamOption::appendType append,
            const fileName& pathname,
            IOstreamOption::streamFormat fmt = IOstreamOption::BINARY
        );

        //- Construct from path / name and initialise the binary header.
        //  Only the name portion is adjusted for valid ensight file naming.
        //  Created as an atomic or in append mode (single-file format).
        ensightGeoFile
        (
            //! Append mode: NO_APPEND or APPEND_ATE
            IOstreamOption::appendType append,
            const fileName& path,
            const fileName& name,
            IOstreamOption::streamFormat fmt = IOstreamOption::BINARY
        );

        //- Construct from path-name and initialise the binary header.
        //  The path-name is adjusted for valid ensight file naming.
        //  Created as an atomic, non-append mode.
        //  \note Since 2024-05 does \em not call beginGeometry() !!
        explicit ensightGeoFile
        (
            const fileName& pathname,
            IOstreamOption::streamFormat fmt = IOstreamOption::BINARY
        )
        :
            ensightGeoFile(IOstreamOption::NO_APPEND, pathname, fmt)
        {}

        //- Construct from path / name and initialise the binary header.
        //  Only the name portion is adjusted for valid ensight file naming.
        //  Created as an atomic, non-append mode.
        //  \note Since 2024-05 does \em not call beginGeometry() !!
        ensightGeoFile
        (
            const fileName& path,
            const fileName& name,
            IOstreamOption::streamFormat fmt = IOstreamOption::BINARY
        )
        :
            ensightGeoFile(IOstreamOption::NO_APPEND, path, name, fmt)
        {}


    //- Destructor
    ~ensightGeoFile() = default;


    // Member Functions

    // Output

        //- Write keyword with trailing newline
        virtual Ostream& writeKeyword(const keyType& key) override;


    // Convenience Output Methods

        //- Start of geometry information
        void beginGeometry();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
