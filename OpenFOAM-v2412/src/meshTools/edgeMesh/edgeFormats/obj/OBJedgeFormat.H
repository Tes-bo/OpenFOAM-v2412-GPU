/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2021-2024 OpenCFD Ltd.
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
    Foam::fileFormats::OBJedgeFormat

Description
    Provide a means of reading/writing Alias/Wavefront OBJ format.

    Does not handle negative vertex indices.

SourceFiles
    OBJedgeFormat.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_OBJedgeFormat_H
#define Foam_OBJedgeFormat_H

#include "edgeMesh.H"
#include "Fstream.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                        Class OBJedgeFormat Declaration
\*---------------------------------------------------------------------------*/

class OBJedgeFormat
:
    public edgeMesh
{
public:

    // Constructors

        //- Construct from file name
        explicit OBJedgeFormat(const fileName& filename);


    // Selectors

        //- Read file and return surface
        static autoPtr<edgeMesh> New(const fileName& name)
        {
            return autoPtr<edgeMesh>::NewFrom<OBJedgeFormat>(name);
        }


    //- Destructor
    virtual ~OBJedgeFormat() = default;


    // Static Functions

        //- Write edge mesh to file in OBJ format
        static void write
        (
            const fileName& filename,
            const edgeMesh& mesh,
            IOstreamOption /*ignored*/ = IOstreamOption(),
            const dictionary& options = dictionary::null
        );


    // Member Functions

        //- Read from file
        virtual bool read(const fileName& filename) override;

        //- Write to file
        virtual void write
        (
            const fileName& name,
            IOstreamOption streamOpt = IOstreamOption(),
            const dictionary& options = dictionary::null
        ) const override
        {
            write(name, *this, streamOpt, options);
        }

        //- Write to file
        virtual void write
        (
            const fileName& name,
            const word& fileType,  /* ignored */
            IOstreamOption streamOpt = IOstreamOption(),
            const dictionary& options = dictionary::null
        ) const override
        {
            write(name, *this, streamOpt, options);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
