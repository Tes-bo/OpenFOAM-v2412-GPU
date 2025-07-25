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
    Foam::fileFormats::edgeMeshFormat

Description
    Provide a means of reading/writing the single-file OpenFOAM edge format.

Note
   This class provides more methods than the regular edge format interface.

SourceFiles
    edgeMeshFormat.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_edgeMeshFormat_H
#define Foam_edgeMeshFormat_H

#include "edgeMesh.H"
#include "Ostream.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                      Class edgeMeshFormat Declaration
\*---------------------------------------------------------------------------*/

class edgeMeshFormat
:
    public edgeMesh
{
protected:

    // Protected Member Functions

        //- Write header information
        static void writeHeader
        (
            Ostream&,
            const pointField&,
            const edgeList&
        );


public:

    // Constructors

        //- Construct from file name
        explicit edgeMeshFormat(const fileName& filename);


    // Selectors

        //- Read file and return edgeMesh
        static autoPtr<edgeMesh> New(const fileName& name)
        {
            return autoPtr<edgeMesh>::NewFrom<edgeMeshFormat>(name);
        }


    //- Destructor
    virtual ~edgeMeshFormat() = default;


    // Static Member Functions

        //- Read edgeMesh components from stream
        static bool read
        (
            Istream&,
            pointField&,
            edgeList&
        );

        //- Write edgeMesh components to stream
        static Ostream& write
        (
            Ostream&,
            const pointField&,
            const edgeList&
        );

        //- Write edge mesh to file with a mimicked IOobject header
        static void write
        (
            const fileName& filename,
            const edgeMesh& mesh,
            IOstreamOption streamOpt = IOstreamOption(),
            const dictionary& options = dictionary::null
        );


    // Member Functions

        //- Read from file
        virtual bool read(const fileName& filename) override;

        //- Write edge mesh to file
        virtual void write
        (
            const fileName& name,
            IOstreamOption streamOpt = IOstreamOption(),
            const dictionary& options = dictionary::null
        ) const override
        {
            write(name, *this, streamOpt, options);
        }

        //- Write edge mesh to file
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
