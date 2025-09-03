/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::fileFormats::STARCDedgeFormat

Description
    Read/write the lines from PROSTAR vrt/cel files.

Note
    Uses the extension \a .inp (input) to denote the format.

See also
    Foam::meshReaders::STARCD

SourceFiles
    STARCDedgeFormat.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_STARCDedgeFormat_H
#define Foam_STARCDedgeFormat_H

#include "edgeMesh.H"
#include "STARCDCore.H"

#include "Fstream.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                      Class STARCDedgeFormat Declaration
\*---------------------------------------------------------------------------*/

class STARCDedgeFormat
:
    public edgeMesh,
    protected STARCDCore
{
    // Private Member Functions

        static inline void writeLines
        (
            Ostream& os,
            const edgeList& edges,
            label starCellId = 1  // 1-based cellId
        );

protected:

    // Protected Member Functions

    static void writeCase
    (
        Ostream& os,
        const pointField& pointLst,
        const label nEdges
    );


public:

    // Constructors

        //- Construct from file name
        explicit STARCDedgeFormat(const fileName& filename);


    // Selectors

        //- Read file and return edgeMesh
        static autoPtr<edgeMesh> New(const fileName& name)
        {
            return autoPtr<edgeMesh>::NewFrom<STARCDedgeFormat>(name);
        }


    //- Destructor
    virtual ~STARCDedgeFormat() = default;


    // Static Functions

        //- Write edge mesh to file in STARCD format
        static void write
        (
            const fileName& filename,
            const edgeMesh& mesh,
            IOstreamOption streamOpt = IOstreamOption(),
            const dictionary& options = dictionary::null
        );


    // Member Functions

        //- Read from file
        virtual bool read(const fileName& name) override;

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
