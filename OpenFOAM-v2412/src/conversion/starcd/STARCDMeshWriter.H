/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2016 OpenCFD Ltd.
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
    Foam::fileFormats::STARCDMeshWriter

Description
    Writes polyMesh in PROSTAR (v4) bnd/cel/vrt format

    The cellTableId and cellTable information are used (if available).
    Otherwise the cellZones are used (if available).

SourceFiles
    STARCDMeshWriter.C

\*---------------------------------------------------------------------------*/

#ifndef STARCDMeshWriter_H
#define STARCDMeshWriter_H

#include "meshWriter.H"
#include "STARCDCore.H"
#include "IOstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                Class fileFormats::STARCDMeshWriter Declaration
\*---------------------------------------------------------------------------*/

class STARCDMeshWriter
:
    public meshWriter,
    protected STARCDCore
{
    // Private data

        //- Write boundary (bnd) file - default true
        bool writeBoundary_;


    // Private Member Functions

        //- Write cells
        void writeCells(const fileName& baseName) const;

        //- Write boundary
        void writeBoundary(const fileName& baseName) const;

        void getCellTable();

        label findDefaultBoundary() const;


        //- No copy construct
        STARCDMeshWriter(const STARCDMeshWriter&) = delete;

        //- No copy assignment
        void operator=(const STARCDMeshWriter&) = delete;


public:

    // Constructors

        //- Write mesh files in PROSTAR format
        STARCDMeshWriter
        (
            const polyMesh&,
            const scalar scaleFactor = 1.0,
            const bool writeBndFile = true
        );


    //- Destructor
    virtual ~STARCDMeshWriter() = default;


    // Member Functions

    // Write

        //- Write volume mesh
        virtual bool write
        (
            const fileName& meshName = fileName::null
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
