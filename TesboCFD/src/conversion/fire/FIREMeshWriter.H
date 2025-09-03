/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::fileFormats::FIREMeshWriter

Description
    Writes polyMesh in AVL/FIRE polyhedra format (fpma, fpmb)

    It is also possible to write compressed formats (fpmaz, fpmbz)

Note
   The fpma, fpmb formats are relatively poorly documented, but are manageable
   to read and write. It is, however, not recommended to import them directly
   into AVL/FIRE (the GUI) since it is generally not robust enough.
   Instead use their file-convertor to reconvert them into their native format.

   In the AVL/FIRE polyhedra format, the faces normals point inwards, whereas
   the OpenFOAM face normals always point outwards.

SourceFiles
    FIREMeshWriter.C

\*---------------------------------------------------------------------------*/

#ifndef FIREMeshWriter_H
#define FIREMeshWriter_H

#include "meshWriter.H"
#include "FIRECore.H"
#include "IOstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                 Class fileFormats::FIREMeshWriter Declaration
\*---------------------------------------------------------------------------*/

class FIREMeshWriter
:
    public meshWriter,
    public FIRECore
{
    // Private Member Functions

        //- No copy construct
        FIREMeshWriter(const FIREMeshWriter&) = delete;

        //- No copy assignment
        void operator=(const FIREMeshWriter&) = delete;

        //- Write points, faces, cells
        bool writeGeometry(OSstream&) const;

        //- Write selections
        bool writeSelections(OSstream&) const;

public:

    // Static data members

        //- Write binary (default ascii)
        static bool binary;

        //- Write with compression (default false)
        static bool compress;

        //- Prefix patches with 'BND_' before writing (default true)
        static bool prefixBoundary;


    // Constructors

        //- Prepare for writing, optionally with scaling.
        //  Treats a zero or negative scale factor as unity scaling.
        FIREMeshWriter(const polyMesh&, const scalar scaleFactor = 1.0);


    //- Destructor
    virtual ~FIREMeshWriter() = default;


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
