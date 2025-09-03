/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2018 OpenCFD Ltd.
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
    Foam::fileFormats::FIREMeshReader

Description
    Read AVL/FIRE fpma, fpmb files.

Note
    Does not handle compressed versions (fpmaz, fpmbz) of these files.

SourceFiles
    FIREMeshReader.C

\*---------------------------------------------------------------------------*/

#ifndef FIREMeshReader_H
#define FIREMeshReader_H

#include "meshReader.H"
#include "FIRECore.H"
#include "labelList.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// Forward declarations
class polyMesh;

namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                 Class fileFormats::FIREMeshReader Declaration
\*---------------------------------------------------------------------------*/

class FIREMeshReader
:
    public meshReader,
    public FIRECore
{

protected:

    // Protected Data

        labelList  owner_;
        labelList  neigh_;

        labelList  faceZoneId_;
        wordList   faceNames_;


    // Protected Member Functions

        //- No copy construct
        FIREMeshReader(const FIREMeshReader&) = delete;

        //- No copy assignment
        void operator=(const FIREMeshReader&) = delete;


        //- Read the mesh from the file(s)
        virtual bool readGeometry(const scalar scaleFactor = 1.0);

        //- Read points from file
        void readPoints(ISstream&, const scalar scaleFactor = 1.0);

        //- Read points from file
        void readFaces(ISstream&);

        //- Read cell connectivities from file
        void readCells(ISstream&);

        //- Read cell/face selections from file
        void readSelections(ISstream&);

        //-
        void reorganize();

        void addPatches(polyMesh&) const;


public:

    // Constructors

        //- Construct by reading file, optionally with scaling
        FIREMeshReader(const fileName&, const scalar scaleFactor = 1.0);


    //- Destructor
    virtual ~FIREMeshReader() = default;


    // Member Functions

        //- Create and return polyMesh
        virtual autoPtr<polyMesh> mesh(const objectRegistry&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
