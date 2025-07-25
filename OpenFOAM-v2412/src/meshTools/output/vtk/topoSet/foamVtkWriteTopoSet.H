/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2018 OpenCFD Ltd.
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

InNamespace
    Foam::vtk

Description
    Write topoSet in VTK format

SourceFiles
    foamVtkWriteTopoSet.C
    foamVtkWriteCellSetFaces.C
    foamVtkWriteFaceSet.C
    foamVtkWritePointSet.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_vtk_WriteTopoSet_H
#define Foam_vtk_WriteTopoSet_H

#include "fileName.H"
#include "UPstream.H"
#include "foamVtkOutputOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class polyMesh;
class topoSet;
class faceSet;
class cellSet;
class pointSet;

namespace vtk
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Dispatch to vtk::writeCellSetFaces, vtk::writeFaceSet, vtk::writePointSet
//  The file name is with/without an extension.
//  \return True on successful dispatch
bool writeTopoSet
(
    const polyMesh& mesh,
    const topoSet& set,
    const vtk::outputOptions opts,
    const fileName& file,
    bool parallel = UPstream::parRun()
);


//- Write faceSet as VTK polydata file.
//  Only one CELL_DATA, which is the original faceID
//  The file name is with/without an extension.
bool writeFaceSet
(
    const polyMesh& mesh,
    const faceSet& set,
    const vtk::outputOptions opts,
    const fileName& file,
    bool parallel = UPstream::parRun()
);


//- Write perimeter faces of cellset to vtk polydata file.
//  The data are the original cell ids
//  The file name is with/without an extension.
bool writeCellSetFaces
(
    const polyMesh& mesh,
    const cellSet& set,
    const vtk::outputOptions opts,
    const fileName& file,
    bool parallel = UPstream::parRun()
);


//- Write pointSet to VTK polydata file.
//  Only one CELL_DATA, which is the original pointID
//  The file name is with/without an extension.
bool writePointSet
(
    const polyMesh& mesh,
    const pointSet& set,
    const vtk::outputOptions opts,
    const fileName& file,
    bool parallel = UPstream::parRun()
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
