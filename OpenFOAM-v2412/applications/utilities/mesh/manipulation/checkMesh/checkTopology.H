/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Description
    Global functions for topology checking

\*---------------------------------------------------------------------------*/

#include "labelList.H"
#include "autoPtr.H"
#include "ZoneMesh.H"
#include "topoSet.H"

namespace Foam
{

// Forward Declarations
class polyMesh;
class pointSet;
class coordSetWriter;
class surfaceWriter;


// Check patch topology.
// In parallel, uses count of edge connections
template<class PatchType>
bool checkPatch
(
    const bool allGeometry,
    const std::string& name,
    const polyMesh& mesh,
    const PatchType& pp,
    const labelUList& meshEdges,
    labelHashSet* pointSetPtr = nullptr,
    labelHashSet* badEdgesPtr = nullptr
);

// OLDER CODE
// ~~~~~~~~~~
// Check patch topology.
// In parallel, uses mesh face ids and global face numbering
template<class PatchType>
bool checkPatch
(
    const bool allGeometry,
    const std::string& name,
    const polyMesh& mesh,
    const PatchType& pp,
    const labelUList& meshFaces,
    const labelUList& meshEdges,
    labelHashSet* pointSetPtr = nullptr,
    labelHashSet* badEdgesPtr = nullptr
);

template<class Zone>
label checkZones
(
    const polyMesh& mesh,
    const ZoneMesh<Zone, polyMesh>& zones,
    topoSet& set
);

label checkTopology
(
    const polyMesh& mesh,
    const bool allTopology,
    const bool allGeometry,
    autoPtr<surfaceWriter>& surfWriter,
    autoPtr<coordSetWriter>& setWriter,
    const bool writeBadEdges = false
);

} // End namespace Foam

// ************************************************************************* //
