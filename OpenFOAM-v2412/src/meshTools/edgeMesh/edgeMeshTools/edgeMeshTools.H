/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenCFD Ltd.
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

Namespace
    Foam::edgeMeshTools

Description
    Collection of static functions related to edgeMesh features.

SourceFiles
    edgeMeshTools.C
    edgeMeshFeatureProximity.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_edgeMeshTools_H
#define Foam_edgeMeshTools_H

#include "tmp.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class edgeMesh;
class extendedEdgeMesh;
class extendedFeatureEdgeMesh;
class triSurface;
class Time;

/*---------------------------------------------------------------------------*\
                      Namespace edgeMeshTools Declaration
\*---------------------------------------------------------------------------*/

namespace edgeMeshTools
{
    //- Write some information
    void writeStats
    (
        Ostream& os,
        const extendedFeatureEdgeMesh& emesh
    );


    //- Calculate proximity of the features to the surface
    tmp<scalarField> featureProximity
    (
        const extendedEdgeMesh& emesh,
        const triSurface& surf,
        const scalar searchDistance
    );


    //- Calculate proximity of the features to the surface and write the field
    tmp<scalarField> writeFeatureProximity
    (
        const Time& runTime,
        const word& basename,
        const extendedEdgeMesh& emesh,
        const triSurface& surf,
        const scalar searchDistance
    );

} // End namespace edgeMeshTools


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
