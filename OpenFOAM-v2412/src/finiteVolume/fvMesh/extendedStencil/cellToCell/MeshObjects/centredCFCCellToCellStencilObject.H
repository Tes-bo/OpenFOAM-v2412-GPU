/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::centredCFCCellToCellStencilObject

Description

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef centredCFCCellToCellStencilObject_H
#define centredCFCCellToCellStencilObject_H

#include "extendedCentredCellToCellStencil.H"
#include "CFCCellToCellStencil.H"
#include "MeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class centredCFCCellToCellStencilObject Declaration
\*---------------------------------------------------------------------------*/

class centredCFCCellToCellStencilObject
:
    public MeshObject
    <
        fvMesh,
        TopologicalMeshObject,
        centredCFCCellToCellStencilObject
    >,
    public extendedCentredCellToCellStencil
{
    // Private Typedefs

        typedef MeshObject
        <
            fvMesh,
            TopologicalMeshObject,
            centredCFCCellToCellStencilObject
        > MeshObject_type;

public:

    TypeName("centredCFCCellToCellStencil");

    // Constructors

        //- Construct from uncompacted cell stencil
        explicit centredCFCCellToCellStencilObject
        (
            const fvMesh& mesh
        )
        :
            MeshObject_type(mesh),
            extendedCentredCellToCellStencil(CFCCellToCellStencil(mesh))
        {}


    //- Destructor
    virtual ~centredCFCCellToCellStencilObject() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
