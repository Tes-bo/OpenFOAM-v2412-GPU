/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::pureUpwindCFCCellToFaceStencilObject

Description

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef pureUpwindCFCCellToFaceStencilObject_H
#define pureUpwindCFCCellToFaceStencilObject_H

#include "extendedUpwindCellToFaceStencil.H"
#include "CFCCellToFaceStencil.H"
#include "MeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class pureUpwindCFCCellToFaceStencilObject Declaration
\*---------------------------------------------------------------------------*/

class pureUpwindCFCCellToFaceStencilObject
:
    public MeshObject
    <
        fvMesh,
        TopologicalMeshObject,
        pureUpwindCFCCellToFaceStencilObject
    >,
    public extendedUpwindCellToFaceStencil
{
    // Private Typedefs

        typedef MeshObject
        <
            fvMesh,
            TopologicalMeshObject,
            pureUpwindCFCCellToFaceStencilObject
        > MeshObject_type;

public:

    TypeName("pureUpwindCFCCellToFaceStencil");


    // Constructors

        //- Construct from uncompacted face stencil
        explicit pureUpwindCFCCellToFaceStencilObject
        (
            const fvMesh& mesh
        )
        :
            MeshObject_type(mesh),
            extendedUpwindCellToFaceStencil(CFCCellToFaceStencil(mesh))
        {
            if (extendedCellToFaceStencil::debug)
            {
                Info<< "Generated pure upwind stencil " << type()
                    << nl << endl;
                writeStencilStats(Info, ownStencil(), ownMap());
            }
        }


    //- Destructor
    virtual ~pureUpwindCFCCellToFaceStencilObject() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
