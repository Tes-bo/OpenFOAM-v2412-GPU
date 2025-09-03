/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
    Copyright (C) 2016,2024 OpenCFD Ltd.
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
    Foam::processorMeshes

Description
    Container for processor mesh addressing.

SourceFiles
    processorMeshes.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_processorMeshes_H
#define Foam_processorMeshes_H

#include "PtrList.H"
#include "fvMesh.H"
#include "IOobjectList.H"
#include "labelIOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class processorMeshes Declaration
\*---------------------------------------------------------------------------*/

class processorMeshes
{
    // Private Data

        const word meshName_;

        //- Processor databases
        PtrList<Time>& databases_;

        //- List of processor meshes
        PtrList<fvMesh> meshes_;

        //- List of processor point addressing lists
        PtrList<labelIOList> pointProcAddressing_;

        //- List of processor face addressing lists
        PtrList<labelIOList> faceProcAddressing_;

        //- List of processor cell addressing lists
        PtrList<labelIOList> cellProcAddressing_;

        //- List of processor boundary addressing lists
        PtrList<labelIOList> boundaryProcAddressing_;


        // pointMesh

            //- List of processor meshes
            PtrList<pointMesh> pMeshes_;

            //- List of processor boundary addressing lists for pointMeshes
            PtrList<labelIOList> pBoundaryProcAddressing_;


    // Private Member Functions

        //- Read all meshes
        void read();

        //- No copy construct
        processorMeshes(const processorMeshes&) = delete;

        //- No copy assignment
        void operator=(const processorMeshes&) = delete;


public:

    //- Runtime type information
    ClassName("processorMeshes");


    // Constructors

        //- Construct from components
        processorMeshes(PtrList<Time>& databases, const word& meshName);


    // Member Functions

        //- Update the meshes based on the mesh files saved in time directories
        polyMesh::readUpdateState readUpdate();

        //- Reconstruct point position after motion in parallel
        void reconstructPoints(fvMesh&);

        const PtrList<fvMesh>& meshes() const noexcept
        {
            return meshes_;
        }

        PtrList<fvMesh>& meshes() noexcept
        {
            return meshes_;
        }

        const PtrList<labelIOList>& pointProcAddressing() const noexcept
        {
            return pointProcAddressing_;
        }

        PtrList<labelIOList>& faceProcAddressing() noexcept
        {
            return faceProcAddressing_;
        }

        const PtrList<labelIOList>& cellProcAddressing() const noexcept
        {
            return cellProcAddressing_;
        }

        const PtrList<labelIOList>& boundaryProcAddressing() const noexcept
        {
            return boundaryProcAddressing_;
        }


        // pointMesh

            const PtrList<pointMesh>& pointMeshes() const noexcept
            {
                return pMeshes_;
            }

            const PtrList<labelIOList>& pointMeshBoundaryProcAddressing()
            const noexcept
            {
                return pBoundaryProcAddressing_;
            }


        //- Helper: remove all procAddressing files from mesh instance
        static void removeFiles(const polyMesh& mesh);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
