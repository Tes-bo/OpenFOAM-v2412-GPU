/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::faMeshBoundaryHalo

Description
    Class for obtaining halo face data for the boundary edges.
    The ordering follows that natural edge ordering of the underlying
    primitive patch.

Note
    The halo faces can be located on-processor or off-processor.

SourceFiles
    faMeshBoundaryHalo.C
    faMeshBoundaryHaloTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef faMeshBoundaryHalo_H
#define faMeshBoundaryHalo_H

#include "mapDistributeBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class faMesh;

/*---------------------------------------------------------------------------*\
                      Class faMeshBoundaryHalo Declaration
\*---------------------------------------------------------------------------*/

class faMeshBoundaryHalo
:
    public mapDistributeBase
{
    // Private Data

        //- List of local input mesh faces required
        labelList inputMeshFaces_;

        //- Internal mapping from boundary index to compact
        labelList boundaryToCompact_;


public:

    // Declare name of the class and its debug switch
    ClassName("faMeshBoundaryHalo");


    // Constructors

        //- Default construct
        explicit faMeshBoundaryHalo(const label comm = UPstream::worldComm);

        //- Construct from mesh
        explicit faMeshBoundaryHalo(const faMesh& mesh);


    // Member Functions

        //- Clear out all parameters
        void clear();

        //- Redefine map and connectivity for a mesh
        void reset(const faMesh& mesh);

        //- The local data size (output)
        label haloSize() const;

        //- List of local input mesh faces required.
        //  \note will not correspond exactly to the boundary inside faces.
        //  Duplicates have been removed and it also contains the
        //  processor-local neighbour faces, which would otherwise not be
        //  handled by the distribute method.
        const labelList& inputMeshFaces() const noexcept
        {
            return inputMeshFaces_;
        }


    // Other

        //- Distribute sparse data.
        //  On output it is adjusted.
        template<class Type>
        void distributeSparse
        (
            List<Type>& fld,
            const labelUList& sparseInputLocations,
            const labelUList& compactOutputMapping
        ) const;

        //- Distribute sparse data.
        //  On output it is adjusted.
        template<class Type>
        void distributeSparse
        (
            List<Type>& fld,
            const labelUList& sparseInputLocations
        ) const;

        //- Distribute sparse data.
        //  The input field one enty per sparse id (inputMeshFaces).
        //  On output it will have for the input sparse
        //  The input field contains  location.
        template<class Type>
        void distributeSparse(List<Type>& fld) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "faMeshBoundaryHaloTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
