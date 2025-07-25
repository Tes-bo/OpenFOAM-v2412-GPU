/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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
    Foam::fvMeshSubsetter

Description
    Extends Foam::fvMeshSubset with two-step subsetting
    (uses polyTopoChange modification).

SourceFiles
    fvMeshSubsetter.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_fvMeshSubsetter_H
#define Foam_fvMeshSubsetter_H

#include "fvMeshSubset.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class fvMeshSubsetter Declaration
\*---------------------------------------------------------------------------*/

class fvMeshSubsetter
:
    public fvMeshSubset
{
    // Private Member Functions

        //- Forwarding to Foam::removeCells
        void removeCellsImpl
        (
            const bitSet& cellsToRemove,
            const labelList& exposedFaces,
            const labelList& patchIDs,
            const bool syncPar
        );

        //- No copy construct
        fvMeshSubsetter(const fvMeshSubsetter&) = delete;

        //- No copy assignment
        void operator=(const fvMeshSubset&) = delete;

public:

    // Constructors

        //- Inherit constructors from fvMeshSubset
        using fvMeshSubset::fvMeshSubset;


    // Member Functions

        //- Inherit all one-step subsetting
        using fvMeshSubset::setCellSubset;


    // Two-step subsetting

        //- Get labels of exposed faces.
        //  These are
        //  - internal faces that become boundary faces
        //  - coupled faces that become uncoupled (since one of the
        //    sides gets deleted)
        labelList getExposedFaces
        (
            const bitSet& selectedCells,
            const bool syncPar = true
        ) const;

        //- Get labels of exposed faces.
        //  These are
        //  - internal faces that become boundary faces
        //  - coupled faces that become uncoupled (since one of the
        //    sides gets deleted)
        labelList getExposedFaces
        (
            const label regioni,
            const labelUList& regions,
            const bool syncPar = true
        ) const;

        //- For every exposed face (from above getExposedFaces)
        //  Uses supplied (existing!) patches
        void setCellSubset
        (
            const bitSet& selectedCells,
            const labelList& exposedFaces,
            const labelList& patchIDs,
            const bool syncPar = true
        );

        //- For every exposed face (from above getExposedFaces)
        //  Uses supplied (existing!) patches
        void setCellSubset
        (
            const label regioni,
            const labelList& regions,
            const labelList& exposedFaces,
            const labelList& patchIDs,
            const bool syncPar = true
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
