/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenCFD Ltd.
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
    Foam::domainDecompositionDryRun

Description
    Testing of domain decomposition for finite-volume meshes

SourceFiles
    domainDecompositionDryRun.C

\*---------------------------------------------------------------------------*/

#ifndef domainDecompositionDryRun_H
#define domainDecompositionDryRun_H

#include "fvMesh.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class domainDecompositionDryRun Declaration
\*---------------------------------------------------------------------------*/

class domainDecompositionDryRun
{
    // Private Data

        //- The mesh
        fvMesh mesh_;

        //- Optional non-standard file for decomposeParDict
        const fileName decompDictFile_;

        //- Optional override for numberOfSubdomains in decomposeParDict
        const label nDomainsOverride_;

        //- Optional override for method in decomposeParDict
        const word methodNameOverride_;


    // Private Members

        //- Write decomposition as volScalarField for visualization
        void writeVolField
        (
            const word& timeName,
            const labelUList& cellToProc
        ) const;

        //- Write cell distribution as VTU file (binary)
        void writeVTK
        (
            const fileName& file,
            const labelUList& cellToProc
        ) const;


public:

    // Constructors

        //- Construct from components.
        //  \param io the IOobject for mesh
        //  \param decompDictFile optional non-standard location for the
        //      decomposeParDict file
        //  \param nDomains optional override value for numberOfSubdomains
        //  \param methodName optional override value for decomposition method
        domainDecompositionDryRun
        (
            const IOobject& io,
            const fileName& decompDictFile = "",
            const label nDomains = 0,
            const word& methodName = ""
        );


    //- Destructor
    ~domainDecompositionDryRun() = default;


    // Member Functions

        //- The mesh
        const fvMesh& mesh() const noexcept
        {
            return mesh_;
        }

        //- Perform dry-run
        void execute(const bool writeCellDist, const bool verbose = false);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
