/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2023 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "lagrangianReconstructor.H"
#include "passivePositionParticleCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::lagrangianReconstructor::verbose_ = 1;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lagrangianReconstructor::lagrangianReconstructor
(
    const fvMesh& mesh,
    const PtrList<fvMesh>& procMeshes,
    const PtrList<labelIOList>& faceProcAddressing,
    const PtrList<labelIOList>& cellProcAddressing
)
:
    mesh_(mesh),
    procMeshes_(procMeshes),
    faceProcAddressing_(faceProcAddressing),
    cellProcAddressing_(cellProcAddressing)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::lagrangianReconstructor::reconstructPositions
(
    const word& cloudName
) const
{
    // Start with empty cloud
    passivePositionParticleCloud lagrangianPositions
    (
        mesh_,
        Foam::zero{},
        cloudName
    );

    forAll(procMeshes_, meshi)
    {
        const labelList& cellMap = cellProcAddressing_[meshi];
        const labelList& faceMap = faceProcAddressing_[meshi];

        // Use a special particle that does not try to find the particle on
        // the mesh. This is to be able to handle particles originating
        // from a different processor. This can happen with some
        // functionObjects - e.g. extractEulerianParticles.
        // These particles should be
        // - written in the old format
        passivePositionParticleCloud lpi(procMeshes_[meshi], cloudName, false);

        for (const passivePositionParticle& ppi : lpi)
        {
            const label mappedCell =
            (
                (ppi.cell() >= 0)
              ? cellMap[ppi.cell()]
              : -1
            );

            // Inverting sign if necessary and subtracting 1 from
            // faceProcAddressing
            const label mappedTetFace =
            (
                (ppi.tetFace() >= 0)
              ? mag(faceMap[ppi.tetFace()]) - 1
              : -1
            );

            if ((ppi.cell() >= 0) && (ppi.tetFace() >= 0))
            {
                // cell,face succesfully mapped. Coordinates inside the cell
                // should be same
                lagrangianPositions.append
                (
                    new passivePositionParticle
                    (
                        mesh_,
                        ppi.coordinates(),
                        mappedCell,
                        mappedTetFace,
                        ppi.procTetPt(mesh_, mappedCell, mappedTetFace)
                    )
                );
            }
            else
            {
                // No valid coordinates. Two choices:
                // - assume reconstructed mesh contains the position so do
                //   a locate with the (reconstructed) mesh
                // - preserve -1 as cell id, maintain the read location
                lagrangianPositions.append
                (

                    //- Option 1: locate on reconstructed mesh
                    //new passivePositionParticle
                    //(
                    //    mesh_,
                    //    ppi.location(),
                    //    mappedCell
                    //)

                    //- Option 2: maintain read location
                    new passivePositionParticle
                    (
                        mesh_,
                        Zero,   // position
                        -1,     // celli
                        -1,     // tetFacei
                        -1,     // tetPti
                        ppi.location()
                    )
                );
            }
        }
    }


    IOPosition<passivePositionParticleCloud>(lagrangianPositions).write();

    // Force writing of "positions" too, if specified via the InfoSwitch
    if (particle::writeLagrangianPositions)
    {
        IOPosition<passivePositionParticleCloud>
        (
            lagrangianPositions,
            cloud::geometryType::POSITIONS
        ).write();
    }

    return lagrangianPositions.size();
}


void Foam::lagrangianReconstructor::reconstructAllFields
(
    const word& cloudName,
    const IOobjectList& cloudObjs,
    const wordRes& selectedFields
)
{
    do
    {
        #undef  doLocalCode
        #define doLocalCode(Type)                                             \
        {                                                                     \
            this->reconstructFields<Type>                                     \
            (                                                                 \
                cloudName,                                                    \
                cloudObjs,                                                    \
                selectedFields                                                \
            );                                                                \
                                                                              \
            this->reconstructFieldFields<Type>                                \
            (                                                                 \
                cloudName,                                                    \
                cloudObjs,                                                    \
                selectedFields                                                \
            );                                                                \
        }

        doLocalCode(label);
        doLocalCode(scalar);
        doLocalCode(vector);
        doLocalCode(sphericalTensor);
        doLocalCode(symmTensor);
        doLocalCode(tensor);

        #undef doLocalCode
    }
    while (false);
}


// ************************************************************************* //
