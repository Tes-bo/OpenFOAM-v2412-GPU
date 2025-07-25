/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 M. Janssens
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

#include "processorColour.H"
#include "processorLduInterface.H"
#include "processorTopologyNew.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(processorColour, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::processorColour::colour
(
    const lduMesh& mesh,
    labelList& procColour
)
{
    procColour.resize_nocopy(Pstream::nProcs(mesh.comm()));
    procColour = -1;

    // Re-use processor-topology analysis

    labelListList procNeighbours(Pstream::nProcs(mesh.comm()));

    // Fill my entry
    {
        const lduInterfacePtrsList patches = mesh.interfaces();

        auto& procToProcs = procNeighbours[Pstream::myProcNo(mesh.comm())];
        label n = 0;
        forAll(patches, patchi)
        {
            if (patches.set(patchi))
            {
                if (isA<processorLduInterface>(patches[patchi]))
                {
                    n++;
                }
            }
        }

        procToProcs.resize_nocopy(n);
        n = 0;
        forAll(patches, patchi)
        {
            if (patches.set(patchi))
            {
                const auto* ppPtr = isA<processorLduInterface>(patches[patchi]);
                if (ppPtr)
                {
                    procToProcs[n++] = ppPtr->neighbProcNo();
                }
            }
        }
    }
    // Send to master
    Pstream::gatherList(procNeighbours, UPstream::msgType(), mesh.comm());


    // Use greedy algorithm for now
    // (see e.g. https://iq.opengenus.org/graph-colouring-greedy-algorithm/)

    if (Pstream::master(mesh.comm()))
    {
        DynamicList<label> front;
        front.append(0);    // start from processor 0

        DynamicList<label> newFront;
        while (front.size())
        {
            //Pout<< "Front:" << front << endl;

            newFront.clear();
            for (const label proci : front)
            {
                if (procColour[proci] == -1)
                {
                    const labelList& nbrs = procNeighbours[proci];
                    const UIndirectList<label> nbrColour(procColour, nbrs);

                    for
                    (
                        label colouri = 0;
                        colouri < Pstream::nProcs();
                        colouri++
                    )
                    {
                        if (!nbrColour.found(colouri))
                        {
                            procColour[proci] = colouri;
                            for (label nbrProci : nbrs)
                            {
                                if (procColour[nbrProci] == -1)
                                {
                                    newFront.append(nbrProci);
                                }
                            }
                            break;
                        }
                    }
                }
            }

            front = std::move(newFront);
        }
    }

    Pstream::broadcast(procColour, mesh.comm());


    //if (false)
    //{
    //    volScalarField volColour
    //    (
    //        IOobject
    //        (
    //            "colour",
    //            mesh.time().timeName(),
    //            mesh,
    //            IOobject::NO_READ,
    //            IOobject::NO_WRITE,
    //            IOobject::NO_REGISTER
    //        ),
    //        mesh,
    //        dimensionedScalar(dimless, procColour[Pstream::myProcNo()]),
    //        fvPatchFieldBase::zeroGradientType()
    //    );
    //    volColour.write();
    //}

    const label nColours = max(procColour)+1;

    if (debug)
    {
        Info<< typeName << " : coloured " << Pstream::nProcs(mesh.comm())
            << " processors with in total " << nColours << " colours" << endl;
    }

    return nColours;
}


void Foam::processorColour::walkFront
(
    const lduMesh& mesh,
    DynamicList<label>& front,
    labelList& cellColour
)
{
    // Colour with the (min) coupled (global) patch

    const lduAddressing& addr = mesh.lduAddr();
    const label* const __restrict__ uPtr = addr.upperAddr().begin();
    const label* const __restrict__ lPtr = addr.lowerAddr().begin();
    const label* const __restrict__ ownStartPtr = addr.ownerStartAddr().begin();
    const label* const __restrict__ losortStartAddrPtr =
        addr.losortStartAddr().begin();
    const label* const __restrict__ losortAddrPtr = addr.losortAddr().begin();

    DynamicList<label> newFront;
    while (true)
    {
        newFront.clear();
        for (const label celli : front)
        {
            const label colouri = cellColour[celli];

            {
                const label fStart = ownStartPtr[celli];
                const label fEnd = ownStartPtr[celli + 1];

                for (label facei=fStart; facei<fEnd; facei++)
                {
                    const label nbr =
                    (
                        lPtr[facei] == celli
                      ? uPtr[facei]
                      : lPtr[facei]
                    );
                    if (cellColour[nbr] == -1)
                    {
                        cellColour[nbr] = colouri;
                        newFront.append(nbr);
                    }
                }
            }
            {
                const label fStart = losortStartAddrPtr[celli];
                const label fEnd = losortStartAddrPtr[celli + 1];

                for (label i=fStart; i<fEnd; i++)
                {
                    label facei = losortAddrPtr[i];
                    const label nbr =
                    (
                        lPtr[facei] == celli
                      ? uPtr[facei]
                      : lPtr[facei]
                    );
                    if (cellColour[nbr] == -1)
                    {
                        cellColour[nbr] = colouri;
                        newFront.append(nbr);
                    }
                }
            }

        }

        if (newFront.empty())
        {
            break;
        }

        front.transfer(newFront);
    }
}


Foam::label Foam::processorColour::cellColour
(
    const lduMesh& mesh,
    labelList& cellColour,
    labelList& patchToColour
)
{
    // Colour with the (min) coupled (global) patch. Compact to have colour 0
    // if a single region. Problematic if same patch on multiple disconnected
    // regions.

    const lduAddressing& addr = mesh.lduAddr();
    const lduInterfacePtrsList patches = mesh.interfaces();

    const label nCells = addr.size();

    patchToColour.resize_nocopy(patches.size());
    patchToColour = -1;
    cellColour.resize_nocopy(nCells);
    cellColour = -1;


    label colouri = 0;

    // Starting front from patch faceCells
    DynamicList<label> front;
    forAll(patches, inti)
    {
        if
        (
            patches.set(inti)
        && !isA<processorLduInterface>(patches[inti])
        )
        {
            // 'global' interface. Seed faceCells with patch index

            if (patchToColour[inti] == -1)
            {
                patchToColour[inti] = colouri++;
            }

            const auto& fc = patches[inti].faceCells();
            forAll(fc, face)
            {
                const label cell = fc[face];
                if (cellColour[cell] != patchToColour[inti])
                {
                    cellColour[cell] = patchToColour[inti];
                    front.append(cell);
                }
            }
        }
    }


    // Walk current mesh
    walkFront(mesh, front, cellColour);


    // Any unvisited cell is still labelMax. Put into colour 0.
    for (auto& colour : cellColour)
    {
        if (colour == -1)
        {
            colour = 0;
        }
    }

    if (debug)
    {
        Info<< typeName << " : coloured " << nCells
            << " cells with in total " << colouri << " colours" << endl;
    }

    return colouri;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::processorColour::processorColour(const lduMesh& mesh)
:
    MeshObject_type(mesh)
{
    nColours_ = colour(mesh, *this);
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

const Foam::processorColour& Foam::processorColour::New(const lduMesh& mesh)
{
    auto* ptr =
        mesh.thisDb().objectRegistry::template getObjectPtr<processorColour>
        (
            processorColour::typeName
        );

    if (!ptr)
    {
        ptr = new processorColour(mesh);

        //regIOobject::store(static_cast<MoveableMeshObject<lduMesh>*>(ptr));
        regIOobject::store(ptr);
    }

    return *ptr;
}


// ************************************************************************* //
