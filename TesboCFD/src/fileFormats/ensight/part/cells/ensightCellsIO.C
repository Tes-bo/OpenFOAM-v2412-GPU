/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2024 OpenCFD Ltd.
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

#include "ensightCells.H"
#include "ensightOutput.H"
#include "InfoProxy.H"
#include "boundBox.H"
#include "polyMesh.H"
#include "cellModel.H"
#include "globalIndex.H"
#include "globalMeshData.H"
#include "manifoldCellsMeshObject.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::ensightCells::writePolysConnectivity
(
    ensightGeoFile& os,
    const polyMesh& mesh,
    const ensightCells& part,
    const labelList& pointToGlobal,
    const bool parallel
)
{
    constexpr ensightCells::elemType etype(ensightCells::elemType::NFACED);

    const label nTotal = part.total(etype);
    const labelUList& addr = part.cellIds(etype);

    if (!nTotal)
    {
        return;
    }

    const IntRange<int> senders =
    (
        parallel
      ? Pstream::subProcs()
      : IntRange<int>()
    );


    if (Pstream::master())
    {
        os.writeKeyword(ensightCells::key(etype));
        os.write(nTotal);
        os.newline();
    }

    // Number of faces per polyhedral (1/line in ASCII)
    {
        labelList send
        (
            ensightOutput::Detail::getPolysNFaces(mesh, addr)
        );

        if (Pstream::master())
        {
            // Main
            os.writeLabels(send);

            // Others
            for (const int proci : senders)
            {
                IPstream fromOther(Pstream::commsTypes::scheduled, proci);
                labelList recv(fromOther);

                os.writeLabels(recv);
            }
        }
        else if (senders)
        {
            OPstream toMaster
            (
                Pstream::commsTypes::scheduled,
                Pstream::masterNo()
            );

            toMaster << send;
        }
    }


    // Number of points for each polyhedral face (1/line in ASCII)
    {
        labelList send
        (
            ensightOutput::Detail::getPolysNPointsPerFace(mesh, addr)
        );

        if (Pstream::master())
        {
            // Main
            os.writeLabels(send);

            // Others
            for (const int proci : senders)
            {
                IPstream fromOther(Pstream::commsTypes::scheduled, proci);
                labelList recv(fromOther);

                os.writeLabels(recv);
            }
        }
        else if (senders)
        {
            OPstream toMaster
            (
                Pstream::commsTypes::scheduled,
                Pstream::masterNo()
            );

            toMaster << send;
        }
    }

    ///const cellList& meshCells = mesh.cells();
    const cellList& meshCells = manifoldCellsMeshObject::New(mesh).cells();

    // List of points id for each face of the above list
    if (Pstream::master())
    {
        // Main
        ensightOutput::writePolysPoints
        (
            os,
            mesh,
            addr,
            pointToGlobal
        );

        // Others
        for (const int proci : senders)
        {
            IPstream fromOther(Pstream::commsTypes::scheduled, proci);
            cellList  cells(fromOther);
            labelList addr(fromOther);
            faceList  faces(fromOther);
            labelList owner(fromOther);

            ensightOutput::writePolysPoints
            (
                os,
                cells,
                addr,
                faces,
                owner
            );
        }
    }
    else if (senders)
    {
        // Renumber faces to use global point numbers
        faceList faces(mesh.faces());
        ListListOps::inplaceRenumber(pointToGlobal, faces);

        OPstream toMaster
        (
            Pstream::commsTypes::scheduled,
            Pstream::masterNo()
        );

        toMaster
            << meshCells
            << addr
            << faces
            << mesh.faceOwner();
    }
}


void Foam::ensightCells::writeShapeConnectivity
(
    ensightGeoFile& os,
    const polyMesh& mesh,
    const ensightCells::elemType etype,
    const ensightCells& part,
    const labelList& pointToGlobal,
    const bool parallel
)
{
    if (etype == ensightCells::elemType::NFACED)
    {
        FatalErrorInFunction
            << "Called for ensight NFACED cell. Programming error\n"
            << exit(FatalError);
    }

    const label nTotal = part.total(etype);
    const labelUList& addr = part.cellIds(etype);

    if (!nTotal)
    {
        return;
    }


    const IntRange<int> senders =
    (
        parallel
      ? Pstream::subProcs()
      : IntRange<int>()
    );


    if (Pstream::master())
    {
        os.writeKeyword(ensightCells::key(etype));
        os.write(nTotal);
        os.newline();
    }


    // Primitive shape - get subset and renumber
    cellShapeList shapes(mesh.cellShapes(), addr);

    ListListOps::inplaceRenumber(pointToGlobal, shapes);

    if (Pstream::master())
    {
        ensightOutput::writeCellShapes(os, shapes);

        for (const int proci : senders)
        {
            IPstream fromOther(Pstream::commsTypes::scheduled, proci);
            cellShapeList recv(fromOther);

            ensightOutput::writeCellShapes(os, recv);
        }
    }
    else if (senders)
    {
        OPstream toMaster
        (
            Pstream::commsTypes::scheduled,
            Pstream::masterNo()
        );

        toMaster << shapes;
    }
}


void Foam::ensightCells::write
(
    ensightGeoFile& os,
    const polyMesh& mesh,
    bool parallel
) const
{
    const ensightCells& part = *this;

    parallel = parallel && Pstream::parRun();

    // Renumber the points/faces into unique points

    label nPoints = 0;  // Total number of points
    labelList pointToGlobal;  // local point to unique global index
    labelList uniqueMeshPointLabels;  // unique global points

    nPoints = meshPointMapppings
    (
        mesh,
        pointToGlobal,
        uniqueMeshPointLabels,
        parallel
    );

    ensightOutput::Detail::writeCoordinates
    (
        os,
        part.index(),
        part.name(),
        nPoints,  // nPoints (global)
        UIndirectList<point>(mesh.points(), uniqueMeshPointLabels),
        parallel  //!< Collective write?
    );


    for (label typei=0; typei < ensightCells::nTypes; ++typei)
    {
        const auto etype = ensightCells::elemType(typei);

        if (etype == ensightCells::elemType::NFACED)
        {
            writePolysConnectivity
            (
                os,
                mesh,
                part,
                pointToGlobal,
                parallel
            );
        }
        else
        {
            writeShapeConnectivity
            (
                os,
                mesh,
                etype,
                part,
                pointToGlobal,
                parallel
            );
        }
    }
}


void Foam::ensightCells::writeBox
(
    ensightGeoFile& os,
    const boundBox& bb,
    const label partIndex,
    const word& partName
)
{
    pointField points;
    cellShapeList shapes;

    if (UPstream::master())
    {
        points = bb.hexCorners();
        shapes.emplace_back(cellModel::HEX, identity(8));
    }

    ensightOutput::Detail::writeCoordinates
    (
        os,
        partIndex,
        partName,
        8,          // nPoints (global)
        points,
        false       // serial only! (parallel=false)
    );

    if (UPstream::master())
    {
        os.writeKeyword(ensightCells::key(ensightCells::elemType::HEXA8));
        os.write(shapes.size());  // one cell (global)
        os.newline();

        ensightOutput::writeCellShapes(os, shapes);
    }
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const InfoProxy<ensightCells>& iproxy
)
{
    const auto& part = *iproxy;

    os << part.name().c_str();

    for (label typei=0; typei < ensightCells::nTypes; ++typei)
    {
        const auto etype = ensightCells::elemType(typei);

        os  << ' ' << ensightCells::elemNames[etype]
            << ':' << part.total(etype);
    }
    os  << nl;

    return os;
}


// ************************************************************************* //
