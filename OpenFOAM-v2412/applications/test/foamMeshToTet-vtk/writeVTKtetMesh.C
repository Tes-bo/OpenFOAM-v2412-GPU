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

\*---------------------------------------------------------------------------*/

#include "polyMesh.H"
#include "Fstream.H"
#include "tetMatcher.H"
#include "foamVtkInternalMeshWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

void writeVTKtetMesh(const fileName& output, const polyMesh& mesh_)
{
    #if WM_LABEL_SIZE == 64
    # define FOAM_VTK_LABEL_NAME "vtktypeint64"
    #else
    # define FOAM_VTK_LABEL_NAME "int"
    #endif

    const faceList& faces = mesh_.faces();
    const labelList& faceOwner = mesh_.faceOwner();

    label nTets = 0;

    for (label facei = 0; facei < mesh_.nInternalFaces(); ++facei)
    {
        nTets += 2 * faces[facei].nTriangles();
    }
    for (const face& f : mesh_.boundaryMesh().faces())
    {
        nTets += f.nTriangles();
    }
    const label nPoints = (mesh_.nPoints() + mesh_.nCells());

    OFstream os(output + ".vtk");

    Info<< "Write: " << os.name() << nl;

    os  << "# vtk DataFile Version 5.1" << nl
        << "tet-mesh" << nl
        << "ASCII" << nl
        << "DATASET UNSTRUCTURED_GRID" << nl
        << nl;

    os  << "POINTS " << nPoints << " float" << nl;
    for (const point& p : mesh_.points())
    {
        os  << float(p.x()) << ' '
            << float(p.y()) << ' '
            << float(p.z()) << nl;
    }
    for (const point& p : mesh_.cellCentres())
    {
        os  << float(p.x()) << ' '
            << float(p.y()) << ' '
            << float(p.z()) << nl;
    }
    os  << nl;

    os  << "CELLS " << (1 + nTets) << ' ' << (4 * nTets) << nl;

    os  << "OFFSETS " << FOAM_VTK_LABEL_NAME << nl
        << 0;  // begin offset = 0
    {
        label offset = 0;
        for (label teti = 0; teti < nTets; ++teti)
        {
            offset += 4;
            os << ' ' << offset;
        }
        os << nl << nl;
    }

    labelList nLocalTets(mesh_.nCells(), Zero);

    os  << nl
        << "CONNECTIVITY " << FOAM_VTK_LABEL_NAME << nl;

    for (label celli = 0; celli < mesh_.nCells(); ++celli)
    {
        const cell& cFaces = mesh_.cells()[celli];

        if (tetMatcher::test(mesh_, celli))
        {
            // Tet: no cell-centre decomposition

            const label facei = cFaces[0];
            const face& f0 = faces[facei];

            // Get the other point from f1. Tbd: check if not duplicate face
            // (ACMI / ignoreBoundaryFaces_).
            const face& f1 = faces[cFaces[1]];
            label apexi = -1;
            forAll(f1, fp)
            {
                apexi = f1[fp];
                if (!f0.found(apexi))
                {
                    break;
                }
            }

            const label p0 = f0[0];
            label p1 = f0[1];
            label p2 = f0[2];

            if (faceOwner[facei] == celli)
            {
                std::swap(p1, p2);
            }

            ++nLocalTets[celli];
            os << p0 << ' ' << p1 << ' ' << p2 << ' ' << apexi << nl;
        }
        else
        {
            for (const label facei : cFaces)
            {
                const face& f = faces[facei];

                label fp0 = mesh_.tetBasePtIs()[facei];

                // Fallback
                if (fp0 < 0)
                {
                    fp0 = 0;
                }

                const label p0 = f[fp0];
                label fp = f.fcIndex(fp0);
                for (label i = 2; i < f.size(); ++i)
                {
                    label p1 = f[fp];
                    fp = f.fcIndex(fp);
                    label p2 = f[fp];

                    if (faceOwner[facei] == celli)
                    {
                        std::swap(p1, p2);
                    }

                    ++nLocalTets[celli];
                    os  << p0 << ' ' << p1 << ' ' << p2 << ' '
                        << (mesh_.nPoints() + celli) << nl;
                }
            }
        }
    }

    os  << nl
        << "CELL_TYPES " << nTets << nl;

    for (label teti = 0; teti < nTets; ++teti)
    {
        if (teti) os  << ' ';
        os  << vtk::cellType::VTK_TETRA;
    }
    os << nl;

    os << nl << "CELL_DATA " << nTets << nl
       << "FIELD FieldData " << 1 << nl;

    os << "cellID " << 1 << ' ' << nTets << " int" << nl;
    forAll(nLocalTets, celli)
    {
        label n = nLocalTets[celli];
        while (n--)
        {
            os << ' ' << celli;
        }
        os << nl;
    }
}

} // End namespace Foam


// ************************************************************************* //
