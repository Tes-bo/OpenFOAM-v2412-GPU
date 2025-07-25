/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

Description
    Test accessors for PDRblock

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "PDRblock.H"

using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::noFunctionObjects();
    argList::addOption("dict", "file", "Alternative PDRblockMeshDict");

    #include "setRootCase.H"
    #include "createTime.H"

    const word dictName("PDRblockMeshDict");

    #include "setSystemRunTimeDictionaryIO.H"

    Info<< "Reading " << dictIO.name() << nl << endl;

    IOdictionary meshDict(dictIO);

    PDRblock mesh(meshDict, true);

    // Write summary
    Info<< "----------------" << nl
        << "Mesh Information" << nl
        << "----------------" << nl
        << "  " << "bounds: " << mesh.bounds() << nl
        << "  " << "nPoints: " << mesh.nPoints()
        << "  " << (mesh.sizes() + labelVector::one) << nl
        << "  " << "nCells: " << mesh.nCells()
        << "  " <<  mesh.sizes() << nl
        << "  " << "nFaces: " << mesh.nFaces() << nl
        << "  " << "nInternalFaces: " << mesh.nInternalFaces() << nl;

    Info<< "----------------" << nl
        << "Addressing" << nl
        << "----------------" << nl
        << "  " << "volume: " << mesh.bounds().volume() << nl
        << "  " << "avg-vol: "<< (mesh.bounds().volume() / mesh.nCells()) << nl;


    Info<< nl << "Check sizes: " << mesh.bounds().span() << nl;

    for (direction cmpt=0; cmpt < vector::nComponents; ++cmpt)
    {
        // Three different ways of getting the same information
        // not all are equally efficient

        scalar totalWidth = 0;

        // Using location
        forAll(mesh.grid()[cmpt], i)
        {
            totalWidth += mesh.grid()[cmpt].width(i);
        }

        Info<< vector::componentNames[cmpt] << "   min/max "
            << minMax(mesh.grid()[cmpt])
            << ' ' << mesh.grid()[cmpt].first()
            << ' ' << mesh.grid()[cmpt].last()
            << "  " << mesh.grid()[cmpt].size() << " cells"  << nl;


        // Use global (i,j,k) accessors, with mid-mesh for other directions

        const label nx = mesh.sizes().x() / (cmpt == vector::X ? 1 : 2);
        const label ny = mesh.sizes().y() / (cmpt == vector::Y ? 1 : 2);
        const label nz = mesh.sizes().z() / (cmpt == vector::Z ? 1 : 2);


        scalar totalSpan = 0;
        scalar totalDelta = 0;

        switch (cmpt)
        {
            case vector::X:
            {
                for (label i=0; i < nx; ++i)
                {
                    totalSpan += mesh.span(i, ny, nz).x();
                    totalDelta += mesh.dx(i);
                }
            }
            break;

            case vector::Y:
            {
                for (label j=0; j < ny; ++j)
                {
                    totalSpan += mesh.span(nx, j, nz).y();
                    totalDelta += mesh.dy(j);
                }
            }
            break;

            case vector::Z:
            {
                for (label k=0; k < nz; ++k)
                {
                    totalSpan += mesh.span(nx, ny, k).z();
                    totalDelta += mesh.dz(k);
                }
            }
            break;
        }

        Info<< "    width = " << totalWidth << " delta = " << totalDelta
            << " span = " << totalSpan << nl;
    }

    {
        const labelVector mid = mesh.sizes() / 2;

        Info<< nl << "Mid-mesh information at " << mid << nl
            << "  centre = " << mesh.C(mid) << nl
            << "  volume = " << mesh.V(mid) << nl
            << "  length = " << mesh.width(mid) << nl;
    }

    // Test findCell
    {
        Info<< nl << "findCell:" << nl;

        for
        (
            const point& pt
          : {
                mesh.bounds().centre(),
                mesh.bounds().min() - 0.1 * mesh.bounds().span(),
                mesh.bounds().max() + 0.1 * mesh.bounds().span()
            }
        )
        {
            labelVector ijk = mesh.findCell(pt);

            Info<< "    " << pt << " = " << ijk;
            if (cmptMin(ijk) < 0) Info<< " [not found]";
            Info<< nl;
        }
    }

    Info<< nl;

    // Fatal with FULLDEBUG
    {
        const bool oldThrowingError = FatalError.throwing(true);

        const label nx = mesh.sizes().x();
        const label ny = mesh.sizes().y();
        const label nz = mesh.sizes().z();

        // Here we get error in x-y-z order
        try
        {
            mesh.checkIndex(nx, ny, nz);
        }
        catch (const Foam::error& err)
        {
            Info<< nl << "Expected: Caught accesor error\n    "
                << err.message().c_str() << nl;
        }

        // No order for the error since it is called parameter-wise
        try
        {
            Info<< "centre at mesh(nx, ny, nz) = "
                << mesh.C(mesh.sizes()) << nl;
        }
        catch (const Foam::error& err)
        {
            Info<< nl << "Expected: Caught accesor error\n    "
                << err.message().c_str() << nl;
        }

        // Sizing error
        try
        {
            mesh.checkSizes(labelVector(nx+1, ny, nz));
        }
        catch (const Foam::error& err)
        {
            Info<< nl << "Expected: Caught sizing error\n    "
                << err.message().c_str() << nl;
        }

        try
        {
            mesh.checkSizes(labelMax/4);
        }
        catch (const Foam::error& err)
        {
            Info<< nl << "Expected: Caught sizing error\n    "
                << err.message().c_str() << nl;
        }

        FatalError.throwing(oldThrowingError);
    }


    Info<< "\nEnd\n" << nl;

    return 0;
}


// ************************************************************************* //
