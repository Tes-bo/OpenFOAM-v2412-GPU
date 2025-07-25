/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023-2024 OpenCFD Ltd.
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
    Test loading of different gravity items

\*---------------------------------------------------------------------------*/

#include "MeshObject.H"
#include "gravityMeshObject.H"
#include "IOobjectList.H"
#include "IOstreams.H"
#include "argList.H"
#include "Time.H"

using namespace Foam;

void printInfo(const meshObjects::gravity& g)
{
    Pout<< "name:" << g.uniformDimensionedVectorField::name()
        << " type:" << g.type()
        << " value:" << g.value() << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::addBoolOption("checkout", "Test checkout with release");
    argList::addBoolOption("release", "Test release instead of delete");

    #include "setRootCase.H"
    #include "createTime.H"

    IOobjectList objects(runTime, runTime.constant());

    Info<< "Found: " << objects << nl << endl;

    for (const IOobject& io : objects.csorted<uniformDimensionedVectorField>())
    {
        if (io.name() == meshObjects::gravity::typeName)
        {
            const auto& g = meshObjects::gravity::New(runTime);
            printInfo(g);
        }
        else
        {
            const auto& g = meshObjects::gravity::New(io.name(), runTime);
            printInfo(g);
        }

        Pout<< "registered:"
            << flatOutput(runTime.sortedToc()) << nl << endl;
    }

    std::unique_ptr<meshObjects::gravity> release1;
    std::unique_ptr<meshObjects::gravity> release2;

    if (args.found("release"))
    {
        // Ugly!
        typedef
            MeshObject<Time, TopologicalMeshObject, meshObjects::gravity>
            parent_type;

        release1 = meshObjects::gravity::Release("g", runTime);
        release2 = meshObjects::gravity::Release("#none#", runTime);

        Info<< "release: " << Switch::name(bool(release1))
            << ", " << Switch::name(bool(release2)) << nl;

        Info<< "after Release: "
            << flatOutput(runTime.sortedToc()) << endl;

        // Do checkout by hand (ugly)
        if (args.found("checkout"))
        {
            if (release1)
            {
                release1->parent_type::checkOut();
            }

            if (release2)
            {
                release2->parent_type::checkOut();
            }

            Info<< "after checkout: "
                << flatOutput(runTime.sortedToc()) << endl;
        }
    }
    else if (args.found("checkout"))
    {
        // Do checkout as part of release
        release1 = meshObjects::gravity::Release("g", runTime, true);
        release2 = meshObjects::gravity::Release("#none#", runTime, true);

        Info<< "release: " << Switch::name(bool(release1))
            << ", " << Switch::name(bool(release2)) << nl;

        Info<< "after Release/Checkout(true) : "
            << flatOutput(runTime.sortedToc()) << endl;
    }
    else
    {
        meshObjects::gravity::Delete("g", runTime);
        meshObjects::gravity::Delete("#none#", runTime);

        Info<< "after Delete: "
            << flatOutput(runTime.sortedToc()) << endl;
    }


    if (meshObjects::gravity::Store(std::move(release1)))
    {
        Info<< "Store pointer" << endl;
    }

    if (release2)
    {
        release2.reset();
        Info<< "Clear pointer" << endl;
    }

    Info<< "Before exit: "
        << flatOutput(runTime.sortedToc()) << endl;

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
