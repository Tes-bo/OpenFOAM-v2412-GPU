/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2020 OpenCFD Ltd.
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

Application
    Test-externalFileCoupler

Description
    Test of master/slave communication etc.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "externalFileCoupler.H"

using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noParallel();
    argList::addOption("sleep", "N", "sleep to add between calls");
    argList::addOption("max", "N", "max number of calls (default: 1000)");
    argList::addBoolOption("slave", "run as slave");

    #include "setRootCase.H"

    const label maxCount = args.getOrDefault<label>("max", 1000);
    const label sleeping = args.getOrDefault<label>("sleep", 0);

    externalFileCoupler coupler;

    if (args.found("slave"))
    {
        const word role = "slave";
        const word other = "master";
        Info<< "Running as " << role << " max=" << maxCount
            << " (sleep " << sleeping << ')' << endl;

        for (label count = 0; count < maxCount; ++count)
        {
            // Wait for master, but stop if status=done was seen
            Info<< role << '[' << count << "] wait for " << other << endl;

            if (!coupler.waitForMaster())
            {
                Info<< role << ": stopping. status=done was detected" << endl;
                break;
            }

            if (sleeping)
            {
                sleep(sleeping);
            }

            // Info<< role << ": switch to " << other << endl;
            coupler.useMaster();
        }

        Info<< role << ": exiting" << endl;
    }
    else
    {
        const word role = "master";
        const word other = "slave";
        Info<< "Running as " << role << " max=" << maxCount
            << " (sleep " << sleeping << ')' << endl;

        for (label count = 0; count < maxCount; ++count)
        {
            // Wait for slave, but stop if status=done was seen

            Info<< role << '[' << count << "] wait for " << other << endl;

            if (!coupler.waitForSlave())
            {
                Info<< role << ": stopping. status=done was detected" << endl;
                break;
            }

            if (sleeping)
            {
                sleep(sleeping);
            }

            // Info<< role << ": switch to " << other << endl;
            coupler.useSlave();
        }

        // shudown - slave should notice and terminate
        Info<< role << ": exiting" << endl;
    }

    return 0;
}


// ************************************************************************* //
