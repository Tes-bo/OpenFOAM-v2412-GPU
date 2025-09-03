/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

Application
    noise

Group
    grpPostProcessingUtilities

Description
    Perform noise analysis of pressure data.

    The utility provides a light wrapper around the run-time selectable
    noise model.  Current options include:
    - point, and
    - surface noise.

Usage
    \verbatim
    noiseModel      surfaceNoise; // pointNoise

    surfaceNoiseCoeffs
    {
        windowModel     Hanning;

        HanningCoeffs
        {
            // Window overlap percentage
            overlapPercent  50;
            symmetric       yes;
            extended        yes;

            // Optional number of windows, default = all available
            nWindow         5;
        }


        // Input files list
        files       ("postProcessing/faceSource1/surface/patch/patch.case";)

        // Surface reader
        reader      ensight;

        // Surface writer
        writer      ensight;

        // Collate times for ensight output
        // - ensures geometry is only written once
        writeOptions
        {
            ensight
            {
                collateTimes true;
            }
        }

        // Number of samples in sampling window, default = 2^16 (=65536)
        N               4096;

        // Write interval for FFT data, default = 1
        fftWriteInterval 100;
    }
    \endverbatim


See also
    noiseModel.H
    windowModel.H

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "noiseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Perform noise analysis of pressure data"
    );

    argList::noCheckProcessorDirectories();

    argList::addOption("dict", "file", "Alternative noiseDict");

    #include "setRootCase.H"

    // As much as possible avoid synchronised operation
    fileHandler().distributed(true);

    #include "createTime.H"

    const word dictName("noiseDict");
    #include "setSystemRunTimeDictionaryIO.H"

    Info<< "Reading " << dictIO.name() << nl << endl;

    IOdictionary dict(dictIO);

    autoPtr<noiseModel> model(noiseModel::New(dict, runTime));
    model->calculate();

    Info<< nl << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
