/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2023 OpenCFD Ltd.
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
    Test-IFstream

Description
    Tests on low-level reading

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Fstream.H"
#include "OSspecific.H"
#include "etcFiles.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noParallel();
    argList::noParallel();
    argList::addOption("ignore", "file", "Test readRaw with ignore");

    #include "setRootCase.H"

    // Test with etc/controlDict (mandatory, from distribution)

    if (!args.found("ignore"))
    {
        const fileName inputFile
        (
            findEtcFile("controlDict", true, 0007)
        );

        Info<< nl << "Test getLine" << nl << inputFile << nl;
        IFstream is(inputFile);

        string buf;

        // Test that buffer resizing works as expected

        Info<< "buf: " << buf.size() << '/' << buf.capacity()
            << nl << nl;

        bool skipDoc = true;

        while (is.good())
        {
            if (skipDoc || is.lineNumber() % 2)
            {
                is.getLine(buf);

                Info<< is.lineNumber() << ": "
                    << buf.size() << '/' << buf.capacity() << ' '
                    << buf << nl;

                if (buf.starts_with("Documentation"))
                {
                    skipDoc = false;

                    Info<< "Found Documentation: skip until closing '}'"
                        << " line number will now be incorrect"
                        << nl;

                    auto n = is.getLine(nullptr, '}');
                    Info<< is.lineNumber() << ": [" << label(n) << ']' << nl;
                }
            }
            else
            {
                auto n = is.getLine(nullptr);  // default: '\n'
                Info<< is.lineNumber() << ": [" << label(n) << ']' << nl;
            }
        }
    }

    fileName testFile;
    if (args.readIfPresent("ignore", testFile))
    {
        if (testFile.has_ext("gz"))
        {
            testFile.remove_ext();
            Info<< "stripping extraneous .gz ending" << endl;
        }

        IFstream is(testFile);
        auto& stdStream = is.stdStream();

        List<char> buffer(1000);

        Info<< "Test readRaw with: " << is.name()
            << " compressed:" << int(is.compression())
            << " file-size:" << is.fileSize() << nl;

        for (int iter = 0; is.good() && iter < 1000; ++iter)
        {
            Info<< "iter:" << iter;
            if (iter % 2)
            {
                Info<< " [read]  ";
                is.readRaw(buffer.data(), buffer.size());
            }
            else
            {
                Info<< " [ignore]";
                is.readRaw(nullptr, buffer.size() / 2);
            }

            Info<< " : " << stdStream.gcount() << endl;
        }
    }


    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
