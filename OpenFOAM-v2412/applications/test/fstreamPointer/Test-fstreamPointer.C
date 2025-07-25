/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2022 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Application
    Test-fstreamPointer

Description
    Low-level fstream tests

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "refPtr.H"
#include "fstreamPointer.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noParallel();

    argList::addOption
    (
        "output",
        "file",
        "Output file name for cat"
    );

    argList::addArgument("file1");
    argList::addArgument("...");
    argList::addArgument("fileN");

    argList::noMandatoryArgs();

    #include "setRootCase.H"

    if (args.size() <= 1)
    {
        InfoErr<< "\nNo input files specified .. stopping\n" << endl;
        return 0;
    }

    refPtr<std::ostream> osRef;

    fileName outputName;
    if (args.readIfPresent("output", outputName))
    {
        InfoErr<< "output: " << outputName;

        IOstreamOption::compressionType comp(IOstreamOption::UNCOMPRESSED);
        if (outputName.has_ext("gz"))
        {
            comp = IOstreamOption::COMPRESSED;
            outputName.remove_ext();

            InfoErr<< " [compress]";
        }
        InfoErr<< nl;

        osRef.reset(ofstreamPointer(outputName, comp).release());
    }
    else
    {
        osRef.ref(std::cout);
        InfoErr<< "output: stdout" << nl;
    }
    auto& os = osRef.ref();

    for (label argi = 1; argi < args.size(); ++argi)
    {
        const auto inputName = args.get<fileName>(argi);

        InfoErr<< "input: " << inputName;

        ifstreamPointer isPtr(inputName);

        if (!isPtr.get() || !isPtr->good())
        {
            InfoErr<< " (not good)" << nl;
            continue;
        }
        InfoErr<< nl;

        auto& is = *isPtr;

        // Loop getting single characters
        // - not efficient, but that is not the test here anyhow

        char c;
        while (is.get(c))
        {
            os << c;
        }
    }

    InfoErr<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
