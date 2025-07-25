/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Application
    Test-fft

Description
    Very simple fft tests

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fft.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

// Simple forward transform
tmp<complexField> forward1D(const tmp<complexField>& input)
{
    const label len = input().size();

    return fft::forwardTransform(input, List<int>({len})) / len;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();

    #include "setRootCase.H"

    // Simple ones - http://www.sccon.ca/sccon/fft/fft3.htm
    {
        complexField input(8, Zero);
        input[0] = 1;

        tmp<complexField> toutput = forward1D(input);

        Info<< nl
            << "input = " << input << nl
            << "output = " << toutput << nl;
    }

    {
        complexField input(8, Zero);
        input[1] = 1;

        tmp<complexField> toutput = forward1D(input);

        Info<< nl
            << "input = " << input << nl
            << "output = " << toutput << nl;
    }

    Info<< nl << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
