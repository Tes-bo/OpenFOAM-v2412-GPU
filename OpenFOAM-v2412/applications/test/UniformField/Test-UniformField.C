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
    Test-UniformField

Description
    Test uniform list/field
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "vector.H"
#include "IOstreams.H"
#include "UniformField.H"

using namespace Foam;

template<class T>
void printInfo(const UniformList<T>& list, const label i=0)
{
    Info<< nl
        << "value: " << list.value() << nl
        << "cast:  " << static_cast<const T&>(list) << nl
        << "[" << i << "] = " << list[i] << nl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    {
        UniformField<scalar> fld(3.14159);

        printInfo(fld, -100);

        // Change value
        fld.value() *= 0.5;

        Info<< nl << "/= 2 " << nl;

        printInfo(fld, -100);
    }

    {
        UniformField<vector> fld(vector(1, 2, -1));

        printInfo(fld);

        // Change value
        fld.value() *= 0.5;

        Info<< nl << "/= 2 " << nl;

        printInfo(fld);
    }


    Info<< "\nEnd\n" << endl;

    return 0;
}

// ************************************************************************* //
