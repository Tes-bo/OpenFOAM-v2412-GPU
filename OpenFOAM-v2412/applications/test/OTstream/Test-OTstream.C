/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "ITstream.H"
#include "OTstream.H"
#include "primitiveFields.H"
#include "argList.H"

using namespace Foam;

void printTokens(const UList<token>& toks)
{
    label count = 0;
    for (const token& t : toks)
    {
        Info<< "token: " << t.info() << nl;
        ++count;
    }

    Info<< count << " tokens" << nl << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    // Test fields

    {
        scalarField fld1({1, 2, 3, 4});

        OTstream os;

        fld1.writeEntry("field", os);

        Info<< "Field: " << fld1 << nl
            << "Tokens: " << flatOutput(os) << nl;

        printTokens(os);
    }

    {
        scalarField fld1(10, scalar(5));

        OTstream os;

        fld1.writeEntry("field", os);

        Info<< "Field: " << fld1 << nl
            << "Tokens: " << flatOutput(os) << nl;

        printTokens(os);
    }

    {
        vector val(1,2, 3);

        OTstream os;

        os << val;

        Info<< "Value: " << val << nl
            << "Tokens: " << flatOutput(os) << nl;

        printTokens(os);
    }

    {
        bool val(true);

        OTstream os;

        os << val;

        Info<< "Value: " << val << nl
            << "Tokens: " << flatOutput(os) << nl;

        printTokens(os);
    }

    {
        tensorField fld1(1, tensor::I);

        OTstream os;

        fld1.writeEntry("field", os);

        Info<< "Field: " << fld1 << nl
            << "Tokens: " << flatOutput(os) << nl;

        printTokens(os);
    }

    {
        labelList fld1(identity(5));

        OTstream os;

        fld1.writeEntry("field", os);

        Info<< "Field: " << fld1 << nl
            << "Tokens: " << flatOutput(os) << nl;

        printTokens(os);
    }


    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
