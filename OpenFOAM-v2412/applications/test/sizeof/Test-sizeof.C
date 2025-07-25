/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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
    Test the sizeof various classes.

\*---------------------------------------------------------------------------*/

#include "bool.H"
#include "Switch.H"
#include "dictionary.H"
#include "zero.H"
#include "IOstreams.H"
#include "PstreamBuffers.H"
#include "argList.H"
#include "Time.H"
#include "IOobject.H"
#include "scalarField.H"

namespace Foam
{
   class hasBoolClass
   {
   public:
      bool b_;

      hasBoolClass(const bool val=false)
      :
          b_(false)
      {}
   };

}


using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    cout<<"sizeof\n------\n";

    if (true)
    {
        cout<<"IOstreamOption:" << sizeof(IOstreamOption) << nl;
        cout<<"IOstream:" << sizeof(IOstream) << nl;
        cout<<"token:" << sizeof(token) << nl;
        cout<<"Istream:" << sizeof(Istream) << nl;
        cout<<"Ostream:" << sizeof(Ostream) << nl;

        cout<<"ISstream:" << sizeof(ISstream) << nl;
        cout<<"OSstream:" << sizeof(OSstream) << nl;

        cout<<"IPstream:" << sizeof(IPstream) << nl;
        cout<<"OPstream:" << sizeof(OPstream) << nl;
    }

    {
        zero x;
        cout<<"zero:" << sizeof(x) << nl;
    }
    #if 0
    {
        argList x(argc, argv);
        cout<<"argList:" << sizeof(x) << nl;

        TimePaths y(x);
        cout<<"TimePaths:" << sizeof(y) << nl;
    }
    #endif
    {
        zero x;
        cout<<"zero:" << sizeof(x) << nl;
    }
    {
        bool x(0);
        cout<<"bool:" << sizeof(x) << nl;
    }
    {
        hasBoolClass x(true);
        cout<<"hasBoolClass:" << sizeof(x) << nl;
    }

    {
        Switch x("n");
        cout<<"Switch:" << sizeof(x) << nl;
        cout<<"Switch::switchType=" << sizeof(Switch::switchType) << nl;
    }

    {
        scalar x(0);
        cout<<"scalar:" << sizeof(x) << nl;
    }

    {
        label x(0);
        cout<<"label:" << sizeof(x) << nl;
    }

    {
        cout<<"short:" << sizeof(short) << nl;
        cout<<"int:" << sizeof(int) << nl;
        cout<<"long:" << sizeof(long) << nl;
        cout<<"float:" << sizeof(float) << nl;
        cout<<"double:" << sizeof(double) << nl;
    }

    {
        cout<<"string:" << sizeof(Foam::string) << nl;
    }

    cout<<"IOobjectOption:" << sizeof(Foam::IOobjectOption) << nl;
    cout<<"IOobject:" << sizeof(Foam::IOobject) << nl;
    cout<<"IOstream:" << sizeof(Foam::IOstream) << nl;
    cout<<"PstreamBuffers:" << sizeof(Foam::PstreamBuffers) << nl;
    cout<<"Time:" << sizeof(Foam::Time) << nl;

    cout<<"tmp<>:" << sizeof(tmp<scalarField>) << nl;

    Info << "---\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
