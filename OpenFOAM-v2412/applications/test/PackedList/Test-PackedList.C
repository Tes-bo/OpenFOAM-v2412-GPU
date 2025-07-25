/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

Description

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "uLabel.H"
#include "IOobject.H"
#include "IOstreams.H"
#include "IFstream.H"
#include "bitSet.H"
#include <climits>


using namespace Foam;

template<unsigned Width>
inline void reportInfo()
{
    const unsigned offset = PackedList<Width>::elem_per_block;

    unsigned useSHL = ((1u << (Width * offset)) - 1);
    unsigned useSHR = (~0u >> (sizeof(unsigned)*CHAR_BIT - Width * offset));

    Info<< nl
        << "PackedList<" << Width << ">" << nl
        << " max_value: " << PackedList<Width>::max_value << nl
        << " packing: " << PackedList<Width>::elem_per_block << nl
        << " utilization: " << (Width * offset) << nl;

    Info<< " Masking:" << nl
        << "  shift << "
        << unsigned(Width * offset) << nl
        << "  shift >> "
        << unsigned((sizeof(unsigned)*CHAR_BIT) - Width * offset)
        << nl;

    hex(Info);
    Info<< "   maskLower: "
        << PackedList<Width>::mask_lower(PackedList<Width>::elem_per_block)
        << nl
        << "      useSHL: " << useSHL << nl
        << "      useSHR: " << useSHR << nl;

    if (useSHL != useSHR)
    {
        Info<< "WARNING:  different results for SHL and SHR" << nl;
    }

    Info<< nl;
    dec(Info);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::addArgument("file .. fileN");

    argList::addBoolOption("mask", "report information about the bit masks");
    argList::addBoolOption("count", "test the count() method");
    argList::addBoolOption
    (
        "info",
        "print an ascii representation of the storage"
    );

    argList args(argc, argv, false, true);


    if (args.found("mask"))
    {
        Info<< "bit width: " << unsigned(sizeof(unsigned)*CHAR_BIT) << endl;
        reportInfo<1>();
        reportInfo<2>();
        reportInfo<3>();
        reportInfo<4>();
        reportInfo<5>();
        reportInfo<6>();
        reportInfo<7>();
        reportInfo<8>();
        reportInfo<9>();
        reportInfo<10>();
        reportInfo<11>();
        reportInfo<12>();
        reportInfo<13>();
        reportInfo<14>();
        reportInfo<15>();
        reportInfo<16>();

        return 0;
    }
    else if (args.size() <= 1)
    {
        args.printUsage();
    }


    for (label argi=1; argi < args.size(); ++argi)
    {
        const auto srcFile = args.get<fileName>(argi);
        Info<< nl << "reading " << srcFile << nl;

        IFstream ifs(srcFile);
        List<label> rawLst(ifs);

        bitSet packLst(rawLst);

        Info<< "size: " << packLst.size() << nl;

        if (args.found("count"))
        {
            unsigned int rawCount = 0;
            forAll(rawLst, elemI)
            {
                if (rawLst[elemI])
                {
                    rawCount++;
                }
            }
            Info<< "raw count: " << rawCount << nl
                << "packed count: " << packLst.count() << nl;
        }

        if (args.found("info"))
        {
            Info<< packLst.info();
        }

        Info<< nl;
        IOobject::writeDivider(Info);
    }

    return 0;
}

// ************************************************************************* //
