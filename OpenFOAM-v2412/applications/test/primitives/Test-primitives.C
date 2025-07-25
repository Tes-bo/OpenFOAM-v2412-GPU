/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2023 OpenCFD Ltd.
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
    Test-primitives

Description
    Parsing etc for primitives.

\*---------------------------------------------------------------------------*/

#include "scalar.H"
#include "FlatOutput.H"
#include "SpanStream.H"
#include "NASCore.H"
#include "parsing.H"
#include "Tuple2.H"
#include "Switch.H"
#include "dictionary.H"
#include "primitiveFields.H"
#include "labelVector.H"

using namespace Foam;

// Shadow fileFormats::NASCore::readNasScalar
inline scalar readNasScalar(const std::string& str)
{
    return fileFormats::NASCore::readNasScalar(str);
}


// As a function
inline Switch readSwitch(const std::string& str)
{
    Switch sw(Switch::find(str));

    // Trap bad input and raise as exit error, not abort
    if (sw.bad())
    {
        FatalErrorInFunction
            << "Unknown switch " << str << nl
            << exit(FatalError);
    }

    return sw;
}


void printInfo(const Switch& sw)
{
    Info<< "Switch " << sw.c_str() << " (enum=" << label(sw.type()) << ")\n";
}


Ostream& toString(Ostream& os, const UList<char>& list)
{
    os << '"';
    for (const char c : list)
    {
        os << c;
    }
    os << '"';

    return os;
}


template<class T1, class T2>
void printValPair(const char* desc, const T1& val1, const T2& val2)
{
    Info<< desc << ' ' << val1 << ' ' << val2 << nl;
}


// The 'naive' implementation. Not exact at end points
inline scalar naive_lerp(const scalar a, const scalar b, const scalar t)
{
    return a + t*(b - a);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class T>
bool hadParsingError
(
    const std::pair<bool, std::string>& input,
    const std::pair<bool, T>& result,
    std::string errMsg
)
{
    if (result.first)
    {
        if (input.first)
        {
            Info<< "(pass) parsed "
                << input.second << " = " << result.second << nl;
        }
        else
        {
            Info<< "(fail) unexpected success for " << input.second << nl;
        }
    }
    else
    {
        if (input.first)
        {
            Info<< "(fail) unexpected";
        }
        else
        {
            Info<< "(pass) expected";
        }

        Info<< " failure " << input.second << "  >> " << errMsg.c_str() << nl;
    }

    return (input.first != result.first);
}


template<class T>
unsigned testParsing
(
    T (*function)(const std::string&),
    std::initializer_list<std::pair<bool, std::string>> tests
)
{
    unsigned nFail = 0;
    string errMsg;

    // Expect some failures
    const bool oldThrowingError = FatalError.throwing(true);
    const bool oldThrowingIOErr = FatalIOError.throwing(true);

    for (const std::pair<bool, std::string>& test : tests)
    {
        std::pair<bool, T> result(false, T());

        try
        {
            result.second = function (test.second);
            result.first = true;
        }
        catch (const Foam::error& err)
        {
            errMsg = err.message();
        }

        if (test.first != result.first)
        {
            ++nFail;
        }

        hadParsingError(test, result, errMsg);
    }

    FatalError.throwing(oldThrowingError);
    FatalIOError.throwing(oldThrowingIOErr);

    return nFail;
}


int main(int argc, char *argv[])
{
    unsigned nFail = 0;

    {
        const float a = 1e8f, b = 1.0f;

        Info<< "lerp exact: "
            << (a == lerp(a, b, 0.0f)) << " "
            << (b == lerp(a, b, 1.0f)) << nl;

        Info<< "naive lerp exact: "
            << (a == naive_lerp(a, b, 0.0f)) << " "
            << (b == naive_lerp(a, b, 1.0f)) << nl;
    }

    {
        const scalar a = 1e24, b = 1.0;

        Info<< "lerp exact: "
            << (a == lerp(a, b, 0.0)) << " "
            << (b == lerp(a, b, 1.0)) << nl;

        Info<< "naive lerp exact: "
            << (a == naive_lerp(a, b, 0.0)) << " "
            << (b == naive_lerp(a, b, 1.0)) << nl;
    }

    {
        const vector a(vector::uniform(1e24)), b(vector::uniform(1.0));

        Info<<"lerp exact: "
            << (a == lerp(a, b, 0.0f)) << " "
            << (b == lerp(a, b, 1.0f)) << nl;

        Info<< "lerp: "
            << lerp(vector::uniform(0), vector::uniform(100), 0.5) << nl;
    }

    {
        const lerpOp1<vector> half(0.5);

        const vector a(vector::uniform(20));
        const vector b(vector::uniform(100));

        Info<< "lerp half: "
            << a << " : " << b << " => " << half(a, b) << nl;
    }

    {
        const labelVector a(labelVector::uniform(10000));
        const labelVector b(labelVector::uniform(1));

        Info<< "lerp (labelVector) = "
            << lerp(a, b, 0.1) << nl;
    }

    {
        const scalar a(0);
        const scalar b(100);

        Info<< "lerp of " << a << " : " << b << nl;

        for (const double t : { 0.0, 0.5, 1.0, -0.5, 1.5 })
        {
            Info<< "    " << t << " = " << lerp(a, b, t) << nl;
        }
    }


    // No yet
    #if 0
    {
        const label a(10000);
        const label b(1);

        Info<<"lerp (label) = "
            << label(lerp(a, b, 0.1)) << nl;
    }

    {
        const bool a(true);
        const bool b(false);

        Info<<"lerp (bool) = "
            << (lerp(a, b, 0.5)) << nl;
    }
    #endif

    {
        const sphericalTensor a(10), b(20);

        Info<<"lerp exact: "
            << (a == lerp(a, b, 0.0f)) << " "
            << (b == lerp(a, b, 1.0f)) << nl;
//        Info<< "lerp: "
//            << lerp(vector::uniform(0), vector::uniform(100), 0.5) << nl;
    }

    {
        const tensor a(tensor::uniform(1e24));
        const tensor b(tensor::uniform(0));

        Info<<"lerp exact: "
            << (a == lerp(a, b, 0.0f)) << " "
            << (b == lerp(a, b, 1.0f)) << nl;
//        Info<< "lerp: "
//            << lerp(vector::uniform(0), vector::uniform(100), 0.5) << nl;
    }


    {
        Info<< nl << "Test Switch parsing:" << nl;
        nFail += testParsing
        (
            &readSwitch,
            {
                { false, "True" },
                { true,  "false" },
                { true,  "on" },
                { false, "None" },
                { true,  "yes" },
                { true,  "none" },
                { false, "default" },
            }
        );

        Info<< nl << "Test Switch defaults:" << nl;

        dictionary dict;
        dict.add("key1" , "true");
        dict.add("key2" , "off");
        dict.add("key3" , "any");

        for (const word& k : { "key", "key1", "key2" , "key3" })
        {
            Switch sw1(k, dict, Switch::YES);
            Switch sw2(k, dict, Switch::NO);

            Info<< nl;
            printInfo(sw1);
            printInfo(sw2);
            Info<< "bool " << bool(sw1) << nl;

            sw1.negate();
            sw2.negate();
            Info<< "negated" << nl;
            printInfo(sw1);
            printInfo(sw2);
        }
    }

    {
        Info<< nl << "Test readDouble: (small=" << doubleScalarVSMALL
            << " great=" << doubleScalarVSMALL << "):" << nl;
        nFail += testParsing
        (
            &readDouble,
            {
                { false, "" },
                { false, "  " },
                { false, " xxx " },
                { false, " 1234E-" },
                { false, " 1234E junk" },
                { true,  " 3.14159 " },
                { true,  " 31.4159E-1 "  },
                { false, " 100E1000 "  },
                { true,  " 1E-40 "  },
                { true,  " 1E-305 "  },
                { true,  " 1E-37 "  },
                { true,  " 1E-300 "  },
            }
        );
    }

    {
        Info<< nl << "Test readFloat: (small=" << floatScalarVSMALL
            << " great=" << floatScalarVGREAT << "):" << nl;

        nFail += testParsing
        (
            &readFloat,
            {
                { true,  " 3.14159 " },
                { true,  " 31.4159E-1 "  },
                { false, " 31.4159E200 "  },
                { true,  " 31.4159E20 "  },
                { true,  " 1E-40 "  },
                { true,  " 1E-305 "  },
                { true,  " 1E-37 "  },
                { true,  " 1E-300 "  },
            }
        );
    }

    {
        Info<< nl << "Test readNasScalar:" << nl;
        nFail += testParsing
        (
            &readNasScalar,
            {
                { true,  " 3.14159 " },
                { true,  " 31.4159E-1 "  },
                { true,  " 314.159-2 "  },
                { true,  " 31.4159E200 "  },
                { true,  " 31.4159E20 "  },
                { true,  " 1E-40 "  },
                { true,  " 1E-305 "  },
                { true,  " 1E-37 "  },
                { true,  " 1E-300 "  },
            }
        );
    }

    {
        Info<< nl << "Test readInt32 (max=" << INT32_MAX << "):" << nl;
        nFail += testParsing
        (
            &readInt32,
            {
                { false, " 3.14159 " },
                { false, " 31E1 " },
                { false, " 31.4159E-1 "  },
                { true,  "100"  },
                { true,  "\t2147483644"  },
                { false, "   2147483700  "  },
            }
        );
    }

    {
        Info<< nl << "Test readUint32 (max="
            << unsigned(UINT32_MAX) << "):" << nl;
        nFail += testParsing
        (
            &readUint32,
            {
                { true,  "\t2147483644"  },
                { true,  " 2147483700  "  },
                { true,  " 4294967295  "  },
                { false, " 4294968000  "  },
            }
        );
    }

    if (true)
    {
        #ifdef WM_DP
        typedef float otherType;
        #else
        typedef double otherType;
        #endif

        Info<< nl << "Test raw binary read of scalar list:"
            << " write " << sizeof(otherType)
            << " read " << sizeof(scalar) << nl;

        List<otherType> srcList(15);
        List<scalar> dstList;  // Read back

        forAll(srcList, i)
        {
            srcList[i] = 1 + 10*i;
        }

        OCharStream os(IOstreamOption::BINARY);
        os << srcList;

        // Recover written contents
        ICharStream is(os.release(), IOstreamOption::BINARY);
        is.setScalarByteSize(sizeof(otherType));

        Info<< "Stream scalar-size ("
            << label(is.scalarByteSize()) << ") is native: "
            << Switch(is.checkScalarSize<otherType>()) << nl;


        token firstToken(is);
        Info<< "List has " << firstToken.info() << " scalar items" << nl;

        dstList.resize(firstToken.labelToken(), 3.14159);

        is.beginRawRead();
        // for (scalar& val : dstList)
        // {
        //     val = readRawScalar(is);
        // }
        readRawScalar(is, dstList.data(), dstList.size());
        is.endRawRead();

        Info<< "Wrote " << flatOutput(srcList) << nl
            << "Read " << flatOutput(dstList) << nl;
    }

    #ifdef COMPAT_OPENFOAM_ORG
    Info<< nl << "compatibility sizes" << nl
        << "name  com  org" << nl
        << "----  ---  ---" << nl;

    printValPair("SMALL", SMALL, small);
    printValPair("GREAT", GREAT, great);
    printValPair("VSMALL", VSMALL, vSmall);
    printValPair("VGREAT", VGREAT, vGreat);
    printValPair("ROOTSMALL", ROOTSMALL, rootSmall);
    printValPair("ROOTGREAT", ROOTGREAT, rootGreat);
    #else
    Info<< nl << "no compatibility sizes" << nl;
    #endif

    if (nFail)
    {
        Info<< nl << "failed " << nFail << " tests" << nl;
        return 1;
    }

    Info<< nl << "passed all tests" << nl;
    return 0;
}

// ************************************************************************* //
