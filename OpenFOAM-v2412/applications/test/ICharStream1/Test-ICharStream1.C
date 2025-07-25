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

Description

\*---------------------------------------------------------------------------*/

#include "SpanStream.H"
#include "wordList.H"
#include "IOstreams.H"
#include "argList.H"

#include <cctype>
#include <cstdio>

using namespace Foam;

Ostream& printString(Ostream& os, const char* first, const char* last)
{
    os << '"';
    for (; first != last; (void)++first)
    {
        os << *first;
    }
    os << '"';

    return os;
}


Ostream& printView(Ostream& os, const char* first, const char* last)
{
    char buf[4];
    os << label(last-first) << '(';

    for (; first != last; (void)++first)
    {
        const char c = *first;

        if (isprint(c))
        {
            os << c;
        }
        else if (c == '\t')
        {
            os << "\\t";
        }
        else if (c == '\n')
        {
            os << "\\n";
        }
        else
        {
            ::snprintf(buf, 4, "%02X", c);
            os << "\\x" << buf;
        }
    }
    os << ')';

    return os;
}


#if __cplusplus >= 201703L
Ostream& printView(Ostream& os, std::string_view s)
{
    return printView(os, s.begin(), s.end());
}
#endif


Ostream& printView(Ostream& os, stdFoam::span<char> s)
{
    return printView(os, s.begin(), s.end());
}


Ostream& printView(Ostream& os, const UList<char>& list)
{
    return printView(os, list.begin(), list.end());
}


Ostream& writeList(Ostream& os, const UList<char>& list)
{
    return printView(os, list.begin(), list.end());
}


Ostream& toString(Ostream& os, const UList<char>& list)
{
    return printString(os, list.begin(), list.end());
}


template<class BufType>
void printInfo(const BufType& buf)
{
    Info<< nl << "=========================" << endl;
    buf.print(Info);
    Info<< "addr: " << Foam::name(buf.list().cdata()) << nl;
    toString(Info, buf.list());
    Info<< nl << "=========================" << endl;
}


template<>
void printInfo(const List<char>& buf)
{
    Info<< nl << "=========================" << endl;
    toString(Info, buf);
    Info<< nl << "=========================" << endl;
}


void printTokens(Istream& is)
{
    label count = 0;
    token t;
    while (is.good())
    {
        is >> t;
        if (t.good())
        {
            ++count;
            Info<<"token: " << t << endl;
        }
    }

    Info<< count << " tokens" << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    // Buffer storage
    DynamicList<char> storage(16);

    OCharStream obuf(std::move(storage));
    obuf << 1002 << " " << "abcd" << " " << "def" << " " << 3.14159 << ";\n";

    // Move contents to output buffer
    printInfo(obuf);

    Info<<nl << "as string: ";
    toString(Info, obuf.list()) << endl;

    Info<< "transfer contents to a List" << endl;

    // Reclaim data storage from OCharStream -> ICharStream
    ICharStream ibuf(std::move(obuf));

    // OLD
    // ICharStream ibuf;
    // {
    //     List<char> data;
    //     obuf.swap(data);
    //     ibuf.swap(data);
    // }

    Info<< nl;
    Info<< nl << "input string:";
    printInfo(ibuf);

    Info<< nl << "orig output:";
    printInfo(obuf);

    printTokens(ibuf);

    Info<<nl << "after:";
    printInfo(ibuf);

    // This should also work
    ibuf.list() = 'X';

    Info<<nl << "overwritten with const value:";
    printInfo(ibuf);

    // Can also change content like this:
    {
        const int n = min(26, ibuf.size());

        for (int i=0; i<n; ++i)
        {
            ibuf.list()[i] = 'A' + i;
        }
    }

    Info<<nl << "directly written:";
    printInfo(ibuf);

    // Swap in/out an entirely new list storage:
    List<char> newvalues(52);
    {
        for (int i=0; i<26; ++i)
        {
            newvalues[2*i+0] = char('a' + i);
            newvalues[2*i+1] = char('A' + i);
        }
    }
    ibuf.swap(newvalues);

    Info<<nl << "after swap:";
    printInfo(ibuf);

    Info<<nl << "swapped out:";
    printInfo(newvalues);

    {
        icharstream is(std::move(newvalues));

        char c = 0;

        Info<< nl
            << "getting values from icharstream of "
            << is.list() << endl;

        // Info<< " (" << is.tellg() << " " << is.remaining() << ")";
        // Info<< "get:";
        while (is.get(c))
        {
            Info<< ' ' << c;
            // Info<< " (" << is.tellg() << " " << is.remaining() << ")";
        }
        Info<< " - end" << nl;

        // Info<< "remaining: " << is.list() << endl;
        // Info<< "remaining: " << is.remaining() << endl;

        // Manipulate the list view
        {
            UList<char> chars(is.list());
            Foam::reverse(chars);
        }

        is.rewind();

        Info<< "get:";
        while (is.get(c))
        {
            Info<< ' ' << c;
        }
        Info<< " - end" << nl;
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
