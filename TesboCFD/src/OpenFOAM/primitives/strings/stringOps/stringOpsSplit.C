/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "stringOps.H"
#include "Pair.H"
#include "Tuple2.H"
#include <vector>

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

inline Foam::word validateVariableName(const std::string& str)
{
    return Foam::stringOps::validate<Foam::word>
    (
        str,
        [](char c)
        {
            return (Foam::word::valid(c) || c == '/' || c == '{' || c == '}');
        }
    );
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

Foam::label Foam::stringOps::splitFunctionArgs
(
    const std::string& str,
    wordRes& args,
    List<Tuple2<word, string>>& namedArgs
)
{
    args.clear();
    namedArgs.clear();

    // Similar to code originally in functionObjectList (v2012 and earlier)
    // except that the function-name handling is now done prior to calling

    // (U)
    //     -> named = ()
    //     -> unnamed = (U)
    //
    // (patch=inlet, p)
    //     -> named = ((patch inlet))
    //     -> unnamed = (p)
    //
    // start=100, stop=200
    //     -> named = ((start 100) (stop 200) )
    //     -> unnamed = ()
    //
    // origin=(0 0 0) , scale=2 , normal=(0 0 1)


    // Use begin/end parse positions
    typedef Pair<std::string::size_type> rangeType;

    // For unnamed: beg/end range of each arg
    std::vector<rangeType> unnamed;

    // For named: list of beg/end ranges for (name, arg)
    std::vector<rangeType> named;

    // The beg/end range of the argument name
    rangeType argName(0, 0);

    // If argName is valid
    bool isNamed = false;

    // The depth of the argument parsing
    int argLevel = 0;

    const auto strLen = str.length();

    // Pass 1: parsing begin/end parse positions.

    for (std::string::size_type pos = 0, beg = 0; pos < strLen; ++pos)
    {
        const bool penultimate = ((pos + 1) == strLen);
        const char c = str[pos];

        if (c == ')')
        {
            --argLevel;
        }

        if (c == '=')
        {
            // Introducer for named argument
            argName = rangeType(beg, pos);
            isNamed = true;
            beg = pos + 1;
        }
        else if (c == '(')
        {
            ++argLevel;
        }
        else if (penultimate || (c == ','))  // OR: (c == ',' || c == ';')
        {
            if (penultimate && (c != ','))   // OR: (c != ',' && c != ';')
            {
                ++pos;  // Until the end, but do not include comma
            }

            if (argLevel == 0)
            {
                if (isNamed)
                {
                    named.push_back(argName);
                    named.push_back(rangeType(beg, pos));
                }
                else
                {
                    unnamed.push_back(rangeType(beg, pos));
                }
                isNamed = false;
                beg = pos + 1;
            }
        }
    }


    // Stage 2: Convert to concrete string and store


    // unnamed
    {
        const label nInputArgs = static_cast<label>(unnamed.size());
        args.resize(nInputArgs);

        label ngood = 0;
        for (label i = 0; i < nInputArgs; ++i)
        {
            const auto& arg = unnamed[i];

            args[ngood] = wordRe
            (
                word::validate
                (
                    str.substr(arg.first(), (arg.second()-arg.first()))
                )
            );

            // Only retain if non-empty
            if (!args[ngood].empty())
            {
                ++ngood;
            }
        }

        args.resize(ngood);
    }

    // named
    {
        const label nInputArgs = static_cast<label>(named.size());
        namedArgs.resize(nInputArgs/2);

        label ngood = 0;
        for (label i = 0; i < nInputArgs; i += 2)
        {
            const auto& name = named[i];
            const auto& arg = named[i+1];

            namedArgs[ngood].first() =
                validateVariableName
                (
                    str.substr(name.first(), (name.second()-name.first()))
                );

            namedArgs[ngood].second() =
                stringOps::trim
                (
                    str.substr(arg.first(), (arg.second()-arg.first()))
                );

            // Only retain if name is non-empty
            if (!namedArgs[ngood].first().empty())
            {
                ++ngood;
            }
        }

        namedArgs.resize(ngood);
    }

    // Return total number of arguments
    return (args.size() + namedArgs.size());
}


// ************************************************************************* //
