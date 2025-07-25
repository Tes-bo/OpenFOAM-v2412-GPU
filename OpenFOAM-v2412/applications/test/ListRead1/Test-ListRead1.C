/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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
    Test-ListRead1

Description
    List reading

\*---------------------------------------------------------------------------*/

#include "OSspecific.H"
#include "argList.H"
#include "wordRes.H"

#include "IOstreams.H"
#include "Fstream.H"
#include "scalar.H"
#include "vector.H"

#include "labelRange.H"
#include "scalarList.H"
#include "HashOps.H"
#include "ListOps.H"
#include "IndirectList.H"
#include "SubList.H"
#include "SliceList.H"
#include "ListPolicy.H"

#include <list>
#include <numeric>
#include <functional>

using namespace Foam;

label chunkSize = 128;

template<class T>
bool readBracketList(List<T>& list, Istream& is)
{
    is.fatalCheck(FUNCTION_NAME);

    token tok(is);

    is.fatalCheck
    (
        "List<T>::readBracketList(Istream&) : reading first token"
    );

    if (!tok.isPunctuation(token::BEGIN_LIST))
    {
        is.putBack(tok);
        return false;
    }

    {
        // "(...)" : read element-wise.
        // Uses chunk-wise reading to avoid too many re-allocations
        // and avoids relocation of contiguous memory until all of the reading
        // is completed. Chunks are wrapped as unique_ptr to ensure proper
        // cleanup on failure.

        // The choice of chunk-size is somewhat arbitrary...
        // constexpr label chunkSize = 128;
        typedef std::unique_ptr<List<T>> chunkType;

        is >> tok;
        is.fatalCheck(FUNCTION_NAME);

        if (tok.isPunctuation(token::END_LIST))
        {
            // Trivial case, an empty list
            list.clear();
            return true;
        }

        // Use all storage
        //private:// list.resize(list.capacity());

        // Start with a few slots, recover current memory where possible
        List<chunkType> chunks(16);
        if (list.empty())
        {
            chunks[0] = chunkType(new List<T>(chunkSize));
        }
        else
        {
            chunks[0] = chunkType(new List<T>(std::move(list)));
        }

        label nChunks = 1;      // Active number of chunks
        label totalCount = 0;   // Total number of elements
        label localIndex = 0;   // Chunk-local index

        InfoErr
            << nl << "initial chunk: " << chunks[0]->size() << endl;

        while (!tok.isPunctuation(token::END_LIST))
        {
            is.putBack(tok);

            if (chunks[nChunks-1]->size() <= localIndex)
            {
                // Increase number of slots (doubling)
                if (nChunks >= chunks.size())
                {
                    chunks.resize(2*chunks.size());
                }

                InfoErr<< "new chunk" << endl;
                chunks[nChunks] = chunkType(new List<T>(chunkSize));
                ++nChunks;
                localIndex = 0;
            }

            is  >> chunks[nChunks-1]->operator[](localIndex);
            ++localIndex;
            ++totalCount;

            InfoErr
                << "    chunk=" << nChunks
                << " index=" << localIndex
                << " total=" << totalCount << nl;

            is.fatalCheck
            (
                "List<T>::readBracketList(Istream&) : "
                "reading entry"
            );

            is >> tok;
            is.fatalCheck(FUNCTION_NAME);
        }

        // Simple case
        if (nChunks == 1)
        {
            list = std::move(*(chunks[0]));
            list.resize(totalCount);
            return true;
        }

        // Destination
        //private:// list.setCapacity_nocopy(totalCount);
        list.resize_nocopy(totalCount);
        auto dest = list.begin();

        for (label chunki = 0; chunki < nChunks; ++chunki)
        {
            List<T> currChunk(std::move(*(chunks[chunki])));
            chunks[chunki].reset(nullptr);

            const label localLen = min(currChunk.size(), totalCount);

            dest = std::move
            (
                currChunk.begin(),
                currChunk.begin(localLen),
                dest
            );

            totalCount -= localLen;
        }
    }

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noParallel();
    argList::noFunctionObjects();
    argList::addOption("chunk-size",  "value", "change read chunk size");
    argList::addArgument("file1 .. fileN");

    argList args(argc, argv, false, true);

    args.readIfPresent("chunk-size", chunkSize);

    Info<< "chunk-size: " << chunkSize << nl;

    if (args.size() <= 1)
    {
        InfoErr<< "Provide a file or files to test" << nl;
    }
    else
    {
        for (label argi=1; argi < args.size(); ++argi)
        {
            const auto input = args.get<fileName>(argi);
            IFstream is(input);

            while (!is.eof())
            {
                labelList list;

                readBracketList(list, is);
                Info<< "read: " << flatOutput(list) << endl;
            }
        }
    }

    return 0;
}

// ************************************************************************* //
