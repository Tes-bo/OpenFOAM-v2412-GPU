/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2023 OpenCFD Ltd.
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

#include "List.H"
#include "Istream.H"
#include "token.H"
#include "contiguous.H"
#include <memory>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T, int SizeMin>
Foam::DynamicList<T, SizeMin>::DynamicList(Istream& is)
:
    List<T>(),
    capacity_(0)
{
    this->readList(is);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T, int SizeMin>
bool Foam::DynamicList<T, SizeMin>::readBracketList(Istream& is)
{
    DynamicList<T, SizeMin>& list = *this;

    is.fatalCheck(FUNCTION_NAME);

    token tok(is);

    is.fatalCheck
    (
        "DynamicList<T>::readBracketList(Istream&) : reading first token"
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
        constexpr label chunkSize = 128;
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
        list.resize(list.capacity());

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

                chunks[nChunks] = chunkType(new List<T>(chunkSize));
                ++nChunks;
                localIndex = 0;
            }

            is  >> chunks[nChunks-1]->operator[](localIndex);
            ++localIndex;
            ++totalCount;

            is.fatalCheck
            (
                "DynamicList<T>::readBracketList(Istream&) : "
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
        list.setCapacity_nocopy(totalCount);
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


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class T, int SizeMin>
Foam::Istream& Foam::DynamicList<T, SizeMin>::readList(Istream& is)
{
    DynamicList<T, SizeMin>& list = *this;

    is.fatalCheck(FUNCTION_NAME);

    token tok(is);

    is.fatalCheck("DynamicList<T>::readList(Istream&) : reading first token");

    if (tok.isCompound())
    {
        // Compound: simply transfer contents

        list.clearStorage();  // Remove old contents
        list.transfer
        (
            tok.transferCompoundToken<List<T>>(is)
        );
    }
    else if (tok.isLabel())
    {
        // Label: could be int(..), int{...} or just a plain '0'

        const label len = tok.labelToken();

        // Resize to length required
        list.resize_nocopy(len);

        if (is.format() == IOstreamOption::BINARY && is_contiguous<T>::value)
        {
            // Binary and contiguous

            if (len)
            {
                Detail::readContiguous<T>
                (
                    is,
                    list.data_bytes(),
                    list.size_bytes()
                );

                is.fatalCheck
                (
                    "DynamicList<T>::readList(Istream&) : "
                    "reading binary block"
                );
            }
        }
        else if (std::is_same<char, typename std::remove_cv<T>::type>::value)
        {
            // Special treatment for char data (binary I/O only)
            const auto oldFmt = is.format(IOstreamOption::BINARY);

            if (len)
            {
                // read(...) includes surrounding start/end delimiters
                is.read(list.data_bytes(), list.size_bytes());

                is.fatalCheck
                (
                    "DynamicList<char>::readList(Istream&) : "
                    "reading binary block"
                );
            }

            is.format(oldFmt);
        }
        else
        {
            // Begin of contents marker
            const char delimiter = is.readBeginList("List");

            if (len)
            {
                if (delimiter == token::BEGIN_LIST)
                {
                    auto iter = list.begin();
                    const auto last = list.end();

                    // Contents
                    for (/*nil*/; (iter != last); (void)++iter)
                    {
                        is >> *iter;

                        is.fatalCheck
                        (
                            "DynamicList<T>::readList(Istream&) : "
                            "reading entry"
                        );
                    }
                }
                else
                {
                    // Uniform content (delimiter == token::BEGIN_BLOCK)

                    T elem;
                    is >> elem;

                    is.fatalCheck
                    (
                        "DynamicList<T>::readList(Istream&) : "
                        "reading the single entry"
                    );

                    // Fill with the value
                    UList<T>::operator=(elem);
                }
            }

            // End of contents marker
            is.readEndList("List");
        }
    }
    else if (tok.isPunctuation(token::BEGIN_LIST))
    {
        // "(...)" : read read as bracketed list

        is.putBack(tok);
        this->readBracketList(is);

        // Could also simply be done with emplace_back for each element
        // but prefer the same mechanism as List::readList to avoid
        // intermediate resizing

        // // list.clear();  // Clear addressing, leave storage intact
        // //
        // // is >> tok;
        // // is.fatalCheck(FUNCTION_NAME);
        // //
        // // while (!tok.isPunctuation(token::END_LIST))
        // // {
        // //     is.putBack(tok);
        // //     is >> list.emplace_back();
        // //
        // //     is.fatalCheck
        // //     (
        // //         "DynamicList<T>::readList(Istream&) : "
        // //         "reading entry"
        // //     );
        // //
        // //     is >> tok;
        // //     is.fatalCheck(FUNCTION_NAME);
        // // }
    }
    else
    {
        list.clear();  // Clear old contents

        FatalIOErrorInFunction(is)
            << "incorrect first token, expected <int> or '(', found "
            << tok.info() << nl
            << exit(FatalIOError);
    }

    return is;
}


// ************************************************************************* //
