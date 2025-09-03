/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2023 OpenCFD Ltd.
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
    Hand-written parsing of STL ASCII format

\*---------------------------------------------------------------------------*/

#include "STLAsciiParse.H"
#include "STLReader.H"
#include "stringOps.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

static inline std::string perrorEOF(std::string expected)
{
    return "Premature EOF while reading '" + expected + "'";
}


static inline std::string perrorParse(std::string expected, std::string found)
{
    return "Parse error. Expecting '" + expected + "' found '" + found + "'";
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Foam
{
namespace Detail
{

//- A lexer for parsing STL ASCII files.
//  Returns DynamicList(s) of points and facets (zoneIds).
//  The facets are within a solid/endsolid grouping
class STLAsciiParseManual
:
    public Detail::STLAsciiParse
{
    enum scanState
    {
        scanSolid = 0,
        scanFacet,
        scanLoop,
        scanVerts,
        scanEndLoop,
        scanEndFacet,
        scanEndSolid
    };

    scanState state_;

    std::string errMsg_;

    //- Like std:csub_match
    typedef std::pair<const char*, const char*> tokenType;

    // Tokenized line
    DynamicList<tokenType> tokens_;

    //- Tokenize
    inline std::string::size_type tokenize(const char *p, const char *pe)
    {
        const char* start = p;
        tokens_.clear();

        // Find not space
        while (p < pe && isspace(*p))
        {
            if (*p == '\n' && lineNum_)
            {
                ++lineNum_;
            }
            ++p;
        }

        while (p != pe)
        {
            const char* beg = p;

            // Find space
            while (p < pe && !isspace(*p))
            {
                ++p;
            }
            tokens_.emplace_back(beg, p);

            // Find next
            while (p < pe && isspace(*p))
            {
                if (*p == '\n')
                {
                    ++lineNum_;
                    return (p - start);
                }
                ++p;
            }
        }

        return (p - start);
    }


public:

    //- Construct with the estimated number of triangles in the STL
    STLAsciiParseManual(const label nTrisEstimated)
    :
        Detail::STLAsciiParse(nTrisEstimated)
    {}

    //- Execute parser
    void execute(std::istream& is);
};

} // End namespace Detail
} // End namespace Foam


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Length of the input read buffer
#define INBUFLEN 16384

void Foam::Detail::STLAsciiParseManual::execute(std::istream& is)
{
    if (!is)
    {
        return;
    }

    // Buffering
    char inbuf[INBUFLEN];
    std::streamsize pending = 0;

    lineNum_ = 0;

    state_ = scanSolid;
    errMsg_.clear();

    // Line-oriented processing loop
    while (is)
    {
        if (pending >= INBUFLEN)
        {
            // We overfilled the buffer while trying to scan a token...
            FatalErrorInFunction
                << "buffer full while scanning near line " << lineNum_ << nl;
            break;
        }

        char *data = inbuf + pending;   // current data buffer
        const std::streamsize buflen = INBUFLEN - pending; // space in buffer

        is.read(data, buflen);
        const std::streamsize gcount = is.gcount();

        if (gcount <= 0)
        {
            // EOF
            // If scanning for next "solid" this is a valid way to exit, but
            // an error if scanning for the initial "solid" or any other token

            switch (state_)
            {
                case scanSolid:
                {
                    if (!lineNum_) errMsg_ = perrorEOF("solid");
                    break;
                }
                case scanFacet:    { errMsg_ = perrorEOF("facet"); break; }
                case scanLoop:     { errMsg_ = perrorEOF("outer loop"); break; }
                case scanVerts:    { errMsg_ = perrorEOF("vertex"); break; }
                case scanEndLoop:  { errMsg_ = perrorEOF("endloop"); break; }
                case scanEndFacet: { errMsg_ = perrorEOF("endfacet"); break; }
                case scanEndSolid: { errMsg_ = perrorEOF("endsolid"); break; }
            }

            // Terminate the parsing loop
            break;
        }

        // p,pe = Ragel parsing point and parsing end (default naming)
        // eof  = Ragel EOF point (default naming)

        char *p = inbuf;
        char *pe = data + gcount;

        // Line-oriented: search backwards to find last newline
        {
            --pe;
            while (*pe != '\n' && pe >= inbuf)
            {
                --pe;
            }
            ++pe;
        }

        std::string cmd;
        do
        {
            // Tokenize
            const auto parsedLen = tokenize(p, pe);
            p += parsedLen;
            if (!parsedLen || tokens_.empty())
            {
                break;
            }

            // Ensure consistent case on the first token
            cmd.assign(tokens_[0].first, tokens_[0].second);
            stringOps::lower(cmd);

            // Handle all expected parse states
            switch (state_)
            {
                case scanSolid:
                {
                    if (cmd == "solid")
                    {
                        if (tokens_.empty())
                        {
                            beginSolid(word::null);
                        }
                        else
                        {
                            beginSolid
                            (
                                word::validate
                                (
                                    tokens_[1].first,
                                    tokens_[1].second
                                )
                            );
                        }

                        state_ = scanFacet;  // Next state
                    }
                    else
                    {
                        errMsg_ = perrorParse("solid", cmd);
                    }
                    break;
                }
                case scanFacet:
                {
                    if (cmd == "color")
                    {
                        // Optional 'color' entry (after solid)
                        // - continue looking for 'facet'
                        continue;
                    }
                    else if (cmd == "facet")
                    {
                        beginFacet();
                        state_ = scanLoop;  // Next state
                    }
                    else if (cmd == "endsolid")
                    {
                        // Finished with 'endsolid' - find next solid
                        state_ = scanSolid;
                    }
                    else
                    {
                        errMsg_ = perrorParse("facet", cmd);
                    }
                    break;
                }
                case scanLoop:
                {
                    if (cmd == "outer")
                    {
                        // More pedantic would with (tokens_[1] == "loop") too
                        state_ = scanVerts;  // Next state
                    }
                    else
                    {
                        errMsg_ = perrorParse("outer loop", cmd);
                    }
                    break;
                }
                case scanVerts:
                {
                    if (cmd == "vertex")
                    {
                        if (tokens_.size() > 3)
                        {
                            // The tokens are space-delimited and thus okay
                            // for atof()
                            addVertexComponent(tokens_[1].first);
                            addVertexComponent(tokens_[2].first);
                            addVertexComponent(tokens_[3].first);
                        }
                        else
                        {
                            errMsg_ = "Error parsing vertex value";
                        }
                    }
                    else if (cmd == "endloop")
                    {
                        state_ = scanEndFacet;  // Next state
                    }
                    else
                    {
                        errMsg_ = perrorParse("vertex", cmd);
                    }
                    break;
                }
                case scanEndLoop:
                {
                    if (cmd == "endloop")
                    {
                        state_ = scanEndFacet;  // Next state
                    }
                    else
                    {
                        errMsg_ = perrorParse("endloop", cmd);
                    }
                    break;
                }
                case scanEndFacet:
                {
                    if (cmd == "endfacet")
                    {
                        endFacet();
                        state_ = scanFacet;  // Next facet, or endsolid
                    }
                    else
                    {
                        errMsg_ = perrorParse("endfacet", cmd);
                    }
                    break;
                }
                case scanEndSolid:
                {
                    if (cmd == "endsolid")
                    {
                        state_ = scanSolid;  // Start over again
                    }
                    else
                    {
                        errMsg_ = perrorParse("endsolid", cmd);
                    }
                    break;
                }
            }
        }
        while (errMsg_.empty());

        // How much still in the buffer?
        pending = data + gcount - pe;

        if (pending)
        {
            memmove(inbuf, pe, pending);
        }

        if (gcount < buflen)
        {
            break; // done
        }

        if (!errMsg_.empty())
        {
            break;
        }
    }

    if (!errMsg_.empty())
    {
        FatalErrorInFunction
            << errMsg_ << nl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fileFormats::STLReader::readAsciiManual
(
    const fileName& filename
)
{
    IFstream is(filename);
    if (!is)
    {
        FatalErrorInFunction
            << "file " << filename << " not found"
            << exit(FatalError);
    }

    // Create with estimated number of triangles in the STL.
    // 180 bytes / triangle. For simplicity, ignore compression

    const auto fileLen = is.fileSize();

    const label nTrisEstimated =
    (
        (fileLen > 0)
      ? max(label(100), label(fileLen/180))
      : label(100)
    );


    Detail::STLAsciiParseManual lexer(nTrisEstimated);
    lexer.execute(is.stdStream());

    transfer(lexer);

    return true;
}


// ************************************************************************* //
