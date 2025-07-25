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
    Ragel-based parsing of STL ASCII format.
    The goto-based finite state machine (FSM) is generated with

        ragel -G2 -o STLAsciiParseRagel.cc STLAsciiParseRagel.rl

\*---------------------------------------------------------------------------*/

#include "STLAsciiParse.H"
#include "STLReader.H"

#pragma GCC diagnostic ignored "-Wimplicit-fallthrough"
#pragma GCC diagnostic ignored "-Wunused-const-variable"

// https://en.wikipedia.org/wiki/STL_%28file_format%29#ASCII_STL
//
// Format
//
// solid [name]
//
// * where name is an optional string.
// * The file continues with any number of triangles,
//   each represented as follows:
//
// [color ...]
// facet normal ni nj nk
//     outer loop
//         vertex v1x v1y v1z
//         vertex v2x v2y v2z
//         vertex v3x v3y v3z
//     endloop
// endfacet
//
// * where each n or v is a floating-point number.
// * The file concludes with
//
// endsolid [name]

// We take some parsing shortcuts.
// - Ignore 'color' lines
// - Only look for initial 'facet '. Ignore 'normal ...'
// - Ignore name for 'endsolid'
//
// ------------------------------------------------------------------------- //

// Ragel machine definition
// Ragel variables (p, pe, eof, cs, top, stack, ts, te, act) defined later...
//
// Can use 'variable p xxx;' etc to change these names


// Define the machine actions
%%{
    machine stlAscii;

    action line     { ++lineNum_; /* Count '\n' occurrences */ }
    action buffer   { tok = p; /* Local token start */ }

    action bsolid   { beginSolid(word::validate(tok, p)); }
    action bfacet   { beginFacet(); }
    action efacet   { endFacet(); }

    action bvertex  { resetVertex(); }
    action vertexCmpt
    {
        const char saveC = *p;
        *p = '\0';  // Make nul-terminated

        addVertexComponent(tok);
        *p = saveC; // Restore previous character
    }

    ## Error handling. Should be less verbose with ragel-7.
    action err_bsolid { die("solid", p, pe); }
    action err_esolid { die("endsolid", p, pe); }
    action err_bfacet { die("facet", p, pe); }
    action err_efacet { die("endfacet", p, pe); }
    action err_bloop  { die("loop", p, pe); }
    action err_eloop  { die("endloop", p, pe); }
    action err_vertex { die("vertex", p, pe); }

}%%


%%{
    machine stlAscii;

    white   = [ \t\f\r];            # Horizontal whitespace
    nl      = (white* '\n' %line);  # Newline
    dnl     = ([^\n]* '\n' %line);  # Discard up to and including newline
    anyspace = ((space -- '\n') | '\n' $line)*;  # Any space, count newlines

    decimal = ((digit* '.' digit+) | (digit+ '.'?)) ;
    number  = [\-+]? (digit+ | decimal) ([Ee][\-+]? digit+)? ;

    bfacet  = anyspace (/facet/i) white %bfacet dnl;
    efacet  = anyspace (/endfacet/i) %efacet dnl;

    bsolid  = anyspace /solid/i (('' | (white+ [^\n]*)) >buffer nl @bsolid);
    esolid  = anyspace /endsolid/i dnl;

    color   = anyspace /color/i dnl;

    bloop   = anyspace (/outer/i white+ /loop/i) dnl;
    eloop   = anyspace (/endloop/i) dnl;

    vertex  = anyspace /vertex/i
        ((white+ (number > buffer %vertexCmpt)){3} nl);

    main :=
    anyspace
    (
        bsolid $!err_bsolid
        color?
        (
            bfacet $!err_bfacet
            bloop  $!err_bloop
            (vertex $!err_vertex){3}   ## Triangles only
            eloop  $!err_eloop
            efacet $!err_efacet
        )*
        esolid $!err_esolid
    )+ anyspace;

}%%


//
// FSM globals
//

%% write data nofinal;


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Foam
{
namespace Detail
{

//- A lexer for parsing STL ASCII files.
//  Returns DynamicList(s) of points and facets (zoneIds).
//  The facets are within a solid/endsolid grouping
class STLAsciiParseRagel
:
    public Detail::STLAsciiParse
{
    //- Error handling
    void die(const char *what, const char *parsing, const char *pe) const;

public:

    //- From input stream, with the estimated number of triangles in the STL
    STLAsciiParseRagel(const label nTrisEstimated)
    :
        Detail::STLAsciiParse(nTrisEstimated)
    {}

    //- Execute lexer
    void execute(std::istream& is);
};

} // End namespace Detail
} // End namespace Foam


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Length of the input read buffer
#define INBUFLEN 16384

void Foam::Detail::STLAsciiParseRagel::execute(std::istream& is)
{
    if (!is)
    {
        return;
    }

    // cs = code state
    int cs;

    // Initialize FSM variables
    %%{write init;}%%   /* ^^^ FSM initialization here ^^^ */;

    // Local token start
    char *tok = nullptr;

    // Buffering
    char inbuf[INBUFLEN];
    std::streamsize pending = 0;

    // Line-oriented processing loop (as per Ragel pdf example)

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
            break;
        }

        // p,pe = Ragel parsing point, parsing end (default naming)
        // eof  = Ragel EOF point (default naming)

        char *p = inbuf;
        const char *pe = data + gcount;
        const char *eof = nullptr;
        if (!is)
        {
            eof = pe;   // Tag 'pe' as being the EOF for the FSM as well
        }

        // Line-oriented: search backwards to find last newline
        {
            --pe;
            while (*pe != '\n' && pe >= inbuf)
            {
                --pe;
            }
            ++pe;
        }

        %%{write exec;}%%       /* ^^^ FSM execution here ^^^ */;

        if (%%{write error;}%% == cs)
        {
            // FSM failed before finding a token
            FatalErrorInFunction
                << "parse error while scanning near line " << lineNum_ << nl;

            if (p)
            {
                std::string::size_type errLen = (pe - p);
                if (errLen > 80)
                {
                    errLen = 80;
                }

                FatalErrorInFunction
                    << "context: " << std::string(p, errLen) << nl;
            }
            break;
        }

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
    }
}


void Foam::Detail::STLAsciiParseRagel::die
(
    const char *what,
    const char *parsing,
    const char *pe
) const
{
    FatalErrorInFunction
        << nl
        << "Parsing error at or near line " << lineNum_
        <<", while parsing for " << what << nl
        << "    Found text '";

    if (parsing)
    {
        // Output first until newline or 80 chars, or end of parsing content
        for (unsigned i=0; i < 80; ++i)
        {
            if (*parsing == '\n' || parsing == pe) break;
            FatalError << *parsing;
            ++parsing;
        }
    }

    FatalError
        << "'\n"
        << exit(FatalError);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fileFormats::STLReader::readAsciiRagel
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


    Detail::STLAsciiParseRagel lexer(nTrisEstimated);
    lexer.execute(is.stdStream());

    transfer(lexer);

    return true;
}


// ************************************************************************* //
