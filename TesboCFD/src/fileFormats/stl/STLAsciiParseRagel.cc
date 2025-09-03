
#line 1 "STLAsciiParseRagel.rl"
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

#line 107 "STLAsciiParseRagel.rl"




#line 150 "STLAsciiParseRagel.rl"



//
// FSM globals
//


#line 96 "STLAsciiParseRagel.cc"
static const int stlAscii_start = 1;
static const int stlAscii_error = 0;

static const int stlAscii_en_main = 1;


#line 158 "STLAsciiParseRagel.rl"


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
    
#line 156 "STLAsciiParseRagel.cc"
	{
	cs = stlAscii_start;
	}

#line 209 "STLAsciiParseRagel.rl"
   /* ^^^ FSM initialization here ^^^ */;

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

        
#line 215 "STLAsciiParseRagel.cc"
	{
	if ( p == pe )
		goto _test_eof;
	switch ( cs )
	{
tr2:
#line 81 "STLAsciiParseRagel.rl"
	{ ++lineNum_; /* Count '\n' occurrences */ }
	goto st1;
st1:
	if ( ++p == pe )
		goto _test_eof1;
case 1:
#line 229 "STLAsciiParseRagel.cc"
	switch( (*p) ) {
		case 10: goto tr2;
		case 32: goto st1;
		case 83: goto st2;
		case 115: goto st2;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st1;
	goto tr0;
tr0:
#line 99 "STLAsciiParseRagel.rl"
	{ die("solid", p, pe); }
	goto st0;
tr12:
#line 99 "STLAsciiParseRagel.rl"
	{ die("solid", p, pe); }
#line 101 "STLAsciiParseRagel.rl"
	{ die("facet", p, pe); }
#line 100 "STLAsciiParseRagel.rl"
	{ die("endsolid", p, pe); }
	goto st0;
tr17:
#line 101 "STLAsciiParseRagel.rl"
	{ die("facet", p, pe); }
#line 100 "STLAsciiParseRagel.rl"
	{ die("endsolid", p, pe); }
	goto st0;
tr30:
#line 100 "STLAsciiParseRagel.rl"
	{ die("endsolid", p, pe); }
	goto st0;
tr39:
#line 101 "STLAsciiParseRagel.rl"
	{ die("facet", p, pe); }
	goto st0;
tr49:
#line 101 "STLAsciiParseRagel.rl"
	{ die("facet", p, pe); }
#line 103 "STLAsciiParseRagel.rl"
	{ die("loop", p, pe); }
	goto st0;
tr52:
#line 103 "STLAsciiParseRagel.rl"
	{ die("loop", p, pe); }
	goto st0;
tr65:
#line 103 "STLAsciiParseRagel.rl"
	{ die("loop", p, pe); }
#line 105 "STLAsciiParseRagel.rl"
	{ die("vertex", p, pe); }
	goto st0;
tr68:
#line 105 "STLAsciiParseRagel.rl"
	{ die("vertex", p, pe); }
	goto st0;
tr184:
#line 105 "STLAsciiParseRagel.rl"
	{ die("vertex", p, pe); }
#line 104 "STLAsciiParseRagel.rl"
	{ die("endloop", p, pe); }
	goto st0;
tr187:
#line 104 "STLAsciiParseRagel.rl"
	{ die("endloop", p, pe); }
	goto st0;
tr197:
#line 104 "STLAsciiParseRagel.rl"
	{ die("endloop", p, pe); }
#line 102 "STLAsciiParseRagel.rl"
	{ die("endfacet", p, pe); }
	goto st0;
tr200:
#line 102 "STLAsciiParseRagel.rl"
	{ die("endfacet", p, pe); }
	goto st0;
tr214:
#line 101 "STLAsciiParseRagel.rl"
	{ die("facet", p, pe); }
#line 102 "STLAsciiParseRagel.rl"
	{ die("endfacet", p, pe); }
#line 100 "STLAsciiParseRagel.rl"
	{ die("endsolid", p, pe); }
	goto st0;
tr233:
#line 99 "STLAsciiParseRagel.rl"
	{ die("solid", p, pe); }
#line 100 "STLAsciiParseRagel.rl"
	{ die("endsolid", p, pe); }
	goto st0;
#line 319 "STLAsciiParseRagel.cc"
st0:
cs = 0;
	goto _out;
tr235:
#line 81 "STLAsciiParseRagel.rl"
	{ ++lineNum_; /* Count '\n' occurrences */ }
	goto st2;
st2:
	if ( ++p == pe )
		goto _test_eof2;
case 2:
#line 331 "STLAsciiParseRagel.cc"
	switch( (*p) ) {
		case 79: goto st3;
		case 111: goto st3;
	}
	goto tr0;
st3:
	if ( ++p == pe )
		goto _test_eof3;
case 3:
	switch( (*p) ) {
		case 76: goto st4;
		case 108: goto st4;
	}
	goto tr0;
st4:
	if ( ++p == pe )
		goto _test_eof4;
case 4:
	switch( (*p) ) {
		case 73: goto st5;
		case 105: goto st5;
	}
	goto tr0;
st5:
	if ( ++p == pe )
		goto _test_eof5;
case 5:
	switch( (*p) ) {
		case 68: goto st6;
		case 100: goto st6;
	}
	goto tr0;
st6:
	if ( ++p == pe )
		goto _test_eof6;
case 6:
	switch( (*p) ) {
		case 9: goto tr8;
		case 10: goto tr9;
		case 32: goto tr8;
	}
	if ( 12 <= (*p) && (*p) <= 13 )
		goto tr8;
	goto tr0;
tr8:
#line 82 "STLAsciiParseRagel.rl"
	{ tok = p; /* Local token start */ }
	goto st7;
st7:
	if ( ++p == pe )
		goto _test_eof7;
case 7:
#line 384 "STLAsciiParseRagel.cc"
	if ( (*p) == 10 )
		goto tr11;
	goto st7;
tr9:
#line 82 "STLAsciiParseRagel.rl"
	{ tok = p; /* Local token start */ }
#line 84 "STLAsciiParseRagel.rl"
	{ beginSolid(word::validate(tok, p)); }
	goto st8;
tr11:
#line 84 "STLAsciiParseRagel.rl"
	{ beginSolid(word::validate(tok, p)); }
	goto st8;
st8:
	if ( ++p == pe )
		goto _test_eof8;
case 8:
#line 402 "STLAsciiParseRagel.cc"
	switch( (*p) ) {
		case 32: goto tr13;
		case 67: goto tr14;
		case 69: goto tr15;
		case 70: goto tr16;
		case 99: goto tr14;
		case 101: goto tr15;
		case 102: goto tr16;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto tr13;
	goto tr12;
tr13:
#line 81 "STLAsciiParseRagel.rl"
	{ ++lineNum_; /* Count '\n' occurrences */ }
	goto st9;
st9:
	if ( ++p == pe )
		goto _test_eof9;
case 9:
#line 423 "STLAsciiParseRagel.cc"
	switch( (*p) ) {
		case 10: goto tr13;
		case 32: goto st9;
		case 67: goto st10;
		case 69: goto st17;
		case 70: goto st25;
		case 99: goto st10;
		case 101: goto st17;
		case 102: goto st25;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st9;
	goto tr17;
tr14:
#line 81 "STLAsciiParseRagel.rl"
	{ ++lineNum_; /* Count '\n' occurrences */ }
	goto st10;
st10:
	if ( ++p == pe )
		goto _test_eof10;
case 10:
#line 445 "STLAsciiParseRagel.cc"
	switch( (*p) ) {
		case 79: goto st11;
		case 111: goto st11;
	}
	goto st0;
st11:
	if ( ++p == pe )
		goto _test_eof11;
case 11:
	switch( (*p) ) {
		case 76: goto st12;
		case 108: goto st12;
	}
	goto st0;
st12:
	if ( ++p == pe )
		goto _test_eof12;
case 12:
	switch( (*p) ) {
		case 79: goto st13;
		case 111: goto st13;
	}
	goto st0;
st13:
	if ( ++p == pe )
		goto _test_eof13;
case 13:
	switch( (*p) ) {
		case 82: goto st14;
		case 114: goto st14;
	}
	goto st0;
st14:
	if ( ++p == pe )
		goto _test_eof14;
case 14:
	if ( (*p) == 10 )
		goto st15;
	goto st14;
st15:
	if ( ++p == pe )
		goto _test_eof15;
case 15:
	switch( (*p) ) {
		case 32: goto tr28;
		case 69: goto tr15;
		case 70: goto tr16;
		case 101: goto tr15;
		case 102: goto tr16;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto tr28;
	goto tr17;
tr28:
#line 81 "STLAsciiParseRagel.rl"
	{ ++lineNum_; /* Count '\n' occurrences */ }
	goto st16;
st16:
	if ( ++p == pe )
		goto _test_eof16;
case 16:
#line 507 "STLAsciiParseRagel.cc"
	switch( (*p) ) {
		case 10: goto tr28;
		case 32: goto st16;
		case 69: goto st17;
		case 70: goto st25;
		case 101: goto st17;
		case 102: goto st25;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st16;
	goto tr17;
tr15:
#line 81 "STLAsciiParseRagel.rl"
	{ ++lineNum_; /* Count '\n' occurrences */ }
	goto st17;
st17:
	if ( ++p == pe )
		goto _test_eof17;
case 17:
#line 527 "STLAsciiParseRagel.cc"
	switch( (*p) ) {
		case 78: goto st18;
		case 110: goto st18;
	}
	goto tr30;
st18:
	if ( ++p == pe )
		goto _test_eof18;
case 18:
	switch( (*p) ) {
		case 68: goto st19;
		case 100: goto st19;
	}
	goto tr30;
st19:
	if ( ++p == pe )
		goto _test_eof19;
case 19:
	switch( (*p) ) {
		case 83: goto st20;
		case 115: goto st20;
	}
	goto tr30;
st20:
	if ( ++p == pe )
		goto _test_eof20;
case 20:
	switch( (*p) ) {
		case 79: goto st21;
		case 111: goto st21;
	}
	goto tr30;
st21:
	if ( ++p == pe )
		goto _test_eof21;
case 21:
	switch( (*p) ) {
		case 76: goto st22;
		case 108: goto st22;
	}
	goto tr30;
st22:
	if ( ++p == pe )
		goto _test_eof22;
case 22:
	switch( (*p) ) {
		case 73: goto st23;
		case 105: goto st23;
	}
	goto tr30;
st23:
	if ( ++p == pe )
		goto _test_eof23;
case 23:
	switch( (*p) ) {
		case 68: goto st24;
		case 100: goto st24;
	}
	goto tr30;
st24:
	if ( ++p == pe )
		goto _test_eof24;
case 24:
	if ( (*p) == 10 )
		goto st164;
	goto st24;
st164:
	if ( ++p == pe )
		goto _test_eof164;
case 164:
	switch( (*p) ) {
		case 32: goto tr234;
		case 83: goto tr235;
		case 115: goto tr235;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto tr234;
	goto tr233;
tr234:
#line 81 "STLAsciiParseRagel.rl"
	{ ++lineNum_; /* Count '\n' occurrences */ }
	goto st165;
st165:
	if ( ++p == pe )
		goto _test_eof165;
case 165:
#line 614 "STLAsciiParseRagel.cc"
	switch( (*p) ) {
		case 10: goto tr234;
		case 32: goto st165;
		case 83: goto st2;
		case 115: goto st2;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st165;
	goto tr0;
tr16:
#line 81 "STLAsciiParseRagel.rl"
	{ ++lineNum_; /* Count '\n' occurrences */ }
	goto st25;
st25:
	if ( ++p == pe )
		goto _test_eof25;
case 25:
#line 632 "STLAsciiParseRagel.cc"
	switch( (*p) ) {
		case 65: goto st26;
		case 97: goto st26;
	}
	goto tr39;
st26:
	if ( ++p == pe )
		goto _test_eof26;
case 26:
	switch( (*p) ) {
		case 67: goto st27;
		case 99: goto st27;
	}
	goto tr39;
st27:
	if ( ++p == pe )
		goto _test_eof27;
case 27:
	switch( (*p) ) {
		case 69: goto st28;
		case 101: goto st28;
	}
	goto tr39;
st28:
	if ( ++p == pe )
		goto _test_eof28;
case 28:
	switch( (*p) ) {
		case 84: goto st29;
		case 116: goto st29;
	}
	goto tr39;
st29:
	if ( ++p == pe )
		goto _test_eof29;
case 29:
	switch( (*p) ) {
		case 9: goto st30;
		case 32: goto st30;
	}
	if ( 12 <= (*p) && (*p) <= 13 )
		goto st30;
	goto tr39;
st30:
	if ( ++p == pe )
		goto _test_eof30;
case 30:
	if ( (*p) == 10 )
		goto tr46;
	goto tr45;
tr45:
#line 85 "STLAsciiParseRagel.rl"
	{ beginFacet(); }
	goto st31;
st31:
	if ( ++p == pe )
		goto _test_eof31;
case 31:
#line 691 "STLAsciiParseRagel.cc"
	if ( (*p) == 10 )
		goto st32;
	goto st31;
tr46:
#line 85 "STLAsciiParseRagel.rl"
	{ beginFacet(); }
	goto st32;
st32:
	if ( ++p == pe )
		goto _test_eof32;
case 32:
#line 703 "STLAsciiParseRagel.cc"
	switch( (*p) ) {
		case 32: goto tr50;
		case 79: goto tr51;
		case 111: goto tr51;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto tr50;
	goto tr49;
tr50:
#line 81 "STLAsciiParseRagel.rl"
	{ ++lineNum_; /* Count '\n' occurrences */ }
	goto st33;
st33:
	if ( ++p == pe )
		goto _test_eof33;
case 33:
#line 720 "STLAsciiParseRagel.cc"
	switch( (*p) ) {
		case 10: goto tr50;
		case 32: goto st33;
		case 79: goto st34;
		case 111: goto st34;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st33;
	goto tr52;
tr51:
#line 81 "STLAsciiParseRagel.rl"
	{ ++lineNum_; /* Count '\n' occurrences */ }
	goto st34;
st34:
	if ( ++p == pe )
		goto _test_eof34;
case 34:
#line 738 "STLAsciiParseRagel.cc"
	switch( (*p) ) {
		case 85: goto st35;
		case 117: goto st35;
	}
	goto tr52;
st35:
	if ( ++p == pe )
		goto _test_eof35;
case 35:
	switch( (*p) ) {
		case 84: goto st36;
		case 116: goto st36;
	}
	goto tr52;
st36:
	if ( ++p == pe )
		goto _test_eof36;
case 36:
	switch( (*p) ) {
		case 69: goto st37;
		case 101: goto st37;
	}
	goto tr52;
st37:
	if ( ++p == pe )
		goto _test_eof37;
case 37:
	switch( (*p) ) {
		case 82: goto st38;
		case 114: goto st38;
	}
	goto tr52;
st38:
	if ( ++p == pe )
		goto _test_eof38;
case 38:
	switch( (*p) ) {
		case 9: goto st39;
		case 32: goto st39;
	}
	if ( 12 <= (*p) && (*p) <= 13 )
		goto st39;
	goto tr52;
st39:
	if ( ++p == pe )
		goto _test_eof39;
case 39:
	switch( (*p) ) {
		case 9: goto st39;
		case 32: goto st39;
		case 76: goto st40;
		case 108: goto st40;
	}
	if ( 12 <= (*p) && (*p) <= 13 )
		goto st39;
	goto tr52;
st40:
	if ( ++p == pe )
		goto _test_eof40;
case 40:
	switch( (*p) ) {
		case 79: goto st41;
		case 111: goto st41;
	}
	goto tr52;
st41:
	if ( ++p == pe )
		goto _test_eof41;
case 41:
	switch( (*p) ) {
		case 79: goto st42;
		case 111: goto st42;
	}
	goto tr52;
st42:
	if ( ++p == pe )
		goto _test_eof42;
case 42:
	switch( (*p) ) {
		case 80: goto st43;
		case 112: goto st43;
	}
	goto tr52;
st43:
	if ( ++p == pe )
		goto _test_eof43;
case 43:
	if ( (*p) == 10 )
		goto st44;
	goto st43;
st44:
	if ( ++p == pe )
		goto _test_eof44;
case 44:
	switch( (*p) ) {
		case 32: goto tr66;
		case 86: goto tr67;
		case 118: goto tr67;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto tr66;
	goto tr65;
tr66:
#line 81 "STLAsciiParseRagel.rl"
	{ ++lineNum_; /* Count '\n' occurrences */ }
	goto st45;
st45:
	if ( ++p == pe )
		goto _test_eof45;
case 45:
#line 849 "STLAsciiParseRagel.cc"
	switch( (*p) ) {
		case 10: goto tr66;
		case 32: goto st45;
		case 86: goto st46;
		case 118: goto st46;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st45;
	goto tr68;
tr67:
#line 81 "STLAsciiParseRagel.rl"
	{ ++lineNum_; /* Count '\n' occurrences */ }
	goto st46;
st46:
	if ( ++p == pe )
		goto _test_eof46;
case 46:
#line 867 "STLAsciiParseRagel.cc"
	switch( (*p) ) {
		case 69: goto st47;
		case 101: goto st47;
	}
	goto tr68;
st47:
	if ( ++p == pe )
		goto _test_eof47;
case 47:
	switch( (*p) ) {
		case 82: goto st48;
		case 114: goto st48;
	}
	goto tr68;
st48:
	if ( ++p == pe )
		goto _test_eof48;
case 48:
	switch( (*p) ) {
		case 84: goto st49;
		case 116: goto st49;
	}
	goto tr68;
st49:
	if ( ++p == pe )
		goto _test_eof49;
case 49:
	switch( (*p) ) {
		case 69: goto st50;
		case 101: goto st50;
	}
	goto tr68;
st50:
	if ( ++p == pe )
		goto _test_eof50;
case 50:
	switch( (*p) ) {
		case 88: goto st51;
		case 120: goto st51;
	}
	goto tr68;
st51:
	if ( ++p == pe )
		goto _test_eof51;
case 51:
	switch( (*p) ) {
		case 9: goto st52;
		case 32: goto st52;
	}
	if ( 12 <= (*p) && (*p) <= 13 )
		goto st52;
	goto tr68;
st52:
	if ( ++p == pe )
		goto _test_eof52;
case 52:
	switch( (*p) ) {
		case 9: goto st52;
		case 32: goto st52;
		case 43: goto tr77;
		case 45: goto tr77;
		case 46: goto tr78;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr79;
	} else if ( (*p) >= 12 )
		goto st52;
	goto tr68;
tr77:
#line 82 "STLAsciiParseRagel.rl"
	{ tok = p; /* Local token start */ }
	goto st53;
st53:
	if ( ++p == pe )
		goto _test_eof53;
case 53:
#line 945 "STLAsciiParseRagel.cc"
	if ( (*p) == 46 )
		goto st54;
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st163;
	goto tr68;
tr78:
#line 82 "STLAsciiParseRagel.rl"
	{ tok = p; /* Local token start */ }
	goto st54;
st54:
	if ( ++p == pe )
		goto _test_eof54;
case 54:
#line 959 "STLAsciiParseRagel.cc"
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st55;
	goto tr68;
st55:
	if ( ++p == pe )
		goto _test_eof55;
case 55:
	switch( (*p) ) {
		case 9: goto tr83;
		case 32: goto tr83;
		case 69: goto st160;
		case 101: goto st160;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st55;
	} else if ( (*p) >= 12 )
		goto tr83;
	goto tr68;
tr83:
#line 90 "STLAsciiParseRagel.rl"
	{
        const char saveC = *p;
        *p = '\0';  // Make nul-terminated

        addVertexComponent(tok);
        *p = saveC; // Restore previous character
    }
	goto st56;
st56:
	if ( ++p == pe )
		goto _test_eof56;
case 56:
#line 993 "STLAsciiParseRagel.cc"
	switch( (*p) ) {
		case 9: goto st56;
		case 32: goto st56;
		case 43: goto tr86;
		case 45: goto tr86;
		case 46: goto tr87;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr88;
	} else if ( (*p) >= 12 )
		goto st56;
	goto tr68;
tr86:
#line 82 "STLAsciiParseRagel.rl"
	{ tok = p; /* Local token start */ }
	goto st57;
st57:
	if ( ++p == pe )
		goto _test_eof57;
case 57:
#line 1015 "STLAsciiParseRagel.cc"
	if ( (*p) == 46 )
		goto st58;
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st159;
	goto tr68;
tr87:
#line 82 "STLAsciiParseRagel.rl"
	{ tok = p; /* Local token start */ }
	goto st58;
st58:
	if ( ++p == pe )
		goto _test_eof58;
case 58:
#line 1029 "STLAsciiParseRagel.cc"
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st59;
	goto tr68;
st59:
	if ( ++p == pe )
		goto _test_eof59;
case 59:
	switch( (*p) ) {
		case 9: goto tr92;
		case 32: goto tr92;
		case 69: goto st156;
		case 101: goto st156;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st59;
	} else if ( (*p) >= 12 )
		goto tr92;
	goto tr68;
tr92:
#line 90 "STLAsciiParseRagel.rl"
	{
        const char saveC = *p;
        *p = '\0';  // Make nul-terminated

        addVertexComponent(tok);
        *p = saveC; // Restore previous character
    }
	goto st60;
st60:
	if ( ++p == pe )
		goto _test_eof60;
case 60:
#line 1063 "STLAsciiParseRagel.cc"
	switch( (*p) ) {
		case 9: goto st60;
		case 32: goto st60;
		case 43: goto tr95;
		case 45: goto tr95;
		case 46: goto tr96;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr97;
	} else if ( (*p) >= 12 )
		goto st60;
	goto tr68;
tr95:
#line 82 "STLAsciiParseRagel.rl"
	{ tok = p; /* Local token start */ }
	goto st61;
st61:
	if ( ++p == pe )
		goto _test_eof61;
case 61:
#line 1085 "STLAsciiParseRagel.cc"
	if ( (*p) == 46 )
		goto st62;
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st155;
	goto tr68;
tr96:
#line 82 "STLAsciiParseRagel.rl"
	{ tok = p; /* Local token start */ }
	goto st62;
st62:
	if ( ++p == pe )
		goto _test_eof62;
case 62:
#line 1099 "STLAsciiParseRagel.cc"
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st63;
	goto tr68;
st63:
	if ( ++p == pe )
		goto _test_eof63;
case 63:
	switch( (*p) ) {
		case 9: goto tr101;
		case 10: goto tr102;
		case 32: goto tr101;
		case 69: goto st152;
		case 101: goto st152;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st63;
	} else if ( (*p) >= 12 )
		goto tr101;
	goto tr68;
tr101:
#line 90 "STLAsciiParseRagel.rl"
	{
        const char saveC = *p;
        *p = '\0';  // Make nul-terminated

        addVertexComponent(tok);
        *p = saveC; // Restore previous character
    }
	goto st64;
st64:
	if ( ++p == pe )
		goto _test_eof64;
case 64:
#line 1134 "STLAsciiParseRagel.cc"
	switch( (*p) ) {
		case 9: goto st64;
		case 10: goto st65;
		case 32: goto st64;
	}
	if ( 12 <= (*p) && (*p) <= 13 )
		goto st64;
	goto tr68;
tr102:
#line 90 "STLAsciiParseRagel.rl"
	{
        const char saveC = *p;
        *p = '\0';  // Make nul-terminated

        addVertexComponent(tok);
        *p = saveC; // Restore previous character
    }
	goto st65;
st65:
	if ( ++p == pe )
		goto _test_eof65;
case 65:
#line 1157 "STLAsciiParseRagel.cc"
	switch( (*p) ) {
		case 32: goto tr106;
		case 86: goto tr107;
		case 118: goto tr107;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto tr106;
	goto tr68;
tr106:
#line 81 "STLAsciiParseRagel.rl"
	{ ++lineNum_; /* Count '\n' occurrences */ }
	goto st66;
st66:
	if ( ++p == pe )
		goto _test_eof66;
case 66:
#line 1174 "STLAsciiParseRagel.cc"
	switch( (*p) ) {
		case 10: goto tr106;
		case 32: goto st66;
		case 86: goto st67;
		case 118: goto st67;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st66;
	goto tr68;
tr107:
#line 81 "STLAsciiParseRagel.rl"
	{ ++lineNum_; /* Count '\n' occurrences */ }
	goto st67;
st67:
	if ( ++p == pe )
		goto _test_eof67;
case 67:
#line 1192 "STLAsciiParseRagel.cc"
	switch( (*p) ) {
		case 69: goto st68;
		case 101: goto st68;
	}
	goto tr68;
st68:
	if ( ++p == pe )
		goto _test_eof68;
case 68:
	switch( (*p) ) {
		case 82: goto st69;
		case 114: goto st69;
	}
	goto tr68;
st69:
	if ( ++p == pe )
		goto _test_eof69;
case 69:
	switch( (*p) ) {
		case 84: goto st70;
		case 116: goto st70;
	}
	goto tr68;
st70:
	if ( ++p == pe )
		goto _test_eof70;
case 70:
	switch( (*p) ) {
		case 69: goto st71;
		case 101: goto st71;
	}
	goto tr68;
st71:
	if ( ++p == pe )
		goto _test_eof71;
case 71:
	switch( (*p) ) {
		case 88: goto st72;
		case 120: goto st72;
	}
	goto tr68;
st72:
	if ( ++p == pe )
		goto _test_eof72;
case 72:
	switch( (*p) ) {
		case 9: goto st73;
		case 32: goto st73;
	}
	if ( 12 <= (*p) && (*p) <= 13 )
		goto st73;
	goto tr68;
st73:
	if ( ++p == pe )
		goto _test_eof73;
case 73:
	switch( (*p) ) {
		case 9: goto st73;
		case 32: goto st73;
		case 43: goto tr116;
		case 45: goto tr116;
		case 46: goto tr117;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr118;
	} else if ( (*p) >= 12 )
		goto st73;
	goto tr68;
tr116:
#line 82 "STLAsciiParseRagel.rl"
	{ tok = p; /* Local token start */ }
	goto st74;
st74:
	if ( ++p == pe )
		goto _test_eof74;
case 74:
#line 1270 "STLAsciiParseRagel.cc"
	if ( (*p) == 46 )
		goto st75;
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st151;
	goto tr68;
tr117:
#line 82 "STLAsciiParseRagel.rl"
	{ tok = p; /* Local token start */ }
	goto st75;
st75:
	if ( ++p == pe )
		goto _test_eof75;
case 75:
#line 1284 "STLAsciiParseRagel.cc"
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st76;
	goto tr68;
st76:
	if ( ++p == pe )
		goto _test_eof76;
case 76:
	switch( (*p) ) {
		case 9: goto tr122;
		case 32: goto tr122;
		case 69: goto st148;
		case 101: goto st148;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st76;
	} else if ( (*p) >= 12 )
		goto tr122;
	goto tr68;
tr122:
#line 90 "STLAsciiParseRagel.rl"
	{
        const char saveC = *p;
        *p = '\0';  // Make nul-terminated

        addVertexComponent(tok);
        *p = saveC; // Restore previous character
    }
	goto st77;
st77:
	if ( ++p == pe )
		goto _test_eof77;
case 77:
#line 1318 "STLAsciiParseRagel.cc"
	switch( (*p) ) {
		case 9: goto st77;
		case 32: goto st77;
		case 43: goto tr125;
		case 45: goto tr125;
		case 46: goto tr126;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr127;
	} else if ( (*p) >= 12 )
		goto st77;
	goto tr68;
tr125:
#line 82 "STLAsciiParseRagel.rl"
	{ tok = p; /* Local token start */ }
	goto st78;
st78:
	if ( ++p == pe )
		goto _test_eof78;
case 78:
#line 1340 "STLAsciiParseRagel.cc"
	if ( (*p) == 46 )
		goto st79;
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st147;
	goto tr68;
tr126:
#line 82 "STLAsciiParseRagel.rl"
	{ tok = p; /* Local token start */ }
	goto st79;
st79:
	if ( ++p == pe )
		goto _test_eof79;
case 79:
#line 1354 "STLAsciiParseRagel.cc"
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st80;
	goto tr68;
st80:
	if ( ++p == pe )
		goto _test_eof80;
case 80:
	switch( (*p) ) {
		case 9: goto tr131;
		case 32: goto tr131;
		case 69: goto st144;
		case 101: goto st144;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st80;
	} else if ( (*p) >= 12 )
		goto tr131;
	goto tr68;
tr131:
#line 90 "STLAsciiParseRagel.rl"
	{
        const char saveC = *p;
        *p = '\0';  // Make nul-terminated

        addVertexComponent(tok);
        *p = saveC; // Restore previous character
    }
	goto st81;
st81:
	if ( ++p == pe )
		goto _test_eof81;
case 81:
#line 1388 "STLAsciiParseRagel.cc"
	switch( (*p) ) {
		case 9: goto st81;
		case 32: goto st81;
		case 43: goto tr134;
		case 45: goto tr134;
		case 46: goto tr135;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr136;
	} else if ( (*p) >= 12 )
		goto st81;
	goto tr68;
tr134:
#line 82 "STLAsciiParseRagel.rl"
	{ tok = p; /* Local token start */ }
	goto st82;
st82:
	if ( ++p == pe )
		goto _test_eof82;
case 82:
#line 1410 "STLAsciiParseRagel.cc"
	if ( (*p) == 46 )
		goto st83;
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st143;
	goto tr68;
tr135:
#line 82 "STLAsciiParseRagel.rl"
	{ tok = p; /* Local token start */ }
	goto st83;
st83:
	if ( ++p == pe )
		goto _test_eof83;
case 83:
#line 1424 "STLAsciiParseRagel.cc"
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st84;
	goto tr68;
st84:
	if ( ++p == pe )
		goto _test_eof84;
case 84:
	switch( (*p) ) {
		case 9: goto tr140;
		case 10: goto tr141;
		case 32: goto tr140;
		case 69: goto st140;
		case 101: goto st140;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st84;
	} else if ( (*p) >= 12 )
		goto tr140;
	goto tr68;
tr140:
#line 90 "STLAsciiParseRagel.rl"
	{
        const char saveC = *p;
        *p = '\0';  // Make nul-terminated

        addVertexComponent(tok);
        *p = saveC; // Restore previous character
    }
	goto st85;
st85:
	if ( ++p == pe )
		goto _test_eof85;
case 85:
#line 1459 "STLAsciiParseRagel.cc"
	switch( (*p) ) {
		case 9: goto st85;
		case 10: goto st86;
		case 32: goto st85;
	}
	if ( 12 <= (*p) && (*p) <= 13 )
		goto st85;
	goto tr68;
tr141:
#line 90 "STLAsciiParseRagel.rl"
	{
        const char saveC = *p;
        *p = '\0';  // Make nul-terminated

        addVertexComponent(tok);
        *p = saveC; // Restore previous character
    }
	goto st86;
st86:
	if ( ++p == pe )
		goto _test_eof86;
case 86:
#line 1482 "STLAsciiParseRagel.cc"
	switch( (*p) ) {
		case 32: goto tr145;
		case 86: goto tr146;
		case 118: goto tr146;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto tr145;
	goto tr68;
tr145:
#line 81 "STLAsciiParseRagel.rl"
	{ ++lineNum_; /* Count '\n' occurrences */ }
	goto st87;
st87:
	if ( ++p == pe )
		goto _test_eof87;
case 87:
#line 1499 "STLAsciiParseRagel.cc"
	switch( (*p) ) {
		case 10: goto tr145;
		case 32: goto st87;
		case 86: goto st88;
		case 118: goto st88;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st87;
	goto tr68;
tr146:
#line 81 "STLAsciiParseRagel.rl"
	{ ++lineNum_; /* Count '\n' occurrences */ }
	goto st88;
st88:
	if ( ++p == pe )
		goto _test_eof88;
case 88:
#line 1517 "STLAsciiParseRagel.cc"
	switch( (*p) ) {
		case 69: goto st89;
		case 101: goto st89;
	}
	goto tr68;
st89:
	if ( ++p == pe )
		goto _test_eof89;
case 89:
	switch( (*p) ) {
		case 82: goto st90;
		case 114: goto st90;
	}
	goto tr68;
st90:
	if ( ++p == pe )
		goto _test_eof90;
case 90:
	switch( (*p) ) {
		case 84: goto st91;
		case 116: goto st91;
	}
	goto tr68;
st91:
	if ( ++p == pe )
		goto _test_eof91;
case 91:
	switch( (*p) ) {
		case 69: goto st92;
		case 101: goto st92;
	}
	goto tr68;
st92:
	if ( ++p == pe )
		goto _test_eof92;
case 92:
	switch( (*p) ) {
		case 88: goto st93;
		case 120: goto st93;
	}
	goto tr68;
st93:
	if ( ++p == pe )
		goto _test_eof93;
case 93:
	switch( (*p) ) {
		case 9: goto st94;
		case 32: goto st94;
	}
	if ( 12 <= (*p) && (*p) <= 13 )
		goto st94;
	goto tr68;
st94:
	if ( ++p == pe )
		goto _test_eof94;
case 94:
	switch( (*p) ) {
		case 9: goto st94;
		case 32: goto st94;
		case 43: goto tr155;
		case 45: goto tr155;
		case 46: goto tr156;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr157;
	} else if ( (*p) >= 12 )
		goto st94;
	goto tr68;
tr155:
#line 82 "STLAsciiParseRagel.rl"
	{ tok = p; /* Local token start */ }
	goto st95;
st95:
	if ( ++p == pe )
		goto _test_eof95;
case 95:
#line 1595 "STLAsciiParseRagel.cc"
	if ( (*p) == 46 )
		goto st96;
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st139;
	goto tr68;
tr156:
#line 82 "STLAsciiParseRagel.rl"
	{ tok = p; /* Local token start */ }
	goto st96;
st96:
	if ( ++p == pe )
		goto _test_eof96;
case 96:
#line 1609 "STLAsciiParseRagel.cc"
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st97;
	goto tr68;
st97:
	if ( ++p == pe )
		goto _test_eof97;
case 97:
	switch( (*p) ) {
		case 9: goto tr161;
		case 32: goto tr161;
		case 69: goto st136;
		case 101: goto st136;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st97;
	} else if ( (*p) >= 12 )
		goto tr161;
	goto tr68;
tr161:
#line 90 "STLAsciiParseRagel.rl"
	{
        const char saveC = *p;
        *p = '\0';  // Make nul-terminated

        addVertexComponent(tok);
        *p = saveC; // Restore previous character
    }
	goto st98;
st98:
	if ( ++p == pe )
		goto _test_eof98;
case 98:
#line 1643 "STLAsciiParseRagel.cc"
	switch( (*p) ) {
		case 9: goto st98;
		case 32: goto st98;
		case 43: goto tr164;
		case 45: goto tr164;
		case 46: goto tr165;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr166;
	} else if ( (*p) >= 12 )
		goto st98;
	goto tr68;
tr164:
#line 82 "STLAsciiParseRagel.rl"
	{ tok = p; /* Local token start */ }
	goto st99;
st99:
	if ( ++p == pe )
		goto _test_eof99;
case 99:
#line 1665 "STLAsciiParseRagel.cc"
	if ( (*p) == 46 )
		goto st100;
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st135;
	goto tr68;
tr165:
#line 82 "STLAsciiParseRagel.rl"
	{ tok = p; /* Local token start */ }
	goto st100;
st100:
	if ( ++p == pe )
		goto _test_eof100;
case 100:
#line 1679 "STLAsciiParseRagel.cc"
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st101;
	goto tr68;
st101:
	if ( ++p == pe )
		goto _test_eof101;
case 101:
	switch( (*p) ) {
		case 9: goto tr170;
		case 32: goto tr170;
		case 69: goto st132;
		case 101: goto st132;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st101;
	} else if ( (*p) >= 12 )
		goto tr170;
	goto tr68;
tr170:
#line 90 "STLAsciiParseRagel.rl"
	{
        const char saveC = *p;
        *p = '\0';  // Make nul-terminated

        addVertexComponent(tok);
        *p = saveC; // Restore previous character
    }
	goto st102;
st102:
	if ( ++p == pe )
		goto _test_eof102;
case 102:
#line 1713 "STLAsciiParseRagel.cc"
	switch( (*p) ) {
		case 9: goto st102;
		case 32: goto st102;
		case 43: goto tr173;
		case 45: goto tr173;
		case 46: goto tr174;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr175;
	} else if ( (*p) >= 12 )
		goto st102;
	goto tr68;
tr173:
#line 82 "STLAsciiParseRagel.rl"
	{ tok = p; /* Local token start */ }
	goto st103;
st103:
	if ( ++p == pe )
		goto _test_eof103;
case 103:
#line 1735 "STLAsciiParseRagel.cc"
	if ( (*p) == 46 )
		goto st104;
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st131;
	goto tr68;
tr174:
#line 82 "STLAsciiParseRagel.rl"
	{ tok = p; /* Local token start */ }
	goto st104;
st104:
	if ( ++p == pe )
		goto _test_eof104;
case 104:
#line 1749 "STLAsciiParseRagel.cc"
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st105;
	goto tr68;
st105:
	if ( ++p == pe )
		goto _test_eof105;
case 105:
	switch( (*p) ) {
		case 9: goto tr179;
		case 10: goto tr180;
		case 32: goto tr179;
		case 69: goto st128;
		case 101: goto st128;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st105;
	} else if ( (*p) >= 12 )
		goto tr179;
	goto tr68;
tr179:
#line 90 "STLAsciiParseRagel.rl"
	{
        const char saveC = *p;
        *p = '\0';  // Make nul-terminated

        addVertexComponent(tok);
        *p = saveC; // Restore previous character
    }
	goto st106;
st106:
	if ( ++p == pe )
		goto _test_eof106;
case 106:
#line 1784 "STLAsciiParseRagel.cc"
	switch( (*p) ) {
		case 9: goto st106;
		case 10: goto st107;
		case 32: goto st106;
	}
	if ( 12 <= (*p) && (*p) <= 13 )
		goto st106;
	goto tr68;
tr180:
#line 90 "STLAsciiParseRagel.rl"
	{
        const char saveC = *p;
        *p = '\0';  // Make nul-terminated

        addVertexComponent(tok);
        *p = saveC; // Restore previous character
    }
	goto st107;
st107:
	if ( ++p == pe )
		goto _test_eof107;
case 107:
#line 1807 "STLAsciiParseRagel.cc"
	switch( (*p) ) {
		case 32: goto tr185;
		case 69: goto tr186;
		case 101: goto tr186;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto tr185;
	goto tr184;
tr185:
#line 81 "STLAsciiParseRagel.rl"
	{ ++lineNum_; /* Count '\n' occurrences */ }
	goto st108;
st108:
	if ( ++p == pe )
		goto _test_eof108;
case 108:
#line 1824 "STLAsciiParseRagel.cc"
	switch( (*p) ) {
		case 10: goto tr185;
		case 32: goto st108;
		case 69: goto st109;
		case 101: goto st109;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st108;
	goto tr187;
tr186:
#line 81 "STLAsciiParseRagel.rl"
	{ ++lineNum_; /* Count '\n' occurrences */ }
	goto st109;
st109:
	if ( ++p == pe )
		goto _test_eof109;
case 109:
#line 1842 "STLAsciiParseRagel.cc"
	switch( (*p) ) {
		case 78: goto st110;
		case 110: goto st110;
	}
	goto tr187;
st110:
	if ( ++p == pe )
		goto _test_eof110;
case 110:
	switch( (*p) ) {
		case 68: goto st111;
		case 100: goto st111;
	}
	goto tr187;
st111:
	if ( ++p == pe )
		goto _test_eof111;
case 111:
	switch( (*p) ) {
		case 76: goto st112;
		case 108: goto st112;
	}
	goto tr187;
st112:
	if ( ++p == pe )
		goto _test_eof112;
case 112:
	switch( (*p) ) {
		case 79: goto st113;
		case 111: goto st113;
	}
	goto tr187;
st113:
	if ( ++p == pe )
		goto _test_eof113;
case 113:
	switch( (*p) ) {
		case 79: goto st114;
		case 111: goto st114;
	}
	goto tr187;
st114:
	if ( ++p == pe )
		goto _test_eof114;
case 114:
	switch( (*p) ) {
		case 80: goto st115;
		case 112: goto st115;
	}
	goto tr187;
st115:
	if ( ++p == pe )
		goto _test_eof115;
case 115:
	if ( (*p) == 10 )
		goto st116;
	goto st115;
st116:
	if ( ++p == pe )
		goto _test_eof116;
case 116:
	switch( (*p) ) {
		case 32: goto tr198;
		case 69: goto tr199;
		case 101: goto tr199;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto tr198;
	goto tr197;
tr198:
#line 81 "STLAsciiParseRagel.rl"
	{ ++lineNum_; /* Count '\n' occurrences */ }
	goto st117;
st117:
	if ( ++p == pe )
		goto _test_eof117;
case 117:
#line 1920 "STLAsciiParseRagel.cc"
	switch( (*p) ) {
		case 10: goto tr198;
		case 32: goto st117;
		case 69: goto st118;
		case 101: goto st118;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st117;
	goto tr200;
tr199:
#line 81 "STLAsciiParseRagel.rl"
	{ ++lineNum_; /* Count '\n' occurrences */ }
	goto st118;
st118:
	if ( ++p == pe )
		goto _test_eof118;
case 118:
#line 1938 "STLAsciiParseRagel.cc"
	switch( (*p) ) {
		case 78: goto st119;
		case 110: goto st119;
	}
	goto tr200;
st119:
	if ( ++p == pe )
		goto _test_eof119;
case 119:
	switch( (*p) ) {
		case 68: goto st120;
		case 100: goto st120;
	}
	goto tr200;
st120:
	if ( ++p == pe )
		goto _test_eof120;
case 120:
	switch( (*p) ) {
		case 70: goto st121;
		case 102: goto st121;
	}
	goto tr200;
st121:
	if ( ++p == pe )
		goto _test_eof121;
case 121:
	switch( (*p) ) {
		case 65: goto st122;
		case 97: goto st122;
	}
	goto tr200;
st122:
	if ( ++p == pe )
		goto _test_eof122;
case 122:
	switch( (*p) ) {
		case 67: goto st123;
		case 99: goto st123;
	}
	goto tr200;
st123:
	if ( ++p == pe )
		goto _test_eof123;
case 123:
	switch( (*p) ) {
		case 69: goto st124;
		case 101: goto st124;
	}
	goto tr200;
st124:
	if ( ++p == pe )
		goto _test_eof124;
case 124:
	switch( (*p) ) {
		case 84: goto st125;
		case 116: goto st125;
	}
	goto tr200;
st125:
	if ( ++p == pe )
		goto _test_eof125;
case 125:
	if ( (*p) == 10 )
		goto tr211;
	goto tr210;
tr210:
#line 86 "STLAsciiParseRagel.rl"
	{ endFacet(); }
	goto st126;
st126:
	if ( ++p == pe )
		goto _test_eof126;
case 126:
#line 2013 "STLAsciiParseRagel.cc"
	if ( (*p) == 10 )
		goto st127;
	goto st126;
tr211:
#line 86 "STLAsciiParseRagel.rl"
	{ endFacet(); }
	goto st127;
st127:
	if ( ++p == pe )
		goto _test_eof127;
case 127:
#line 2025 "STLAsciiParseRagel.cc"
	switch( (*p) ) {
		case 32: goto tr28;
		case 69: goto tr15;
		case 70: goto tr16;
		case 101: goto tr15;
		case 102: goto tr16;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto tr28;
	goto tr214;
st128:
	if ( ++p == pe )
		goto _test_eof128;
case 128:
	switch( (*p) ) {
		case 43: goto st129;
		case 45: goto st129;
	}
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st130;
	goto tr68;
st129:
	if ( ++p == pe )
		goto _test_eof129;
case 129:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st130;
	goto tr68;
st130:
	if ( ++p == pe )
		goto _test_eof130;
case 130:
	switch( (*p) ) {
		case 9: goto tr179;
		case 10: goto tr180;
		case 32: goto tr179;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st130;
	} else if ( (*p) >= 12 )
		goto tr179;
	goto tr68;
tr175:
#line 82 "STLAsciiParseRagel.rl"
	{ tok = p; /* Local token start */ }
	goto st131;
st131:
	if ( ++p == pe )
		goto _test_eof131;
case 131:
#line 2077 "STLAsciiParseRagel.cc"
	switch( (*p) ) {
		case 9: goto tr179;
		case 10: goto tr180;
		case 32: goto tr179;
		case 46: goto st105;
		case 69: goto st128;
		case 101: goto st128;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st131;
	} else if ( (*p) >= 12 )
		goto tr179;
	goto tr68;
st132:
	if ( ++p == pe )
		goto _test_eof132;
case 132:
	switch( (*p) ) {
		case 43: goto st133;
		case 45: goto st133;
	}
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st134;
	goto tr68;
st133:
	if ( ++p == pe )
		goto _test_eof133;
case 133:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st134;
	goto tr68;
st134:
	if ( ++p == pe )
		goto _test_eof134;
case 134:
	switch( (*p) ) {
		case 9: goto tr170;
		case 32: goto tr170;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st134;
	} else if ( (*p) >= 12 )
		goto tr170;
	goto tr68;
tr166:
#line 82 "STLAsciiParseRagel.rl"
	{ tok = p; /* Local token start */ }
	goto st135;
st135:
	if ( ++p == pe )
		goto _test_eof135;
case 135:
#line 2132 "STLAsciiParseRagel.cc"
	switch( (*p) ) {
		case 9: goto tr170;
		case 32: goto tr170;
		case 46: goto st101;
		case 69: goto st132;
		case 101: goto st132;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st135;
	} else if ( (*p) >= 12 )
		goto tr170;
	goto tr68;
st136:
	if ( ++p == pe )
		goto _test_eof136;
case 136:
	switch( (*p) ) {
		case 43: goto st137;
		case 45: goto st137;
	}
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st138;
	goto tr68;
st137:
	if ( ++p == pe )
		goto _test_eof137;
case 137:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st138;
	goto tr68;
st138:
	if ( ++p == pe )
		goto _test_eof138;
case 138:
	switch( (*p) ) {
		case 9: goto tr161;
		case 32: goto tr161;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st138;
	} else if ( (*p) >= 12 )
		goto tr161;
	goto tr68;
tr157:
#line 82 "STLAsciiParseRagel.rl"
	{ tok = p; /* Local token start */ }
	goto st139;
st139:
	if ( ++p == pe )
		goto _test_eof139;
case 139:
#line 2186 "STLAsciiParseRagel.cc"
	switch( (*p) ) {
		case 9: goto tr161;
		case 32: goto tr161;
		case 46: goto st97;
		case 69: goto st136;
		case 101: goto st136;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st139;
	} else if ( (*p) >= 12 )
		goto tr161;
	goto tr68;
st140:
	if ( ++p == pe )
		goto _test_eof140;
case 140:
	switch( (*p) ) {
		case 43: goto st141;
		case 45: goto st141;
	}
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st142;
	goto tr68;
st141:
	if ( ++p == pe )
		goto _test_eof141;
case 141:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st142;
	goto tr68;
st142:
	if ( ++p == pe )
		goto _test_eof142;
case 142:
	switch( (*p) ) {
		case 9: goto tr140;
		case 10: goto tr141;
		case 32: goto tr140;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st142;
	} else if ( (*p) >= 12 )
		goto tr140;
	goto tr68;
tr136:
#line 82 "STLAsciiParseRagel.rl"
	{ tok = p; /* Local token start */ }
	goto st143;
st143:
	if ( ++p == pe )
		goto _test_eof143;
case 143:
#line 2241 "STLAsciiParseRagel.cc"
	switch( (*p) ) {
		case 9: goto tr140;
		case 10: goto tr141;
		case 32: goto tr140;
		case 46: goto st84;
		case 69: goto st140;
		case 101: goto st140;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st143;
	} else if ( (*p) >= 12 )
		goto tr140;
	goto tr68;
st144:
	if ( ++p == pe )
		goto _test_eof144;
case 144:
	switch( (*p) ) {
		case 43: goto st145;
		case 45: goto st145;
	}
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st146;
	goto tr68;
st145:
	if ( ++p == pe )
		goto _test_eof145;
case 145:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st146;
	goto tr68;
st146:
	if ( ++p == pe )
		goto _test_eof146;
case 146:
	switch( (*p) ) {
		case 9: goto tr131;
		case 32: goto tr131;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st146;
	} else if ( (*p) >= 12 )
		goto tr131;
	goto tr68;
tr127:
#line 82 "STLAsciiParseRagel.rl"
	{ tok = p; /* Local token start */ }
	goto st147;
st147:
	if ( ++p == pe )
		goto _test_eof147;
case 147:
#line 2296 "STLAsciiParseRagel.cc"
	switch( (*p) ) {
		case 9: goto tr131;
		case 32: goto tr131;
		case 46: goto st80;
		case 69: goto st144;
		case 101: goto st144;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st147;
	} else if ( (*p) >= 12 )
		goto tr131;
	goto tr68;
st148:
	if ( ++p == pe )
		goto _test_eof148;
case 148:
	switch( (*p) ) {
		case 43: goto st149;
		case 45: goto st149;
	}
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st150;
	goto tr68;
st149:
	if ( ++p == pe )
		goto _test_eof149;
case 149:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st150;
	goto tr68;
st150:
	if ( ++p == pe )
		goto _test_eof150;
case 150:
	switch( (*p) ) {
		case 9: goto tr122;
		case 32: goto tr122;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st150;
	} else if ( (*p) >= 12 )
		goto tr122;
	goto tr68;
tr118:
#line 82 "STLAsciiParseRagel.rl"
	{ tok = p; /* Local token start */ }
	goto st151;
st151:
	if ( ++p == pe )
		goto _test_eof151;
case 151:
#line 2350 "STLAsciiParseRagel.cc"
	switch( (*p) ) {
		case 9: goto tr122;
		case 32: goto tr122;
		case 46: goto st76;
		case 69: goto st148;
		case 101: goto st148;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st151;
	} else if ( (*p) >= 12 )
		goto tr122;
	goto tr68;
st152:
	if ( ++p == pe )
		goto _test_eof152;
case 152:
	switch( (*p) ) {
		case 43: goto st153;
		case 45: goto st153;
	}
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st154;
	goto tr68;
st153:
	if ( ++p == pe )
		goto _test_eof153;
case 153:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st154;
	goto tr68;
st154:
	if ( ++p == pe )
		goto _test_eof154;
case 154:
	switch( (*p) ) {
		case 9: goto tr101;
		case 10: goto tr102;
		case 32: goto tr101;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st154;
	} else if ( (*p) >= 12 )
		goto tr101;
	goto tr68;
tr97:
#line 82 "STLAsciiParseRagel.rl"
	{ tok = p; /* Local token start */ }
	goto st155;
st155:
	if ( ++p == pe )
		goto _test_eof155;
case 155:
#line 2405 "STLAsciiParseRagel.cc"
	switch( (*p) ) {
		case 9: goto tr101;
		case 10: goto tr102;
		case 32: goto tr101;
		case 46: goto st63;
		case 69: goto st152;
		case 101: goto st152;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st155;
	} else if ( (*p) >= 12 )
		goto tr101;
	goto tr68;
st156:
	if ( ++p == pe )
		goto _test_eof156;
case 156:
	switch( (*p) ) {
		case 43: goto st157;
		case 45: goto st157;
	}
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st158;
	goto tr68;
st157:
	if ( ++p == pe )
		goto _test_eof157;
case 157:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st158;
	goto tr68;
st158:
	if ( ++p == pe )
		goto _test_eof158;
case 158:
	switch( (*p) ) {
		case 9: goto tr92;
		case 32: goto tr92;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st158;
	} else if ( (*p) >= 12 )
		goto tr92;
	goto tr68;
tr88:
#line 82 "STLAsciiParseRagel.rl"
	{ tok = p; /* Local token start */ }
	goto st159;
st159:
	if ( ++p == pe )
		goto _test_eof159;
case 159:
#line 2460 "STLAsciiParseRagel.cc"
	switch( (*p) ) {
		case 9: goto tr92;
		case 32: goto tr92;
		case 46: goto st59;
		case 69: goto st156;
		case 101: goto st156;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st159;
	} else if ( (*p) >= 12 )
		goto tr92;
	goto tr68;
st160:
	if ( ++p == pe )
		goto _test_eof160;
case 160:
	switch( (*p) ) {
		case 43: goto st161;
		case 45: goto st161;
	}
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st162;
	goto tr68;
st161:
	if ( ++p == pe )
		goto _test_eof161;
case 161:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st162;
	goto tr68;
st162:
	if ( ++p == pe )
		goto _test_eof162;
case 162:
	switch( (*p) ) {
		case 9: goto tr83;
		case 32: goto tr83;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st162;
	} else if ( (*p) >= 12 )
		goto tr83;
	goto tr68;
tr79:
#line 82 "STLAsciiParseRagel.rl"
	{ tok = p; /* Local token start */ }
	goto st163;
st163:
	if ( ++p == pe )
		goto _test_eof163;
case 163:
#line 2514 "STLAsciiParseRagel.cc"
	switch( (*p) ) {
		case 9: goto tr83;
		case 32: goto tr83;
		case 46: goto st55;
		case 69: goto st160;
		case 101: goto st160;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st163;
	} else if ( (*p) >= 12 )
		goto tr83;
	goto tr68;
	}
	_test_eof1: cs = 1; goto _test_eof; 
	_test_eof2: cs = 2; goto _test_eof; 
	_test_eof3: cs = 3; goto _test_eof; 
	_test_eof4: cs = 4; goto _test_eof; 
	_test_eof5: cs = 5; goto _test_eof; 
	_test_eof6: cs = 6; goto _test_eof; 
	_test_eof7: cs = 7; goto _test_eof; 
	_test_eof8: cs = 8; goto _test_eof; 
	_test_eof9: cs = 9; goto _test_eof; 
	_test_eof10: cs = 10; goto _test_eof; 
	_test_eof11: cs = 11; goto _test_eof; 
	_test_eof12: cs = 12; goto _test_eof; 
	_test_eof13: cs = 13; goto _test_eof; 
	_test_eof14: cs = 14; goto _test_eof; 
	_test_eof15: cs = 15; goto _test_eof; 
	_test_eof16: cs = 16; goto _test_eof; 
	_test_eof17: cs = 17; goto _test_eof; 
	_test_eof18: cs = 18; goto _test_eof; 
	_test_eof19: cs = 19; goto _test_eof; 
	_test_eof20: cs = 20; goto _test_eof; 
	_test_eof21: cs = 21; goto _test_eof; 
	_test_eof22: cs = 22; goto _test_eof; 
	_test_eof23: cs = 23; goto _test_eof; 
	_test_eof24: cs = 24; goto _test_eof; 
	_test_eof164: cs = 164; goto _test_eof; 
	_test_eof165: cs = 165; goto _test_eof; 
	_test_eof25: cs = 25; goto _test_eof; 
	_test_eof26: cs = 26; goto _test_eof; 
	_test_eof27: cs = 27; goto _test_eof; 
	_test_eof28: cs = 28; goto _test_eof; 
	_test_eof29: cs = 29; goto _test_eof; 
	_test_eof30: cs = 30; goto _test_eof; 
	_test_eof31: cs = 31; goto _test_eof; 
	_test_eof32: cs = 32; goto _test_eof; 
	_test_eof33: cs = 33; goto _test_eof; 
	_test_eof34: cs = 34; goto _test_eof; 
	_test_eof35: cs = 35; goto _test_eof; 
	_test_eof36: cs = 36; goto _test_eof; 
	_test_eof37: cs = 37; goto _test_eof; 
	_test_eof38: cs = 38; goto _test_eof; 
	_test_eof39: cs = 39; goto _test_eof; 
	_test_eof40: cs = 40; goto _test_eof; 
	_test_eof41: cs = 41; goto _test_eof; 
	_test_eof42: cs = 42; goto _test_eof; 
	_test_eof43: cs = 43; goto _test_eof; 
	_test_eof44: cs = 44; goto _test_eof; 
	_test_eof45: cs = 45; goto _test_eof; 
	_test_eof46: cs = 46; goto _test_eof; 
	_test_eof47: cs = 47; goto _test_eof; 
	_test_eof48: cs = 48; goto _test_eof; 
	_test_eof49: cs = 49; goto _test_eof; 
	_test_eof50: cs = 50; goto _test_eof; 
	_test_eof51: cs = 51; goto _test_eof; 
	_test_eof52: cs = 52; goto _test_eof; 
	_test_eof53: cs = 53; goto _test_eof; 
	_test_eof54: cs = 54; goto _test_eof; 
	_test_eof55: cs = 55; goto _test_eof; 
	_test_eof56: cs = 56; goto _test_eof; 
	_test_eof57: cs = 57; goto _test_eof; 
	_test_eof58: cs = 58; goto _test_eof; 
	_test_eof59: cs = 59; goto _test_eof; 
	_test_eof60: cs = 60; goto _test_eof; 
	_test_eof61: cs = 61; goto _test_eof; 
	_test_eof62: cs = 62; goto _test_eof; 
	_test_eof63: cs = 63; goto _test_eof; 
	_test_eof64: cs = 64; goto _test_eof; 
	_test_eof65: cs = 65; goto _test_eof; 
	_test_eof66: cs = 66; goto _test_eof; 
	_test_eof67: cs = 67; goto _test_eof; 
	_test_eof68: cs = 68; goto _test_eof; 
	_test_eof69: cs = 69; goto _test_eof; 
	_test_eof70: cs = 70; goto _test_eof; 
	_test_eof71: cs = 71; goto _test_eof; 
	_test_eof72: cs = 72; goto _test_eof; 
	_test_eof73: cs = 73; goto _test_eof; 
	_test_eof74: cs = 74; goto _test_eof; 
	_test_eof75: cs = 75; goto _test_eof; 
	_test_eof76: cs = 76; goto _test_eof; 
	_test_eof77: cs = 77; goto _test_eof; 
	_test_eof78: cs = 78; goto _test_eof; 
	_test_eof79: cs = 79; goto _test_eof; 
	_test_eof80: cs = 80; goto _test_eof; 
	_test_eof81: cs = 81; goto _test_eof; 
	_test_eof82: cs = 82; goto _test_eof; 
	_test_eof83: cs = 83; goto _test_eof; 
	_test_eof84: cs = 84; goto _test_eof; 
	_test_eof85: cs = 85; goto _test_eof; 
	_test_eof86: cs = 86; goto _test_eof; 
	_test_eof87: cs = 87; goto _test_eof; 
	_test_eof88: cs = 88; goto _test_eof; 
	_test_eof89: cs = 89; goto _test_eof; 
	_test_eof90: cs = 90; goto _test_eof; 
	_test_eof91: cs = 91; goto _test_eof; 
	_test_eof92: cs = 92; goto _test_eof; 
	_test_eof93: cs = 93; goto _test_eof; 
	_test_eof94: cs = 94; goto _test_eof; 
	_test_eof95: cs = 95; goto _test_eof; 
	_test_eof96: cs = 96; goto _test_eof; 
	_test_eof97: cs = 97; goto _test_eof; 
	_test_eof98: cs = 98; goto _test_eof; 
	_test_eof99: cs = 99; goto _test_eof; 
	_test_eof100: cs = 100; goto _test_eof; 
	_test_eof101: cs = 101; goto _test_eof; 
	_test_eof102: cs = 102; goto _test_eof; 
	_test_eof103: cs = 103; goto _test_eof; 
	_test_eof104: cs = 104; goto _test_eof; 
	_test_eof105: cs = 105; goto _test_eof; 
	_test_eof106: cs = 106; goto _test_eof; 
	_test_eof107: cs = 107; goto _test_eof; 
	_test_eof108: cs = 108; goto _test_eof; 
	_test_eof109: cs = 109; goto _test_eof; 
	_test_eof110: cs = 110; goto _test_eof; 
	_test_eof111: cs = 111; goto _test_eof; 
	_test_eof112: cs = 112; goto _test_eof; 
	_test_eof113: cs = 113; goto _test_eof; 
	_test_eof114: cs = 114; goto _test_eof; 
	_test_eof115: cs = 115; goto _test_eof; 
	_test_eof116: cs = 116; goto _test_eof; 
	_test_eof117: cs = 117; goto _test_eof; 
	_test_eof118: cs = 118; goto _test_eof; 
	_test_eof119: cs = 119; goto _test_eof; 
	_test_eof120: cs = 120; goto _test_eof; 
	_test_eof121: cs = 121; goto _test_eof; 
	_test_eof122: cs = 122; goto _test_eof; 
	_test_eof123: cs = 123; goto _test_eof; 
	_test_eof124: cs = 124; goto _test_eof; 
	_test_eof125: cs = 125; goto _test_eof; 
	_test_eof126: cs = 126; goto _test_eof; 
	_test_eof127: cs = 127; goto _test_eof; 
	_test_eof128: cs = 128; goto _test_eof; 
	_test_eof129: cs = 129; goto _test_eof; 
	_test_eof130: cs = 130; goto _test_eof; 
	_test_eof131: cs = 131; goto _test_eof; 
	_test_eof132: cs = 132; goto _test_eof; 
	_test_eof133: cs = 133; goto _test_eof; 
	_test_eof134: cs = 134; goto _test_eof; 
	_test_eof135: cs = 135; goto _test_eof; 
	_test_eof136: cs = 136; goto _test_eof; 
	_test_eof137: cs = 137; goto _test_eof; 
	_test_eof138: cs = 138; goto _test_eof; 
	_test_eof139: cs = 139; goto _test_eof; 
	_test_eof140: cs = 140; goto _test_eof; 
	_test_eof141: cs = 141; goto _test_eof; 
	_test_eof142: cs = 142; goto _test_eof; 
	_test_eof143: cs = 143; goto _test_eof; 
	_test_eof144: cs = 144; goto _test_eof; 
	_test_eof145: cs = 145; goto _test_eof; 
	_test_eof146: cs = 146; goto _test_eof; 
	_test_eof147: cs = 147; goto _test_eof; 
	_test_eof148: cs = 148; goto _test_eof; 
	_test_eof149: cs = 149; goto _test_eof; 
	_test_eof150: cs = 150; goto _test_eof; 
	_test_eof151: cs = 151; goto _test_eof; 
	_test_eof152: cs = 152; goto _test_eof; 
	_test_eof153: cs = 153; goto _test_eof; 
	_test_eof154: cs = 154; goto _test_eof; 
	_test_eof155: cs = 155; goto _test_eof; 
	_test_eof156: cs = 156; goto _test_eof; 
	_test_eof157: cs = 157; goto _test_eof; 
	_test_eof158: cs = 158; goto _test_eof; 
	_test_eof159: cs = 159; goto _test_eof; 
	_test_eof160: cs = 160; goto _test_eof; 
	_test_eof161: cs = 161; goto _test_eof; 
	_test_eof162: cs = 162; goto _test_eof; 
	_test_eof163: cs = 163; goto _test_eof; 

	_test_eof: {}
	if ( p == eof )
	{
	switch ( cs ) {
	case 164: 
#line 81 "STLAsciiParseRagel.rl"
	{ ++lineNum_; /* Count '\n' occurrences */ }
	break;
	case 1: 
	case 2: 
	case 3: 
	case 4: 
	case 5: 
	case 6: 
	case 7: 
#line 99 "STLAsciiParseRagel.rl"
	{ die("solid", p, pe); }
	break;
	case 17: 
	case 18: 
	case 19: 
	case 20: 
	case 21: 
	case 22: 
	case 23: 
	case 24: 
#line 100 "STLAsciiParseRagel.rl"
	{ die("endsolid", p, pe); }
	break;
	case 25: 
	case 26: 
	case 27: 
	case 28: 
	case 29: 
	case 30: 
	case 31: 
#line 101 "STLAsciiParseRagel.rl"
	{ die("facet", p, pe); }
	break;
	case 117: 
	case 118: 
	case 119: 
	case 120: 
	case 121: 
	case 122: 
	case 123: 
	case 124: 
	case 125: 
	case 126: 
#line 102 "STLAsciiParseRagel.rl"
	{ die("endfacet", p, pe); }
	break;
	case 33: 
	case 34: 
	case 35: 
	case 36: 
	case 37: 
	case 38: 
	case 39: 
	case 40: 
	case 41: 
	case 42: 
	case 43: 
#line 103 "STLAsciiParseRagel.rl"
	{ die("loop", p, pe); }
	break;
	case 108: 
	case 109: 
	case 110: 
	case 111: 
	case 112: 
	case 113: 
	case 114: 
	case 115: 
#line 104 "STLAsciiParseRagel.rl"
	{ die("endloop", p, pe); }
	break;
	case 45: 
	case 46: 
	case 47: 
	case 48: 
	case 49: 
	case 50: 
	case 51: 
	case 52: 
	case 53: 
	case 54: 
	case 55: 
	case 56: 
	case 57: 
	case 58: 
	case 59: 
	case 60: 
	case 61: 
	case 62: 
	case 63: 
	case 64: 
	case 65: 
	case 66: 
	case 67: 
	case 68: 
	case 69: 
	case 70: 
	case 71: 
	case 72: 
	case 73: 
	case 74: 
	case 75: 
	case 76: 
	case 77: 
	case 78: 
	case 79: 
	case 80: 
	case 81: 
	case 82: 
	case 83: 
	case 84: 
	case 85: 
	case 86: 
	case 87: 
	case 88: 
	case 89: 
	case 90: 
	case 91: 
	case 92: 
	case 93: 
	case 94: 
	case 95: 
	case 96: 
	case 97: 
	case 98: 
	case 99: 
	case 100: 
	case 101: 
	case 102: 
	case 103: 
	case 104: 
	case 105: 
	case 106: 
	case 128: 
	case 129: 
	case 130: 
	case 131: 
	case 132: 
	case 133: 
	case 134: 
	case 135: 
	case 136: 
	case 137: 
	case 138: 
	case 139: 
	case 140: 
	case 141: 
	case 142: 
	case 143: 
	case 144: 
	case 145: 
	case 146: 
	case 147: 
	case 148: 
	case 149: 
	case 150: 
	case 151: 
	case 152: 
	case 153: 
	case 154: 
	case 155: 
	case 156: 
	case 157: 
	case 158: 
	case 159: 
	case 160: 
	case 161: 
	case 162: 
	case 163: 
#line 105 "STLAsciiParseRagel.rl"
	{ die("vertex", p, pe); }
	break;
	case 9: 
	case 15: 
	case 16: 
#line 101 "STLAsciiParseRagel.rl"
	{ die("facet", p, pe); }
#line 100 "STLAsciiParseRagel.rl"
	{ die("endsolid", p, pe); }
	break;
	case 32: 
#line 101 "STLAsciiParseRagel.rl"
	{ die("facet", p, pe); }
#line 103 "STLAsciiParseRagel.rl"
	{ die("loop", p, pe); }
	break;
	case 44: 
#line 103 "STLAsciiParseRagel.rl"
	{ die("loop", p, pe); }
#line 105 "STLAsciiParseRagel.rl"
	{ die("vertex", p, pe); }
	break;
	case 116: 
#line 104 "STLAsciiParseRagel.rl"
	{ die("endloop", p, pe); }
#line 102 "STLAsciiParseRagel.rl"
	{ die("endfacet", p, pe); }
	break;
	case 107: 
#line 105 "STLAsciiParseRagel.rl"
	{ die("vertex", p, pe); }
#line 104 "STLAsciiParseRagel.rl"
	{ die("endloop", p, pe); }
	break;
	case 8: 
#line 99 "STLAsciiParseRagel.rl"
	{ die("solid", p, pe); }
#line 101 "STLAsciiParseRagel.rl"
	{ die("facet", p, pe); }
#line 100 "STLAsciiParseRagel.rl"
	{ die("endsolid", p, pe); }
	break;
	case 127: 
#line 101 "STLAsciiParseRagel.rl"
	{ die("facet", p, pe); }
#line 102 "STLAsciiParseRagel.rl"
	{ die("endfacet", p, pe); }
#line 100 "STLAsciiParseRagel.rl"
	{ die("endsolid", p, pe); }
	break;
#line 2921 "STLAsciiParseRagel.cc"
	}
	}

	_out: {}
	}

#line 261 "STLAsciiParseRagel.rl"
       /* ^^^ FSM execution here ^^^ */;

        if (0 == cs)
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
