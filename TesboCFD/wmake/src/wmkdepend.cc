
#line 1 "wmkdepend.rl"
/*---------------------------------*- C -*-----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2022 OpenCFD Ltd.
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
    wmkdepend

Description
    A faster dependency list generator that emulates the behaviour and
    output of 'cpp -M'.

    For each (quoted) include directive detected, the specified include
    directories are searched and that file (if found) is also scanned.

    Each include file is visited only once, which helps make this faster
    than cpp. It also handles missing files somewhat more gracefully.

    The scanner is built with Ragel.
    The goto-based finite state machine (FSM) is generated with

        ragel -G2 -o wmkdepend.cc wmkdepend.rl

    The FSM can be visualized (eg, in a browser) with the following command

        ragel -pV wmkdepend.rl | dot -Tsvg -owmkdepend.svg

Usage
    wmkdepend [-Idir..] [-iheader...] [-eENV...] [-oFile] filename

Note
    May not capture all possible corner cases or line continuations such as

        #include \
            "file.H"

\*---------------------------------------------------------------------------*/
/*
 * With cpp:
 *
 * cpp -x c++ -std=c++11
 *     -M -DWM_$(WM_PRECISION_OPTION) -DWM_LABEL_SIZE=$(WM_LABEL_SIZE) |
 * sed -e 's,^$(WM_PROJECT_DIR)/,$$(WM_PROJECT_DIR)/,' \
 *     -e 's,^$(WM_THIRD_PARTY_DIR)/,$$(WM_THIRD_PARTY_DIR)/,'
*/

#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <forward_list>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#pragma GCC diagnostic ignored "-Wimplicit-fallthrough"
#pragma GCC diagnostic ignored "-Wunused-const-variable"

// Length of the input read buffer
#define INBUFLEN 16384

// The executable name (for messages), without requiring access to argv[]
#define EXENAME  "wmkdepend"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

//- Program usage
void usage()
{
    std::cerr
        <<
        "\nUsage: " EXENAME
        " [-Idir...] [-iheader...] [-eENV...] [-oFile]"
        " filename\n\n"
        "  -Idir     Directories to be searched for headers.\n"
        "  -iheader  Headers to be ignored.\n"
        "  -eENV     Environment variable path substitutions.\n"
        "  -oFile    Write output to File.\n"
        "  -q        Suppress 'No such file' warnings.\n"
        "  -v        Report each include file to stderr.\n"
        "\nDependency list generator, similar to 'cpp -M'\n\n";
}

//- Suppress some error messages
bool optQuiet = false;

//- Verbose progress
bool optVerbose = false;

//- The top-level source file being processed
std::string sourceFile;

//- The output stream
FILE* output = stdout;

//- All file opening and writing
namespace Files
{
    //- Include directories
    static std::vector<std::string> dirs;

    //- Set of files already visited
    static std::unordered_set<std::string> visited;

    //- Environment substitution
    struct envEntry
    {
        std::string name;
        std::string value;
        size_t len;

        envEntry(std::string&& k, std::string&& v)
        :
            name(std::move(k)),
            value(std::move(v)),
            len(value.size())
        {}
    };

    //- List of environ variables to substitute
    std::forward_list<envEntry> envlist;

    //- Clear, reset all variables
    void reset(unsigned ndirs = 0u)
    {
        dirs.clear();
        visited.clear();
        envlist.clear();
        if (ndirs)
        {
            dirs.reserve(ndirs);
        }
    }


    //- Add environ replacement
    //
    // Eg for 'WM_PROJECT_DIR' stores the pair
    //     '$(WM_PROJECT_DIR)/'
    //     '/path/openfoam/project/'
    //
    // No-op if the environment does not exists or is empty
    void addEnv(const char* key)
    {
        const size_t keylen = key ? strlen(key) : 0;
        if (!keylen) return;

        const char *val = ::getenv(key);

        if (!val || !*val) return;

        // The "$(ENV)/" portion
        std::string oldText;
        oldText.reserve(keylen+4);
        oldText.append("$(").append(key,keylen).append(")/");

        // The "/env/value/" portion
        std::string newText(val);
        if (newText.back() != '/')
        {
            newText.append(1, '/');
        }

        envlist.emplace_front(std::move(oldText), std::move(newText));
    }


    //- Open a file for reading and emit its qualified name to output.
    //
    //  Uses env substitutions at the beginning of the path
    //
    //  Eg,
    //      /path/openfoam/project/directory/name
    //   -> $(WM_PROJECT_DIR)/directory/name
    //
    //  \return nullptr on failure
    FILE* openAndEmit(const std::string& fileName)
    {
        const auto len = fileName.size();
        const char *fname = fileName.c_str();

        FILE *infile = ::fopen(fname, "r");

        if (infile)
        {
            // Mark as having been visited
            visited.insert(fileName);

            for (const auto& entry : envlist)
            {
                if
                (
                    len > entry.len
                 && !fileName.compare(0, entry.len, entry.value)
                )
                {
                    fname += entry.len;  // Now positioned after the '/'
                    fputs(entry.name.c_str(), output);
                    break;
                }
            }

            fputs(fname, output);
            fputs(" \\\n", output);
        }
        else if (errno == EMFILE)
        {
            std::cerr
                << EXENAME ": too many open files while opening '"
                << fileName << "'\n"
                << "Please change your open descriptor limit\n";
        }

        return infile;
    }


    //- Open a not previously visited file for reading,
    //  using the include dirs as required.
    //
    //  On success, emits the resolved name on stdout
    //
    //  \return nullptr on failure
    FILE* open(const std::string& fileName)
    {
        // Bad file name, or already visited
        if (fileName.empty() || visited.find(fileName) != visited.end())
        {
            return nullptr;
        }

        FILE* infile = openAndEmit(fileName);
        if (!infile)
        {
            std::string fullName;

            for (const auto& d : dirs)
            {
                if (d.empty())
                {
                    continue;
                }

                fullName.clear();
                fullName.reserve(d.size() + fileName.size() + 1);

                fullName.append(d);
                if (fullName.back() != '/')
                {
                    fullName.append("/");
                }
                fullName.append(fileName);

                infile = openAndEmit(fullName);

                if (infile)
                {
                    break;
                }
            }
        }

        // Mark as having been visited.
        // This also avoids re-testing and multiple messages.
        visited.insert(fileName);

        // Report failures
        if (!infile && !optQuiet)
        {
            std::cerr
                << EXENAME ": could not open '"
                << fileName << "' for source file '"
                << sourceFile << "'";

            if (dirs.size())
            {
                std::cerr << ": " << strerror(errno);
            }

            std::cerr << '\n' << std::flush;
        }

        return infile;
    }

} // End namespace Files


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Ragel machine definition
// Ragel variables (p, pe, eof, cs, top, stack, ts, te, act) defined later...
//
// Can use 'variable p xxx;' etc to change these names


#line 320 "wmkdepend.cc"
static const int wmkdep_start = 21;
static const int wmkdep_error = 0;

static const int wmkdep_en_main = 21;


#line 346 "wmkdepend.rl"



/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void processFile(std::string fileName);

//
// Open a file and process.
// The file name is given by the [first,last) range
//
void processFile(const char* first, const char* last)
{
    // Extra safety
    if (first && last && last > first)
    {
        processFile(std::string(first, last));
    }
}


//
// Open a file and process
//
void processFile(std::string fileName)
{
    FILE* infile = Files::open(fileName);
    if (optVerbose)
    {
        std::cerr << fileName << '\n';
    }
    if (!infile) return;

    // ts, te  = Ragel token start/end points (default naming)
    // act, cs = Ragel action, code state, respectively (default naming)
    char *ts, *te;
    int act, cs;

    // Initialize FSM variables
    
#line 368 "wmkdepend.cc"
	{
	cs = wmkdep_start;
	ts = 0;
	te = 0;
	act = 0;
	}

#line 385 "wmkdepend.rl"
   /* ^^^ FSM initialization here ^^^ */;

    // Local token start
    char *tok = nullptr;

    // Buffering
    char inbuf[INBUFLEN];
    size_t pending = 0;

    // Processing loop (as per Ragel pdf example)
    for (bool good = true; good; /*nil*/)
    {
        if (pending >= INBUFLEN)
        {
            // We overfilled the buffer while trying to scan a token...
            std::cerr
                << EXENAME ": buffer full while scanning '"
                << fileName << "'\n";
            break;
        }

        char *data = inbuf + pending;   // current data buffer
        const size_t buflen = INBUFLEN - pending; // space in buffer

        // p,pe = Ragel parsing point and parsing end (default naming)
        // eof  = Ragel EOF point (default naming)

        const size_t gcount = ::fread(data, 1, buflen, infile);

        char *p = data;
        char *pe = p + gcount;
        char *eof = nullptr;

        if (::feof(infile))
        {
            eof = pe;   // Tag 'pe' as being the EOF for the FSM as well
            good = false;
        }
        else if (!gcount)
        {
            break;
        }

        
#line 421 "wmkdepend.cc"
	{
	if ( p == pe )
		goto _test_eof;
	switch ( cs )
	{
tr0:
#line 1 "NONE"
	{	switch( act ) {
	case 0:
	{{goto st0;}}
	break;
	default:
	{{p = ((te))-1;}}
	break;
	}
	}
	goto st21;
tr2:
#line 344 "wmkdepend.rl"
	{te = p+1;}
	goto st21;
tr17:
#line 344 "wmkdepend.rl"
	{{p = ((te))-1;}}
	goto st21;
tr21:
#line 339 "wmkdepend.rl"
	{te = p+1;}
	goto st21;
tr29:
#line 342 "wmkdepend.rl"
	{te = p+1;}
	goto st21;
tr31:
#line 341 "wmkdepend.rl"
	{te = p+1;}
	goto st21;
tr36:
#line 336 "wmkdepend.rl"
	{te = p;p--;}
	goto st21;
tr37:
#line 344 "wmkdepend.rl"
	{te = p;p--;}
	goto st21;
tr38:
#line 342 "wmkdepend.rl"
	{te = p;p--;}
	goto st21;
st21:
#line 1 "NONE"
	{ts = 0;}
#line 1 "NONE"
	{act = 0;}
	if ( ++p == pe )
		goto _test_eof21;
case 21:
#line 1 "NONE"
	{ts = p;}
#line 481 "wmkdepend.cc"
	switch( (*p) ) {
		case 10: goto st23;
		case 11: goto tr34;
		case 32: goto tr32;
		case 35: goto st2;
		case 47: goto st15;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto tr32;
	goto st1;
st1:
	if ( ++p == pe )
		goto _test_eof1;
case 1:
	if ( (*p) == 10 )
		goto tr2;
	goto st1;
tr32:
#line 1 "NONE"
	{te = p+1;}
#line 336 "wmkdepend.rl"
	{act = 1;}
	goto st22;
st22:
	if ( ++p == pe )
		goto _test_eof22;
case 22:
#line 509 "wmkdepend.cc"
	switch( (*p) ) {
		case 10: goto st23;
		case 11: goto tr34;
		case 32: goto tr32;
		case 35: goto st2;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto tr32;
	goto st1;
st23:
	if ( ++p == pe )
		goto _test_eof23;
case 23:
	if ( (*p) == 32 )
		goto st23;
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st23;
	goto tr36;
tr34:
#line 1 "NONE"
	{te = p+1;}
#line 336 "wmkdepend.rl"
	{act = 1;}
	goto st24;
st24:
	if ( ++p == pe )
		goto _test_eof24;
case 24:
#line 538 "wmkdepend.cc"
	switch( (*p) ) {
		case 10: goto st23;
		case 32: goto tr34;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto tr34;
	goto st1;
st2:
	if ( ++p == pe )
		goto _test_eof2;
case 2:
	switch( (*p) ) {
		case 9: goto st2;
		case 10: goto tr2;
		case 32: goto st2;
		case 105: goto st3;
	}
	if ( 12 <= (*p) && (*p) <= 13 )
		goto st2;
	goto st1;
st3:
	if ( ++p == pe )
		goto _test_eof3;
case 3:
	switch( (*p) ) {
		case 10: goto tr2;
		case 110: goto st4;
	}
	goto st1;
st4:
	if ( ++p == pe )
		goto _test_eof4;
case 4:
	switch( (*p) ) {
		case 10: goto tr2;
		case 99: goto st5;
	}
	goto st1;
st5:
	if ( ++p == pe )
		goto _test_eof5;
case 5:
	switch( (*p) ) {
		case 10: goto tr2;
		case 108: goto st6;
	}
	goto st1;
st6:
	if ( ++p == pe )
		goto _test_eof6;
case 6:
	switch( (*p) ) {
		case 10: goto tr2;
		case 117: goto st7;
	}
	goto st1;
st7:
	if ( ++p == pe )
		goto _test_eof7;
case 7:
	switch( (*p) ) {
		case 10: goto tr2;
		case 100: goto st8;
	}
	goto st1;
st8:
	if ( ++p == pe )
		goto _test_eof8;
case 8:
	switch( (*p) ) {
		case 10: goto tr2;
		case 101: goto st9;
	}
	goto st1;
st9:
	if ( ++p == pe )
		goto _test_eof9;
case 9:
	switch( (*p) ) {
		case 9: goto st9;
		case 10: goto tr2;
		case 32: goto st9;
		case 34: goto st10;
	}
	if ( 12 <= (*p) && (*p) <= 13 )
		goto st9;
	goto st1;
st10:
	if ( ++p == pe )
		goto _test_eof10;
case 10:
	switch( (*p) ) {
		case 10: goto tr13;
		case 34: goto st1;
	}
	goto tr12;
tr12:
#line 320 "wmkdepend.rl"
	{ tok = p;  /* Local token start */ }
	goto st11;
st11:
	if ( ++p == pe )
		goto _test_eof11;
case 11:
#line 643 "wmkdepend.cc"
	switch( (*p) ) {
		case 10: goto tr15;
		case 34: goto tr16;
	}
	goto st11;
tr13:
#line 1 "NONE"
	{te = p+1;}
#line 320 "wmkdepend.rl"
	{ tok = p;  /* Local token start */ }
	goto st25;
tr15:
#line 1 "NONE"
	{te = p+1;}
	goto st25;
st25:
	if ( ++p == pe )
		goto _test_eof25;
case 25:
#line 663 "wmkdepend.cc"
	if ( (*p) == 34 )
		goto tr19;
	goto st12;
st12:
	if ( ++p == pe )
		goto _test_eof12;
case 12:
	if ( (*p) == 34 )
		goto tr19;
	goto st12;
tr19:
#line 322 "wmkdepend.rl"
	{
        processFile(tok, p);
        tok = nullptr;          /* Done with buffer */
    }
	goto st13;
st13:
	if ( ++p == pe )
		goto _test_eof13;
case 13:
#line 685 "wmkdepend.cc"
	if ( (*p) == 10 )
		goto tr21;
	goto st13;
tr16:
#line 322 "wmkdepend.rl"
	{
        processFile(tok, p);
        tok = nullptr;          /* Done with buffer */
    }
	goto st14;
st14:
	if ( ++p == pe )
		goto _test_eof14;
case 14:
#line 700 "wmkdepend.cc"
	if ( (*p) == 10 )
		goto tr21;
	goto st14;
st15:
	if ( ++p == pe )
		goto _test_eof15;
case 15:
	switch( (*p) ) {
		case 10: goto tr2;
		case 42: goto st16;
		case 47: goto st20;
	}
	goto st1;
st16:
	if ( ++p == pe )
		goto _test_eof16;
case 16:
	switch( (*p) ) {
		case 10: goto tr25;
		case 42: goto st19;
	}
	goto st16;
tr25:
#line 1 "NONE"
	{te = p+1;}
	goto st26;
st26:
	if ( ++p == pe )
		goto _test_eof26;
case 26:
#line 731 "wmkdepend.cc"
	if ( (*p) == 42 )
		goto st18;
	goto st17;
st17:
	if ( ++p == pe )
		goto _test_eof17;
case 17:
	if ( (*p) == 42 )
		goto st18;
	goto st17;
st18:
	if ( ++p == pe )
		goto _test_eof18;
case 18:
	switch( (*p) ) {
		case 42: goto st18;
		case 47: goto tr29;
	}
	goto st17;
st19:
	if ( ++p == pe )
		goto _test_eof19;
case 19:
	switch( (*p) ) {
		case 10: goto tr25;
		case 42: goto st19;
		case 47: goto tr30;
	}
	goto st16;
tr30:
#line 1 "NONE"
	{te = p+1;}
#line 342 "wmkdepend.rl"
	{act = 4;}
	goto st27;
st27:
	if ( ++p == pe )
		goto _test_eof27;
case 27:
#line 771 "wmkdepend.cc"
	if ( (*p) == 10 )
		goto tr2;
	goto st1;
st20:
	if ( ++p == pe )
		goto _test_eof20;
case 20:
	if ( (*p) == 10 )
		goto tr31;
	goto st20;
st0:
cs = 0;
	goto _out;
	}
	_test_eof21: cs = 21; goto _test_eof; 
	_test_eof1: cs = 1; goto _test_eof; 
	_test_eof22: cs = 22; goto _test_eof; 
	_test_eof23: cs = 23; goto _test_eof; 
	_test_eof24: cs = 24; goto _test_eof; 
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
	_test_eof25: cs = 25; goto _test_eof; 
	_test_eof12: cs = 12; goto _test_eof; 
	_test_eof13: cs = 13; goto _test_eof; 
	_test_eof14: cs = 14; goto _test_eof; 
	_test_eof15: cs = 15; goto _test_eof; 
	_test_eof16: cs = 16; goto _test_eof; 
	_test_eof26: cs = 26; goto _test_eof; 
	_test_eof17: cs = 17; goto _test_eof; 
	_test_eof18: cs = 18; goto _test_eof; 
	_test_eof19: cs = 19; goto _test_eof; 
	_test_eof27: cs = 27; goto _test_eof; 
	_test_eof20: cs = 20; goto _test_eof; 

	_test_eof: {}
	if ( p == eof )
	{
	switch ( cs ) {
	case 1: goto tr0;
	case 22: goto tr36;
	case 23: goto tr36;
	case 24: goto tr36;
	case 2: goto tr0;
	case 3: goto tr0;
	case 4: goto tr0;
	case 5: goto tr0;
	case 6: goto tr0;
	case 7: goto tr0;
	case 8: goto tr0;
	case 9: goto tr0;
	case 10: goto tr0;
	case 11: goto tr0;
	case 25: goto tr37;
	case 12: goto tr17;
	case 13: goto tr17;
	case 14: goto tr0;
	case 26: goto tr37;
	case 17: goto tr17;
	case 18: goto tr17;
	case 27: goto tr38;
	}
	}

	_out: {}
	}

#line 428 "wmkdepend.rl"
       /* ^^^ FSM execution here ^^^ */;

        if (0 == cs)
        {
            // Typically only arises when missing a trailing newline
            std::cerr
                << EXENAME ": parse error while scanning '"
                << fileName << "' ... perhaps missing a final newline\n";
            break;
        }

        if (ts)
        {
            // Preserve incomplete token.
            // We have the normal ragel range (ts, te) but potentially
            // our own local buffer start as 'tok'

            if (tok && tok >= ts)
            {
                tok = inbuf + (tok - ts);
            }
            else
            {
                tok = nullptr;          // safety
            }

            pending = pe - ts;
            memmove(inbuf, ts, pending);
            te = inbuf + (te - ts);     // token end (after memmove)
            ts = inbuf;                 // token start
        }
        else
        {
            pending = 0;
        }
    }
    ::fclose(infile);
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        std::cerr
            << EXENAME ": input file not supplied\n";
        return 1;
    }

    unsigned nIncDirs = 0;

    // Prechecks:
    // - help
    // - count number of -I directories
    for (int i = 1; i < argc; ++i)
    {
        if (argv[i][0] != '-') continue;

        switch (argv[i][1])
        {
            case 'h':           // Option: -h, -help
                usage();
                return 0;
                break;

            case 'q':           // Option: -q (quiet)
                optQuiet = true;
                break;

            case 'v':           // Option: -v (verbose)
                optVerbose = true;
                break;

            case 'I':           // Option: -Idir
                ++nIncDirs;
                break;

            // Could check other options, warn about unknown options...
        }
    }

    sourceFile.assign(argv[argc-1]);
    std::string outputFile;

    // Verify that input file has an extension
    {
        const auto dot = sourceFile.find_last_of("./");
        if (dot == std::string::npos || sourceFile[dot] != '.')
        {
            std::cerr
                << EXENAME ": no file extension in source file name '"
                << sourceFile << "'\n";

            return 1;
        }
    }

    Files::reset(nIncDirs);

    // Build list of -I directories and add -i ignores
    for (int i = 1; i < argc; ++i)
    {
        const size_t optLen = strlen(argv[i]);

        if (!strncmp(argv[i], "-I", 2))
        {
            // Option: -Idir
            if (optLen > 2)
            {
                Files::dirs.emplace_back(argv[i] + 2);
            }
        }
        else if (!strncmp(argv[i], "-i", 2))
        {
            // Option: -iheader
            if (optLen > 2)
            {
                Files::visited.insert(argv[i] + 2);
            }
        }
        else if (!strncmp(argv[i], "-e", 2))
        {
            // Option: -eENV
            if (optLen > 2)
            {
                Files::addEnv(argv[i] + 2);
            }
        }
        else if (!strncmp(argv[i], "-o", 2))
        {
            // Option: -oFile
            if (optLen > 2)
            {
                outputFile.assign(argv[i] + 2);
            }
        }
    }


    // Output to an intermediate file
    std::string outputTmp;
    if (outputFile.length())
    {
        outputTmp = outputFile + ".part";

        output = ::fopen(outputTmp.c_str(), "w");

        if (!output)
        {
            std::cerr
                << EXENAME ": could not open file '"
                << outputTmp << "' for output: " << strerror(errno) << '\n';
            return 1;
        }
    }

    fputs("$(OBJECTS_DIR)/", output);
    fputs(sourceFile.c_str(), output);
    fputs(".dep: \\\n", output);

    processFile(sourceFile);

    fputs("\n#END\n", output);
    fflush(output);

    // Atomic move of intermediate file to final output file
    if (outputFile.length())
    {
        ::fclose(output);
        ::rename(outputTmp.c_str(), outputFile.c_str());
    }

    return 0;
}


/*****************************************************************************/
