%{
/*---------------------------------*- C -*-----------------------------------*\
 =========                   |
 \\      /   F ield          | OpenFOAM: The Open Source CFD Toolbox
  \\    /    O peration      |
   \\  /     A nd            | www.openfoam.com
    \\/      M anipulation   |
------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2020 OpenCFD Ltd.
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
    wmkdep

Description
    A fast dependency list generator that emulates the behaviour and output
    of cpp -M. However, the output contains no duplicates and is thus
    approx. 40% faster than cpp.
    It also handles missing files somewhat more gracefully.

    The algorithm uses flex to scan for includes and searches the files found.
    The files are only visited once (the names of the files visited are hashed),
    which makes this faster than cpp.

Usage
    wmkdep [-Idir..] [-iheader...] [-eENV...] [-oFile] [-q] [-v] filename

\*---------------------------------------------------------------------------*/
/* With cpp:
 *
 * cpp -x c++ -std=c++11 -nostdinc -nostdinc++
 *     -M -DWM_$(WM_PRECISION_OPTION) -DWM_LABEL_SIZE=$(WM_LABEL_SIZE) |
 * sed -e 's,^$(WM_PROJECT_DIR)/,$$(WM_PROJECT_DIR)/,' \
 *     -e 's,^$(WM_THIRD_PARTY_DIR)/,$$(WM_THIRD_PARTY_DIR)/,'
 */

#define FILE_STACK_SIZE 300
#define HASH_TABLE_SIZE 500

#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>  /* POSIX */
#include <dirent.h>     /* POSIX */

/* The executable name (for messages), without requiring access to argv[] */
#define EXENAME  "wmkdep"
#undef yywrap           /* sometimes a macro by default */
#pragma GCC diagnostic ignored "-Wunused-function"

void nextFile(const char* fileName);

/*---------------------------------------------------------------------------*/
%}

%x CMNT CFNAME FFNAME
%%

"//".*\n            ;               /* Remove C++-style comment */

"/*"                BEGIN(CMNT);    /* Begin removing C-style comment */
<CMNT>.|\n          ;
<CMNT>"*/"          BEGIN(INITIAL); /* End removing C-style comment */

^[ \t]*#[ \t]*include[ \t]+\"   BEGIN(CFNAME);  /* C-file name */
<CFNAME>[^"\n ]*    { BEGIN(INITIAL); nextFile(yytext); } /* "-quoted */

"      "include[ \t]+\'   BEGIN(FFNAME);        /* FORTRAN-file name */
<FFNAME>[^']*       { BEGIN(INITIAL); nextFile(yytext); } /* '-quoted */

.|\t|\n             ;

%%

/*---------------------------------------------------------------------------*/


/*
 * A char* entry in hash table
 */
struct HashEntry
{
    char* name;
    struct HashEntry* next;
};


/*
 * Lookup name in hashTable.
 * if found - return 1
 * if not found - insert in hashTable and return 0
 */
static int lookUp(struct HashEntry* hashTable[HASH_TABLE_SIZE], const char* str)
{
    unsigned idx = 0; /* Hash index */
    {
        const char* pp = str;
        while (*pp) idx = idx << 1 ^ *pp++;
        idx %= HASH_TABLE_SIZE;
    }

    /* Search for entry */
    struct HashEntry* entry;
    for (entry = hashTable[idx]; entry; entry = entry->next)
    {
        if (!strcmp(str, entry->name))
        {
            return 1; /* True: entry found */
        }
    }

    /* Entry not found - insert a new entry */
    entry = (struct HashEntry*)malloc(sizeof(struct HashEntry));
    entry->name = strdup(str);
    entry->next = hashTable[idx];
    hashTable[idx] = entry;

    return 0;   /* False: entry did not previously exist */
}


/*
 * Free allocated memory in hashTable.
 */
static void free_hashTable(struct HashEntry* hashTable[HASH_TABLE_SIZE])
{
    int idx;
    for (idx = 0; idx < HASH_TABLE_SIZE; ++idx)
    {
        struct HashEntry* entry = hashTable[idx];
        hashTable[idx] = NULL;

        while (entry)
        {
            struct HashEntry* next = entry->next;
            free(entry->name);
            free(entry);

            entry = next;
        }
    }
}


/*
 * Environment entry - as a linked-list
 */
struct KeyValue
{
    char* name;
    char* value;
    size_t len;
    struct KeyValue* next;
};


/* List of environ variables to substitute */
struct KeyValue* envTable = NULL;

/*
 * Add envTable replacements:
 *
 * Eg,
 *     /openfoam/project/path/directory/xyz
 *  -> $(WM_PROJECT_DIR)/directory/xyz
 */
static void add_env(const char* key)
{
    const char *val = getenv(key);

    if (val && *val)
    {
        const size_t keyLen = strlen(key);
        const size_t valLen = strlen(val);

        /* "$(ENV)/" */
        char *replace = (char*)malloc(keyLen + 5);
        strcpy(replace, "$(");
        strcat(replace, key);
        strcat(replace, ")/");

        /* "/env/value/" */
        char *orig = (char*)malloc(valLen + 2);
        strcpy(orig, val);
        if (val[valLen-1] != '/')
        {
            strcat(orig, "/");
        }

        struct KeyValue* entry =
            (struct KeyValue*)malloc(sizeof(struct KeyValue));

        entry->name  = replace;
        entry->value = orig;
        entry->len   = strlen(orig);
        entry->next  = envTable;
        envTable = entry;
    }
}


/*
 * Free allocated memory in envTable.
 */
static void free_envTable()
{
    struct KeyValue* entry = envTable;

    while (entry)
    {
        struct KeyValue* next = entry->next;
        free(entry->name);
        free(entry->value);
        free(entry);

        entry = next;
    }
}


/*
 * Print fileName to stdout,
 * with envTable substitutions at the beginning of the path
 *
 * Eg,
 *     /openfoam/project/path/directory/xyz
 *  -> $(WM_PROJECT_DIR)/directory/xyz
 */
static void print_fileName(const char* fileName)
{
    const size_t len = strlen(fileName);
    const char *substr = fileName;

    struct KeyValue* entry = envTable;
    while (entry)
    {
        if (len > entry->len && !strncmp(fileName, entry->value, entry->len))
        {
            substr = (fileName + entry->len);
            fputs(entry->name, stdout);
            break;
        }
        entry = entry->next;
    }

    fputs(substr, stdout);
    fputs(" \\\n", stdout);
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int optQuiet = 0;
int optVerbose = 0;
int nDirectories = 0;
char** directories = NULL;
char* sourceFile = NULL;

/* Set of files already visited */
struct HashEntry* visitedFiles[HASH_TABLE_SIZE];

/* Buffer pointer stack counter */
int currentBuffer = 0;

/* Buffer pointer stack */
YY_BUFFER_STATE buffers[FILE_STACK_SIZE];

/* Directory paths for the loaded files */
const char* bufferPaths[FILE_STACK_SIZE];


int main(int argc, char* argv[])
{
    int i;

    if (argc < 2)
    {
        fputs(EXENAME ": input file not supplied\n", stderr);
        return 1;
    }

#if 0
    for (i = 0; i < argc; ++i)
    {
        if (i) fputs(" ", stderr);
        fputs(argv[i], stderr);
    }
#endif

    /* Prechecks:
     * - help
     * - count number of -I directories
     */
    nDirectories = 0;
    for (i = 1; i < argc; ++i)
    {
        if (argv[i][0] != '-') continue;

        switch (argv[i][1])
        {
            case 'h':           /* Option: -h, -help */
                fputs
                (
                    "\nUsage: " EXENAME
                    " [-Idir...] [-iheader...] [-eENV...] [-oFile] [-q] [-v]"
                    " filename\n\n"
                    "  -Idir     Directories to be searched for headers.\n"
                    "  -iheader  Headers to be ignored.\n"
                    "  -eENV     Environment variable path substitutions.\n"
                    "  -oFile    Write output to File.\n"
                    "  -q        Suppress 'No such file' warnings.\n"
                    "  -v        Report each include file to stderr.\n"
                    "\nDependency list generator, similar to 'cpp -M'\n\n",
                    stderr
                );
                return 0;
                break;

            case 'q':           /* Option: -q (quiet) */
                ++optQuiet;
                break;

            case 'v':           /* Option: -v (verbose) */
                ++optVerbose;
                break;

            case 'I':           /* Option: -Idir */
                ++nDirectories;
                break;

            /* Could check other options, warn about unknown options... */
        }
    }

    sourceFile = strdup(argv[argc-1]);

    /* Verify that it has an extension */
    {
        char *base = strrchr(sourceFile, '/');
        if (!base)
        {
            base = sourceFile;
        }
        if (!strrchr(base, '.'))
        {
            fprintf
            (
                stderr,
                EXENAME ": cannot find extension in source file name '%s'\n",
                sourceFile
            );
            exit(1);
        }
    }

    const char *outputFile = NULL;

    /* Build list of -I directories and add -i ignores */
    directories = (char**)malloc(sizeof(char*)*nDirectories);
    nDirectories = 0;
    for (i = 1; i < argc; ++i)
    {
        const size_t optLen = strlen(argv[i]);

        if (!strncmp(argv[i], "-I", 2))
        {
            /* Option: -Idir */
            if (optLen > 2)
            {
                directories[nDirectories++] = strdup(argv[i] + 2);
            }
        }
        else if (!strncmp(argv[i], "-i", 2))
        {
            /* Option: -iheader */
            if (optLen > 2)
            {
                lookUp(visitedFiles, (argv[i] + 2));
            }
        }
        else if (!strncmp(argv[i], "-e", 2))
        {
            /* Option: -eENV */
            if (optLen > 2)
            {
                add_env(argv[i] + 2);
            }
        }
        else if (!strncmp(argv[i], "-o", 2))
        {
            /* Option: -oFile */
            if (optLen > 2)
            {
                outputFile = (argv[i] + 2);
            }
        }
    }

    /* Initialize buffer path for currentBuffer */
    currentBuffer = 0;
    bufferPaths[currentBuffer] = NULL;

    /* Start of output */
    if (outputFile)
    {
        FILE *reopened = freopen(outputFile, "w", stdout);
        if (!reopened)
        {
            fprintf
            (
                stderr,
                EXENAME ": could not open file '%s' for output: %s\n",
                outputFile, strerror(errno)
            );
            fflush(stderr);

            return 1;
        }
    }

    fputs("$(OBJECTS_DIR)/", stdout);
    fputs(sourceFile, stdout);
    fputs(".dep: \\\n", stdout);

    nextFile(sourceFile);
    yylex();

    fputs("\n#END\n", stdout);
    fflush(stdout);

    for (i = nDirectories-1; i >= 0; --i)
    {
        free(directories[i]);
    }
    free(directories);
    free(sourceFile);

    free_hashTable(visitedFiles);
    free_envTable();

    return 0;
}


/*
 * Open a file for reading and print its qualified name
 */
static FILE* fopen_file(const char* dirName, const char* fileName)
{
    FILE *file;

    if (dirName && *dirName)
    {
        const size_t dirLen = strlen(dirName);
        char* fullName = (char*)malloc(dirLen + strlen(fileName) + 2);

        strcpy(fullName, dirName);
        if (dirName[dirLen-1] != '/')
        {
            strcat(fullName, "/");
        }
        strcat(fullName, fileName);

        file = fopen(fullName, "r");

        if (!file && errno == EMFILE)
        {
            fprintf
            (
                stderr,
                EXENAME ": too many open files while opening '%s'\n"
                "Please change your open descriptor limit\n",
                fullName
            );
        }

        if (file)
        {
            print_fileName(fullName);
        }

        free(fullName);
    }
    else
    {
        file = fopen(fileName, "r");

        if (!file && errno == EMFILE)
        {
            fprintf
            (
                stderr,
                EXENAME ": too many open files while opening '%s'\n"
                "Please change your open descriptor limit\n",
                fileName
            );
        }

        if (file)
        {
            print_fileName(fileName);
        }
    }

    if (file && optVerbose)
    {
        fputs(fileName, stderr);
        fputs("\n", stderr);
    }

    return file;
}


/*
 * Open a file and create buffer and put on stack
 */
void nextFile(const char* fileName)
{
    if (lookUp(visitedFiles, fileName))
    {
        return;   /* Already existed (did not insert) */
    }

    if (currentBuffer >= FILE_STACK_SIZE)
    {
        fprintf
        (
            stderr,
            EXENAME ": depth of file search exceeds stack size %d "
            "while opening '%s' for file '%s'\n",
            FILE_STACK_SIZE, fileName, sourceFile
        );

        fflush(stdout);
        fflush(stderr);
        exit(1);
    }

    /* Pointer to new file which is set if the file is successfully opened */
    FILE* newyyin = NULL;

    /* Check if the file has same path as the file in the current buffer */
    if (bufferPaths[currentBuffer] != NULL)
    {
        newyyin = fopen_file(bufferPaths[currentBuffer], fileName);
        if (newyyin)
        {
            buffers[currentBuffer++] = YY_CURRENT_BUFFER;
            bufferPaths[currentBuffer] = bufferPaths[currentBuffer-1];

            yy_switch_to_buffer(yy_create_buffer(newyyin, YY_BUF_SIZE));
            return;
        }
    }

    newyyin = fopen_file(NULL, fileName);
    if (newyyin)
    {
        buffers[currentBuffer++] = YY_CURRENT_BUFFER;
        bufferPaths[currentBuffer] = NULL;

        yy_switch_to_buffer(yy_create_buffer(newyyin, YY_BUF_SIZE));
    }
    else
    {
        int d;
        for (d=0; d<nDirectories; ++d)
        {
            newyyin = fopen_file(directories[d], fileName);
            if (newyyin)
            {
                buffers[currentBuffer++] = YY_CURRENT_BUFFER;
                bufferPaths[currentBuffer] = directories[d];

                yy_switch_to_buffer(yy_create_buffer(newyyin, YY_BUF_SIZE));
                return;
            }
        }

        if (!optQuiet)
        {
            fprintf
            (
                stderr,
                EXENAME ": could not open file '%s' for source file '%s'",
                fileName, sourceFile
            );
            if (nDirectories)
            {
                fprintf(stderr, ": %s", strerror(errno));
            }
            fputs("\n", stderr);
            fflush(stderr);
        }

        fflush(stdout);

        /* Only report the first occurrence */
        lookUp(visitedFiles, fileName);
    }
}


/*
 * The lexer calls yywrap to handle EOF conditions
 */
int yywrap()
{
    /* Close the file for the buffer which has just reached EOF */
    /* This causes strange problems on several systems:
    fclose(yyin);
    yyin = NULL;
    */

    /* Delete the buffer */
    yy_delete_buffer(YY_CURRENT_BUFFER);

    /* Set buffer counter to previous buffer */
    currentBuffer--;

    if (currentBuffer >= 0)
    {
        /* Buffer counter refers to a valid file */

        /* Reset input buffer to the previous buffer on the stack */
        yy_switch_to_buffer(buffers[currentBuffer]);

        /* Return to the normal state for the previous buffer on the stack */
        BEGIN(INITIAL);

        /* Return 0 to inform lex to continue reading */
        return 0;
    }

    /* No more buffers on the stack:
     * Return 1 to inform lex to finish now that all buffers have been read
     */
    return 1;
}


/*****************************************************************************/
