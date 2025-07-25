/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

Class
    Foam::functionEntries::codeStream

Description
    Dictionary entry that contains C++ OpenFOAM code that is compiled to
    generate the entry itself. So
    - codeStream reads three entries: 'code', 'codeInclude' (optional),
    'codeOptions' (optional)
    and uses those to generate library sources inside \c codeStream/
    - these get compiled using 'wmake libso'
    - the resulting library is loaded in executed with as arguments
    \code
        (const dictionary& dict, Ostream& os)
    \endcode
    where the dictionary is the current dictionary.
    - the code has to write into Ostream which is then used to construct
    the actual dictionary entry.


    E.g. to set the internal field of a field:

    \verbatim
    internalField  #codeStream
    {
        code
        #{
            const IOdictionary& d = static_cast<const IOdictionary&>(dict);
            const fvMesh& mesh = refCast<const fvMesh>(d.db());
            scalarField fld(mesh.nCells(), 12.34);
            fld.writeEntry("", os);
        #};

        //- Optional:
        codeInclude
        #{
            #include "fvCFD.H"
        #};

        //- Optional:
        codeOptions
        #{
            -I$(LIB_SRC)/finiteVolume/lnInclude
        #};
    };
    \endverbatim


    Note the \c \#{ ... \c \#} syntax is a 'verbatim' input mode that allows
    inputting strings with embedded newlines.

    Limitations:
    - '~' symbol not allowed inside the code sections.
    - probably some other limitations (uses string::expand which expands
      \c \$ and \c ~ sequences)

Note
    The code to be compiled is stored under the local \c codeStream directory
    with a subdirectory name corresponding to the SHA1 of the contents.

    The corresponding library code is located under the local
    \c codeStream/platforms/$WM_OPTIONS/lib directory in a library
    \c libcodeStream_SHA1.so

SourceFiles
    codeStream.C

\*---------------------------------------------------------------------------*/

#ifndef functionEntries_codeStream_H
#define functionEntries_codeStream_H

#include "functionEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class dlLibraryTable;

namespace functionEntries
{

// Forward Declarations
class calcEntry;

/*---------------------------------------------------------------------------*\
                         Class codeStream Declaration
\*---------------------------------------------------------------------------*/

class codeStream
:
    public functionEntry
{
protected:

        //- Interpreter function type
        typedef void (*streamingFunctionType)(Ostream&, const dictionary&);


    // Protected Member Functions

        //- Helper: access IOobject for master-only-reading functionality
        static bool doingMasterOnlyReading(const dictionary& dict);

        //- Helper function: access to dlLibraryTable of Time
        static dlLibraryTable& libs(const dictionary& dict);

        //- Construct, compile, load and return streaming function
        static streamingFunctionType getFunction
        (
            const dictionary& parentDict,
            const dictionary& codeDict
        );


        //- Evaluate dynamically compiled code, returning result as string
        static string evaluate(const dictionary& parentDict, Istream& is);


public:

    //- Name of the C code template to be used
    static constexpr const char* const codeTemplateC = "codeStreamTemplate.C";

    //- Runtime type information
    ClassName("codeStream");


    // Member Functions

        //- Execute in a primitiveEntry context
        static bool execute
        (
            const dictionary& parentDict,
            primitiveEntry& entry,
            Istream& is
        );

        //- Execute in a sub-dict context
        static bool execute(dictionary& parentDict, Istream& is);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionEntries
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
