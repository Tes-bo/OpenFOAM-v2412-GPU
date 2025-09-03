/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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
    Foam::dynamicCodeContext

Description
    Encapsulation of dynamic code dictionaries

SourceFiles
    dynamicCodeContext.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dynamicCodeContext_H
#define Foam_dynamicCodeContext_H

#include "dictionary.H"
#include "SHA1.H"
#include <functional>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class dynamicCodeContext Declaration
\*---------------------------------------------------------------------------*/

class dynamicCodeContext
{
    // Private Data

        //- The parent dictionary context
        std::reference_wrapper<const dictionary> dict_;

        //- The SHA1 of the contents
        SHA1 sha1_;

        //- The "codeOptions" entry (optional)
        string codeOptions_;

        //- The "codeLibs" entry (optional)
        string codeLibs_;

        //- The "codeInclude" entry (optional)
        string codeInclude_;

        //- The "localCode" entry (optional)
        string localCode_;

        //- The "code" entry (optional)
        string code_;


public:

    // Constructors

        //- Default construct
        dynamicCodeContext();

        //- Construct from a dictionary
        explicit dynamicCodeContext(const dictionary& dict);


    // Static Member Functions

        //- Cleanup string and expand with dictionary parameters
        static void inplaceExpand(string& str, const dictionary& dict);

        //- Prefix a \#line directive to code.
        //  The input lineNum is 0-based.
        //  Is a no-op if any of the arguments are invalid
        //  (lineNum is negative, code or file are empty)
        //
        //  \return The change in string length caused by the directive.
        //  This can potentially be used to recover the substring portions.
        static unsigned addLineDirective
        (
            string& code,
            label lineNum,
            const string& file
        );

        //- Prefix a \#line directive to code.
        //  The name of the dictionary is used for the 'file' name.
        static unsigned addLineDirective
        (
            string& code,
            label lineNum,
            const dictionary& dict
        );


    // Member Functions

        //- Not using dummy code context (dictionary::null)
        bool good() const noexcept;

        //- Same as good()
        bool valid() const noexcept { return good(); }

        //- Set code context from a dictionary
        void setCodeContext(const dictionary& dict);

        //- Return the parent dictionary context
        const dictionary& dict() const noexcept
        {
            return dict_.get();
        }

        //- The code options (Make/options)
        const string& options() const noexcept
        {
            return codeOptions_;
        }

        //- The code libs (LIB_LIBS)
        const string& libs() const noexcept
        {
            return codeLibs_;
        }

        //- The code includes
        const string& include() const noexcept
        {
            return codeInclude_;
        }

        //- The local (file-scope) code
        const string& localCode() const noexcept
        {
            return localCode_;
        }

        //- The code
        const string& code() const noexcept
        {
            return code_;
        }

        //- The SHA1 calculated from options, libs, include, code, etc.
        const SHA1& sha1() const noexcept
        {
            return sha1_;
        }

        //- Add content to SHA1 hashing
        void append(const std::string& str)
        {
            sha1_.append(str);
        }


    // Reading

        //- Locate literal dictionary entry, nullptr if not found
        const entry* findEntry(const word& key) const;

        //- Read string entry from context dictionary
        //- append content to SHA1 hashing and add line number etc.
        //
        //  The string is cleared before reading.
        bool readEntry
        (
            const word& key,
            string& str,
            bool mandatory = true,
            bool withLineNum = true
        );

        //- Read optional string entry from context dictionary,
        //- append content to SHA1 hashing and add line number etc.
        //
        //  The string is cleared before reading.
        bool readIfPresent
        (
            const word& key,
            string& str,
            bool withLineNum = true
        );


    // Member Operators

        //- Cast to dictionary
        operator const dictionary&() const noexcept
        {
            return dict_.get();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
