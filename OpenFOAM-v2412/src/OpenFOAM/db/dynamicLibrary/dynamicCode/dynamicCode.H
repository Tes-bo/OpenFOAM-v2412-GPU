/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
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
    Foam::dynamicCode

Description
    Tools for handling dynamic code compilation

SourceFiles
    dynamicCode.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dynamicCode_H
#define Foam_dynamicCode_H

#include "fileName.H"
#include "DynamicList.H"
#include "HashTable.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class dynamicCodeContext;
class ISstream;
class OSstream;
class SHA1Digest;

/*---------------------------------------------------------------------------*\
                        Class dynamicCode Declaration
\*---------------------------------------------------------------------------*/

class dynamicCode
{
    // Private Data

        //- Root for dynamic code compilation
        fileName codeRoot_;

        //- Subdirectory name for loading libraries
        const fileName libSubDir_;

        //- Name for code
        word codeName_;

        //- Name for code subdirectory
        word codeDirName_;

        //- Files to copy and filter
        DynamicList<fileName> compileFiles_;

        //- Files to copy and filter
        DynamicList<fileName> copyFiles_;

        //- Direct contents for files
        DynamicList<Tuple2<fileName, std::string>> createFiles_;

        //- Variables to use during filtering
        HashTable<string> filterVars_;

        //- Contents for Make/options
        std::string makeOptions_;


protected:

    // Static Data Members

        //- Directory for library targets for Make/files
        static const char* const targetLibDir;

        //- Top-level directory name for copy/compiling
        static const char* const topDirName;


    // Protected Member Functions

        //- Copy lines while expanding variables
        static void copyAndFilter
        (
            ISstream&,
            OSstream&,
            const HashTable<string>& mapping
        );

        //- Resolve code-templates via the codeTemplateEnvName
        //  alternatively in the codeTemplateDirName via Foam::findEtcFile
        static bool resolveTemplates
        (
            const UList<fileName>& templateNames,
            DynamicList<fileName>& resolvedFiles,
            DynamicList<fileName>& badFiles
        );

        //- Write SHA1 value as C-comment
        bool writeCommentSHA1(Ostream&) const;

        //- Copy/create Make/files prior to compilation
        bool createMakeFiles() const;

        //- Copy/create Make/options prior to compilation
        bool createMakeOptions() const;


        //- Write digest to Make/SHA1Digest
        bool writeDigest(const SHA1Digest&) const;

        //- Write digest to Make/SHA1Digest
        bool writeDigest(const std::string&) const;


public:

    // Static Data Members

        //- Name of the code template environment variable
        //  Used to located the codeTemplateName
        static const word codeTemplateEnvName;

        //- Name of the code template sub-directory
        //  Used when locating the codeTemplateName via Foam::findEtcFile
        static const fileName codeTemplateDirName;

        //- Flag if system operations are allowed
        static int allowSystemOperations;


    // Static Member Functions

        //- Check security for creating dynamic code
        static void checkSecurity(const char* title, const dictionary&);


    // Generated Methods

        //- No copy construct
        dynamicCode(const dynamicCode&) = delete;

        //- No copy assignment
        void operator=(const dynamicCode&) = delete;


    // Constructors

        //- Construct for a specified code name and code directory name
        //  Defaults to using the code name for the code directory name
        explicit dynamicCode
        (
            const word& codeName,
            const word& codeDirName = ""
        );


    // Member Functions

        //- Return the code-name
        const word& codeName() const noexcept { return codeName_; }

        //- Return the code-dirname
        const word& codeDirName() const noexcept { return codeDirName_; }

        //- Root for dynamic code compilation
        //  Expanded from \<case\>/dynamicCode
        const fileName& codeRoot() const noexcept { return codeRoot_; }

        //- Subdirectory name for loading libraries
        //  Expanded from platforms/\$WM_OPTIONS/lib
        const fileName& libSubDir() const noexcept { return libSubDir_; }

        //- Path for specified code name
        //  Corresponds to codeRoot()/codeDirName()
        fileName codePath() const
        {
            return codeRoot_/codeDirName_;
        }

        //- Path for specified code name relative to \<case\>
        //  Corresponds to topDirName/codeDirName()
        fileName codeRelPath() const;

        //- Library path for specified code name
        //  Corresponds to codeRoot()/libSubDir()/lib\<codeName\>.so
        fileName libPath() const;

        //- Library path for specified code name relative to \<case\>
        //  Corresponds to
        //  dynamicCode/codeDirName()/libSubDir()/lib\<codeName\>.so
        fileName libRelPath() const;


        //- Path for SHA1Digest
        //  Corresponds to codePath()/Make/SHA1Digest
        fileName digestFile() const
        {
            return codeRoot_/codeDirName_/"Make/SHA1Digest";
        }


        //- Clear files and variables
        void clear();

        //- Clear files and reset variables to specified context
        void reset(const dynamicCodeContext&);


        //- Add a file template name, which will be found and filtered
        void addCompileFile(const fileName& name);

        //- Add a file template name, which will be found and filtered
        void addCopyFile(const fileName& name);

        //- Add a file to create with its contents. Will not be filtered
        void addCreateFile(const fileName& name, const std::string& contents);

        //- Define filter variables for code, codeInclude, SHA1sum
        void setFilterContext(const dynamicCodeContext&);

        //- Define a filter variable
        void setFilterVariable(const word& key, const std::string& value);

        //- Define contents for Make/options
        void setMakeOptions(const std::string& content);


        //- Verify if the copied code is up-to-date, based on Make/SHA1Digest
        bool upToDate(const dynamicCodeContext& context) const;

        //- Verify if the copied code is up-to-date, based on Make/SHA1Digest
        bool upToDate(const SHA1Digest& sha1) const;

        //- Copy/create files prior to compilation
        bool copyOrCreateFiles(const bool verbose = false) const;

        //- Compile a libso
        bool wmakeLibso() const;


    // Synchronisation

        //- Wait for libPath() file to appear on sub-ranks.
        //  Returns immediately if non-parallel or
        //  IOobject::fileModificationSkew is not set
        static void waitForFile
        (
            const fileName& file,
            const dictionary& contextDict
        );


    // Convenience

        //- Define a filter variables TemplateType and FieldType
        template<class Type>
        void setFieldTemplates()
        {
            std::string val(pTraits<Type>::typeName);

            // Template type
            setFilterVariable("TemplateType", val);

            // Field type - eg, ScalarField, VectorField, ...
            val[0] = toupper(val[0]);
            val += "Field";
            setFilterVariable("FieldType", val);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
