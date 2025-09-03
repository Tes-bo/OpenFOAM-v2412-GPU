/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2012 OpenFOAM Foundation
    Copyright (C) 2016-2024 OpenCFD Ltd.
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

Namespace
    Foam::stringOps

Description
    Collection of static functions for various string-related operations

SourceFiles
    stringOps.C
    stringOpsTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_stringOps_H
#define Foam_stringOps_H

#include "scalar.H"
#include "SubStrings.H"
#include "dictionary.H"
#include "HashTable.H"
#include "stringOpsSort.H"
#include "stringOpsEvaluate.H"
#include "word.H"
#include "wordRes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class OSstream;
template<class T1, class T2> class Tuple2;

/*---------------------------------------------------------------------------*\
                        Namespace stringOps Declaration
\*---------------------------------------------------------------------------*/

namespace stringOps
{
    //- Count the number of occurrences of the specified character
    //  Correctly handles nullptr.
    std::string::size_type count(const char* s, const char c);

    //- Count the number of occurrences of the specified character
    std::string::size_type count(const std::string& s, const char c);

    //- True if text matches one of the selector expressions
    inline bool match(const UList<wordRe>& selectors, const std::string& text)
    {
        return wordRes::match(selectors, text);
    }

    //- Quote any meta-characters in given string
    template<class StringType, class UnaryPredicate>
    StringType quotemeta
    (
        const StringType& str,
        const UnaryPredicate& meta,
        const char quote = '\\'
    );

    //- Expand occurrences of variables according to the mapping
    //- and return the expanded string.
    //
    //  \sa stringOps::inplaceExpand() for details
    string expand
    (
        const std::string& s,
        const HashTable<string>& mapping,
        const char sigil = '$'
    );

    //- Inplace expand occurrences of variables according to the mapping.
    //- Does \b not use environment values.
    //
    //  Expansion includes:
    //  -# Variables
    //    - \c $VAR
    //    - \c ${VAR}
    //
    //  -# Default and alternative values as per the POSIX shell:
    //     - \c ${parameter:-defValue}
    //       If parameter is unset or null, the \c defValue is substituted.
    //       Otherwise, the value of parameter is substituted.
    //     - \c ${parameter:+altValue}
    //       If parameter is unset or null, nothing is substituted.
    //       Otherwise the \c altValue is substituted.
    //  .
    //
    // General behaviour:
    //  - Unknown entries are removed silently.
    //  - Malformed entries (eg, brace mismatch, sigil followed by unknown
    //    characters) are left as is.
    //
    //  \param[in,out] s The string to modify inplace.
    //  \param mapping The lookup table
    //  \param sigil The leading sigil. Can be changed to avoid conflict
    //      with other string expansions. (default: '$')
    //
    // \par Expansion behaviour
    // - alternatives = True
    // - environment = False
    // - allow empty = True
    // - subDict = Not applicable
    // .
    //
    // \note Currently only used by Foam::dynamicCode.
    void inplaceExpand
    (
        std::string& s,
        const HashTable<string>& mapping,
        const char sigil = '$'
    );


    //- Inplace expand occurrences of variables according to the dictionary
    //- and (optionally) environment variables.
    //
    //  Expansion includes:
    //  -# Dictionary variables and (optionally) environment variables
    //    - \c $VAR
    //    - \c ${VAR}
    //    - \c ${VAR:-defValue}
    //    - \c ${VAR:+altValue}
    //  -# Mathematical evaluation using stringOps::evaluate
    //    - \c ${{EXPR}}
    //  -# Current directory
    //    - leading "./"
    //          : the current directory - Foam::cwd()
    //  -# Leading tag expansion for commonly used directories
    //    - <b> \<etc\>/ </b>
    //          : user/group/other OpenFOAM etc directory
    //    - <b> \<etc:</b><em>[ugoa]+</em>)<b>\>/ </b>
    //          : user/group/other etc with specified location mode
    //    - <b> \<case\>/ </b>
    //          : The \c $FOAM_CASE directory
    //    - <b> \<constant\>/ </b>
    //          : The \c $FOAM_CASE/constant directory
    //    - <b> \<system\>/ </b>
    //          : The \c $FOAM_CASE/system directory
    //  -# Tilde expansion
    //    - leading "~/" : home directory
    //    - leading "~user" : home directory for specified user
    //  -# Default and alternative values as per the POSIX shell:
    //     - \c ${parameter:-defValue}
    //       If parameter is unset or null, the \c defValue is substituted.
    //       Otherwise, the value of parameter is substituted.
    //     - \c ${parameter:+altValue}
    //       If parameter is unset or null, nothing is substituted.
    //       Otherwise the \c altValue is substituted.
    //  .
    //
    // General behaviour:
    //  - Malformed entries (eg, brace mismatch, sigil followed by unknown
    //    characters) are left as is.
    //  - Supports recursive variable expansions.
    //    For example, "${var${num}}" and "${{100 + ${var}}}"
    //
    //  \param[in,out] s The string to modify inplace
    //  \param dict The dictionary context for the expansion
    //  \param allowEnv Allow use of environment values as fallback
    //  \param allowEmpty Allow empty expansions, or treat as Fatal
    //  \param allowSubDict Allow expansion of subDict entries as well as
    //      primitive entries (default: false)
    //  \param sigil The leading sigil. Can be changed to avoid conflict
    //      with other string expansions. (default: '$')
    //
    //  \sa Foam::findEtcEntry(), Foam::findEtcEntries(), stringOps::evaluate()
    //
    // \par Expansion behaviour
    // - alternatives = True
    // - environment = Given by parameter
    // - allow empty = Given by parameter
    // - subDict = Given by parameter (default: False)
    // .
    //
    //  \note This function has too many parameters and should generally
    //      be avoided in user coding.
    void inplaceExpand
    (
        std::string& s,
        const dictionary& dict,
        const bool allowEnv,
        const bool allowEmpty,
        const bool allowSubDict = false,
        const char sigil = '$'
    );


    //- Expand occurrences of dictionary or environment variables.
    //
    //  Empty expansions are allowed.
    //  Serialization of subDict entries is permitted.
    //
    //  \sa stringOps::inplaceExpand(std::string&, const dictionary&, char)
    string expand
    (
        const std::string& s,
        const dictionary& dict,
        const char sigil = '$'
    );

    //- Inplace expand occurrences of dictionary or environment variables.
    //
    //  Empty expansions are allowed.
    //  Serialization of subDict entries is permitted.
    //
    //  \sa
    //  stringOps::inplaceExpand
    //  (std::string&, const dictionary&, bool, bool, bool, char)
    //
    // \par Expansion behaviour
    // - alternatives = True
    // - environment = True
    // - allow empty = True
    // - subDict = True
    // .
    void inplaceExpand
    (
        std::string& s,
        const dictionary& dict,
        const char sigil = '$'
    );


    //- Expand initial tags, tildes, and all occurrences of environment
    //- variables.
    //
    //  \sa
    //  stringOps::inplaceExpand(std::string&, bool);
    string expand(const std::string& s, const bool allowEmpty = false);


    //- Expand initial tags, tildes, and all occurrences of environment
    //- variables
    //
    //  The expansion behaviour is identical to
    //  stringOps::inplaceExpand
    //  (std::string&, const dictionary&, bool, bool, bool, char)
    //  except that there is no dictionary and the environment variables
    //  are always enabled.
    //
    // \par Expansion behaviour
    // - alternatives = True
    // - environment = True
    // - allow empty = Given by parameter (default: False)
    // - subDict = Not applicable
    // .
    void inplaceExpand(std::string& s, const bool allowEmpty = false);


    //- Replace environment variable contents with its name.
    //  This is essentially the inverse operation for inplaceExpand
    //  for a single element.
    //  Return true if a replacement was successful.
    bool inplaceReplaceVar(std::string& s, const word& varName);

    //- Return a copy of the input string with validated characters
    template<class StringType, class UnaryPredicate>
    StringType validate
    (
        const std::string& str,
        const UnaryPredicate& accept,
        const bool invert=false   //!< Invert the test logic
    );

    //- Find (first, last) non-space locations in string or sub-string.
    //  This may change to std::string_view in the future.
    std::pair<size_t, size_t>
    findTrim
    (
        const std::string& s,
        size_t pos = 0,
        size_t len = std::string::npos
    );

    //- Return string trimmed of leading whitespace
    string trimLeft(const std::string& s);

    //- Trim leading whitespace inplace
    void inplaceTrimLeft(std::string& s);

    //- Return string trimmed of trailing whitespace
    string trimRight(const std::string& s);

    //- Trim trailing whitespace inplace
    void inplaceTrimRight(std::string& s);

    //- Return string trimmed of leading and trailing whitespace
    string trim(const std::string& s);

    //- Trim leading and trailing whitespace inplace
    void inplaceTrim(std::string& s);

    //- Eliminate whitespace inplace
    void inplaceRemoveSpace(std::string& s);


    //- Return string with C/C++ comments removed
    string removeComments(const std::string& s);

    //- Remove C/C++ comments inplace
    void inplaceRemoveComments(std::string& s);


    //- Return string copy transformed with std::tolower on each character
    string lower(const std::string& s);

    //- Inplace transform string with std::tolower on each character
    void inplaceLower(std::string& s);

    //- Return string copy transformed with std::toupper on each character
    string upper(const std::string& s);

    //- Inplace transform string with std::toupper on each character
    void inplaceUpper(std::string& s);

    //- Split out arguments (named or unnamed) from an input string.
    //
    //  For example,
    //  \verbatim
    //  (U)
    //      -> named = ()
    //      -> unnamed = (U)
    //
    //  (patch=inlet, p)
    //      -> named = ((patch inlet))
    //      -> unnamed = (p)
    //
    //  testing, start=100, stop=200
    //      -> named = ((start 100)(stop 200))
    //      -> unnamed = (testing)
    //  \endverbatim
    //
    //  \return total number of arguments
    label splitFunctionArgs
    (
        const std::string& str,
        wordRes& args,
        List<Tuple2<word, string>>& namedArgs
    );

    //- Split string into sub-strings at the delimiter character.
    //  Empty sub-strings are normally suppressed.
    //  Behaviour is ill-defined if delim is a NUL character.
    template<class StringType>
    Foam::SubStrings<StringType> split
    (
        //! The string to split
        const StringType& str,
        //! The delimiter for splitting. Ill-defined if NUL character
        const char delim,
        //! Offset within string to start splitting
        std::string::size_type pos = 0,
        //! Retain empty fields
        const bool keepEmpty = false
    );

    //- Split string into sub-strings using delimiter string.
    //  Empty sub-strings are normally suppressed.
    template<class StringType>
    Foam::SubStrings<StringType> split
    (
        //! The string to split
        const StringType& str,
        //! The delimiters for splitting. Ill-defined if empty
        const std::string& delim,
        //! Offset within string to start splitting
        std::string::size_type pos = 0,
        //! Retain empty fields
        const bool keepEmpty = false
    );

    //- Split string into sub-strings using any characters in delimiter.
    //  Empty sub-strings are normally suppressed.
    //  Behaviour is ill-defined if delim is an empty string.
    template<class StringType>
    Foam::SubStrings<StringType> splitAny
    (
        //! The string to split
        const StringType& str,
        //! The delimiters for splitting. Ill-defined if empty!
        const std::string& delim,
        //! Offset within string to start splitting
        std::string::size_type pos = 0
    );

    //- Split string into sub-strings using a fixed field width.
    //  Behaviour is ill-defined if width is zero.
    template<class StringType>
    Foam::SubStrings<StringType> splitFixed
    (
        //! The string to split
        const StringType& str,
        //! Fixed field width for each sub-string
        const std::string::size_type width,
        //! Offset within string to start splitting
        std::string::size_type pos = 0
    );

    //- Split string into sub-strings at whitespace (TAB, NL, VT, FF, CR, SPC)
    //  Empty sub-strings are suppressed.
    template<class StringType>
    Foam::SubStrings<StringType> splitSpace
    (
        //! The string to split
        const StringType& str,
        //! Offset within string to start splitting
        std::string::size_type pos = 0
    );

    //- Output string with text wrapping.
    //  Always includes a trailing newline, unless the string itself is empty.
    void writeWrapped
    (
        //! The output stream
        OSstream& os,
        //! The text to be output
        const std::string& str,
        //! The max-width before wrapping
        const std::string::size_type width,
        //! Indentation for continued lines
        const std::string::size_type indent = 0,
        //! Escape any backslashes on output
        const bool escape = false
    );

} // End namespace stringOps

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "stringOpsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
