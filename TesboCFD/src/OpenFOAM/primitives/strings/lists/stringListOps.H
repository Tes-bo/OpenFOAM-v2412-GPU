/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd.
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

InNamespace
    Foam

Description
    Operations on lists of strings.

Namespace
    Foam::stringListOps

Description
    Various utility functions to work on lists of strings.

SourceFiles
    stringListOpsTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_stringListOps_H
#define Foam_stringListOps_H

#include "labelList.H"
#include "stringList.H"
#include "wordRes.H"
#include "ops.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    //- Find first list item that matches, -1 on failure
    template<class UnaryMatchPredicate, class StringType>
    label firstMatchingString
    (
        const UnaryMatchPredicate& matcher,
        const UList<StringType>& input,
        const bool invert=false
    );


    //- Extract list indices for all matches.
    //  The unary match predicate has the following signature:
    //  \code
    //  bool operator()(const std::string& text);
    //  \endcode
    //
    //  \return List indices for matching strings
    //  \see IndirectList::subset_if
    template<class UnaryMatchPredicate, class StringType>
    labelList findMatchingStrings
    (
        const UnaryMatchPredicate& matcher,
        const UList<StringType>& input,
        const bool invert=false
    );


    //- Return list indices for strings matching the regular expression
    //  Template partial specialization of findMatchingStrings
    template<class StringType>
    labelList findStrings
    (
        const regExp& matcher,
        const UList<StringType>& input,
        const bool invert=false
    )
    {
        return findMatchingStrings(matcher, input, invert);
    }


    //- Return list indices for strings matching the regular expression
    //  Template partial specialization of findMatchingStrings
    template<class StringType>
    labelList findStrings
    (
        const keyType& matcher,
        const UList<StringType>& input,
        const bool invert=false
    )
    {
        return
        (
            matcher.isPattern()
          ? findMatchingStrings(regExp(matcher), input, invert)
          : findMatchingStrings(matcher, input, invert)
        );
    }


    //- Return list indices for strings matching the regular expression
    //  Template partial specialization of findMatchingStrings
    template<class StringType>
    labelList findStrings
    (
        const wordRe& matcher,
        const UList<StringType>& input,
        const bool invert=false
    )
    {
        return findMatchingStrings(matcher, input, invert);
    }


    //- Return list indices for strings matching one of the regular expression
    //  Template partial specialization of findMatchingStrings
    template<class StringType>
    labelList findStrings
    (
        const wordRes& matcher,
        const UList<StringType>& input,
        const bool invert=false
    )
    {
        return findMatchingStrings(matcher, input, invert);
    }

    //- Return list indices for strings matching one of the regular expression
    //  Template partial specialization of findMatchingStrings
    template<class StringType>
    labelList findStrings
    (
        const UList<wordRe>& patterns,
        const UList<StringType>& input,
        const bool invert=false
    )
    {
        return findMatchingStrings(wordRes::matcher(patterns), input, invert);
    }


    // Subsetting multi-string matches (similar to ListOp):

    //- Extract elements of StringList when regular expression matches
    //  optionally invert the match
    //  eg, to extract all selected elements:
    //  \code
    //    subsetMatchingStrings<regExp, stringList>(myRegExp, list);
    //  \endcode
    //  \see IndirectList::subset_if
    template<class UnaryMatchPredicate, class StringListType>
    StringListType subsetMatchingStrings
    (
        const UnaryMatchPredicate& matcher,
        const StringListType& input,
        const bool invert=false
    );


    //- Extract elements of StringList when regular expression matches
    //  Template partial specialization of subsetMatchingStrings
    template<class StringListType>
    StringListType subsetStrings
    (
        const regExp& matcher,
        const StringListType& input,
        const bool invert=false
    )
    {
        return subsetMatchingStrings(matcher, input, invert);
    }


    //- Extract elements of StringList when regular expression matches
    //  Template partial specialization of subsetMatchingStrings
    template<class StringListType>
    StringListType subsetStrings
    (
        const keyType& matcher,
        const StringListType& input,
        const bool invert=false
    )
    {
        return
        (
            matcher.isPattern()
          ? subsetMatchingStrings(regExp(matcher), input, invert)
          : subsetMatchingStrings(matcher, input, invert)
        );
    }

    //- Extract elements of StringList when regular expression matches
    //  Template partial specialization of subsetMatchingStrings
    template<class StringListType>
    StringListType subsetStrings
    (
        const wordRe& matcher,
        const StringListType& input,
        const bool invert=false
    )
    {
        return subsetMatchingStrings(matcher, input, invert);
    }

    //- Extract elements of StringList when regular expression matches
    //  Template partial specialization of subsetMatchingStrings
    template<class StringListType>
    StringListType subsetStrings
    (
        const wordRes& matcher,
        const StringListType& input,
        const bool invert=false
    )
    {
        return subsetMatchingStrings(matcher, input, invert);
    }


    //- Extract elements of StringList when regular expression matches
    //  Template partial specialization of subsetMatchingStrings
    template<class StringListType>
    StringListType subsetStrings
    (
        const UList<wordRe>& patterns,
        const StringListType& input,
        const bool invert=false
    )
    {
        return subsetMatchingStrings(wordRes::matcher(patterns), input, invert);
    }


    //- Inplace extract elements of StringList when regular expression matches
    //  optionally invert the match
    //  eg, to extract all selected elements:
    //    inplaceSubsetMatchingStrings<regExp, stringList>(myRegExp, lst);
    template<class UnaryMatchPredicate, class StringListType>
    void inplaceSubsetMatchingStrings
    (
        const UnaryMatchPredicate& matcher,
        StringListType& input,
        const bool invert=false
    );

    //- Inplace extract elements of StringList when regular expression matches
    //  Template partial specialization of inplaceSubsetMatchingStrings
    template<class StringListType>
    void inplaceSubsetStrings
    (
        const regExp& matcher,
        StringListType& input,
        const bool invert=false
    )
    {
        inplaceSubsetMatchingStrings(matcher, input, invert);
    }

    //- Extract elements of StringList when regular expression matches
    //  Template partial specialization of subsetMatchingStrings
    template<class StringListType>
    void inplaceSubsetStrings
    (
        const keyType& matcher,
        StringListType& input,
        const bool invert=false
    )
    {
        return
        (
            matcher.isPattern()
          ? inplaceSubsetMatchingStrings(regExp(matcher), input, invert)
          : inplaceSubsetMatchingStrings(matcher, input, invert)
        );
    }

    //- Inplace extract elements of StringList when regular expression matches
    //  Template partial specialization of inplaceSubsetMatchingStrings
    template<class StringListType>
    void inplaceSubsetStrings
    (
        const wordRe& matcher,
        StringListType& input,
        const bool invert=false
    )
    {
        inplaceSubsetMatchingStrings(matcher, input, invert);
    }

    //- Inplace extract elements of StringList when regular expression matches
    //  Template partial specialization of inplaceSubsetMatchingStrings
    template<class StringListType>
    void inplaceSubsetStrings
    (
        const wordRes& matcher,
        StringListType& input,
        const bool invert=false
    )
    {
        inplaceSubsetMatchingStrings(matcher, input, invert);
    }

    //- Inplace extract elements of StringList when regular expression matches
    //  Template partial specialization of inplaceSubsetMatchingStrings
    template<class StringListType>
    void inplaceSubsetStrings
    (
        const UList<wordRe>& regexs,
        StringListType& input,
        const bool invert=false
    )
    {
        inplaceSubsetMatchingStrings(wordRes::matcher(regexs), input, invert);
    }


/*---------------------------------------------------------------------------*\
                      Namespace stringListOps Declaration
\*---------------------------------------------------------------------------*/

namespace stringListOps
{

//- Functor to determine if a string is exists in a list of strings.
//  For example,
//
//  \code
//  reduce(text, stringListOps::foundOp<word>(myNames));
//  \endcode
template<class StringType>
struct foundOp
{
    const UList<StringType>& values;

    foundOp(const UList<StringType>& list)
    :
        values(list)
    {}

    bool operator()(const std::string& text) const
    {
        return values.found(text);
    }
};


//- Return ids for items with matching names.
//  The filter predicate is a combination of allow and deny lists
//
//  \return List indices for matches
template<class StringListType, class AccessOp = identityOp>
labelList findMatching
(
    const StringListType& input,
    const wordRes::filter& pred,
    AccessOp aop = identityOp()
);

//- Return ids for items with matching names.
//  Uses a combination of allow and deny lists
//
//  An empty 'allow' accepts everything not in 'deny'.
//  A literal 'allow' match has higher priority than any 'deny'.
//  A regex 'allow' match has lower priority than any 'deny'.
//
//  Example (when applied to a list of words),
//  \verbatim
//     input:  ( abc apple test other val val1 val2 wall wall1 wall2 )
//     allow:  ( abc def "t.*" other val val1 "wall.*" )
//     deny:   ( "[ab].*" "t.*" other "val[0-9]" wall )
//
//     result:  (abc other val val1 wall1 wall2)
//  \endverbatim
//
//  \return List indices for matches
template<class StringListType, class AccessOp = identityOp>
labelList findMatching
(
    const StringListType& input,
    const wordRes& allow,
    const wordRes& deny = wordRes::null(),
    AccessOp aop = identityOp()
);

} // End namespace stringListOps


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Housekeeping

namespace Foam
{
    //- Deprecated(2018-02) find using C-string as a regex
    //  \deprecated(2018-02) Treating string as regex may be inefficient
    //      and lead to unintended results.
    //      Use regExp, keyType, wordRe instead, or findMatchingStrings()
    template<class StringType>
    labelList findStrings
    (
        const char* disallowed,
        const UList<StringType>& input,
        const bool invert=false
    ) = delete;

    //- Deprecated(2018-02) find using string as a regex
    //  \deprecated(2018-02) Treating string as regex may be inefficient
    //      and lead to unintended results.
    //      Use regExp, keyType, wordRe instead, or findMatchingStrings()
    template<class StringType>
    labelList findStrings
    (
        const std::string& disallowed,
        const UList<StringType>& input,
        const bool invert=false
    ) = delete;

    //- Deprecated(2018-02) subset using C-string as a regex
    //  \deprecated(2018-02) Treating string as regex may be inefficient
    //      and lead to unintended results.
    //      Use regExp, keyType, wordRe instead, or subsetMatchingStrings()
    template<class StringListType>
    StringListType subsetStrings
    (
        const char* disallowed,
        const StringListType& input,
        const bool invert=false
    ) = delete;

    //- Deprecated(2018-02) subset using string as a regex
    //  \deprecated(2018-02) Treating string as regex may be inefficient
    //      and lead to unintended results.
    //      Use regExp, keyType, wordRe instead, or subsetMatchingStrings()
    template<class StringListType>
    StringListType subsetStrings
    (
        const std::string& disallowed,
        const StringListType& input,
        const bool invert=false
    ) = delete;

    //- Deprecated(2018-02) subset using C-string as a regex
    //  \deprecated(2018-02) Treating string as regex may be inefficient
    //      and lead to unintended results.
    //      Use regExp, keyType, wordRe instead, or inplaceSubsetMatchingStrings()
    template<class StringListType>
    void inplaceSubsetStrings
    (
        const char* disallowed,
        StringListType& input,
        const bool invert=false
    ) = delete;

    //- Deprecated(2018-02) subset using string as a regex
    //  \deprecated(2018-02) Treating string as regex may be inefficient
    //      and lead to unintended results.
    //      Use keyType, wordRe instead, or inplaceSubsetMatchingStrings()
    template<class StringListType>
    void inplaceSubsetStrings
    (
        const std::string& disallowed,
        StringListType& input,
        const bool invert=false
    ) = delete;

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "stringListOpsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
