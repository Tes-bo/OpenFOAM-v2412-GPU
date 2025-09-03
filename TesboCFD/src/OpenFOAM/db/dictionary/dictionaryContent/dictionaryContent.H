/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::dictionaryContent

Description
    A wrapper for dictionary content, \em without operators that could
    affect inheritance patterns.

\*---------------------------------------------------------------------------*/

#ifndef Foam_dictionaryContent_H
#define Foam_dictionaryContent_H

#include "dictionary.H"
#include "wordList.H"
#include "wordRes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class dictionaryContent Declaration
\*---------------------------------------------------------------------------*/

class dictionaryContent
{
    // Private Data

        //- The dictionary content
        dictionary dict_;

public:

    // Constructors

        //- Default construct
        dictionaryContent() = default;

        //- Copy construct
        dictionaryContent(const dictionaryContent&) = default;

        //- Move construct
        dictionaryContent(dictionaryContent&&) = default;

        //- Copy construct from dictionary
        explicit dictionaryContent(const dictionary& dict)
        :
            dict_(dict)
        {}

        //- Move construct from dictionary
        explicit dictionaryContent(dictionary&& dict)
        :
            dict_(std::move(dict))
        {}


    //- Destructor
    virtual ~dictionaryContent() = default;


    // Member Functions

        //- Copy construct a dictionary,
        //- filtered by simple allow/deny lists
        //
        //  An empty 'allow' list accepts everything not in the 'deny' list.
        //
        //  \return filtered dictionary copy
        static dictionary copyDict
        (
            const dictionary& input,
            const wordList& allow = wordList(),
            const wordList& deny = wordList()
        );

        //- Copy construct a dictionary,
        //- filtered by a combination of allow/deny lists
        //
        //  An empty 'allow' list accepts everything not in the 'deny' list.
        //  A literal match has higher priority over a regex match.
        //
        //  \verbatim
        //     input:  ( abc apple wall wall1 wall2 )
        //     allow:  ( abc def "wall.*" )
        //     deny:   ( "[ab].*" wall )
        //
        //     result: (abc wall1 wall2)
        //  \endverbatim
        //
        //  \return filtered dictionary copy
        static dictionary copyDict
        (
            const dictionary& input,
            const wordRes& allow,
            const wordRes& deny = wordRes()
        );


        //- Read-access to the content
        const dictionary& dict() const noexcept
        {
            return dict_;
        }

        //- Copy assign new content
        void dict(const dictionary& dict)
        {
            dict_ = dict;
        }

        //- Move assign new content
        void dict(dictionary&& dict)
        {
            dict_ = std::move(dict);
        }


    // Operators

        //- No copy assignment
        void operator=(const dictionaryContent&) = delete;

        //- No move assignment
        void operator=(dictionaryContent&&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
