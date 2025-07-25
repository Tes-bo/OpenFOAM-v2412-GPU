/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2022 OpenCFD Ltd.
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
    Foam::bitSetOrBoolList

Description
    Simple wrapper for handling test() on bitSet or boolList
    without a templating layer or lambda expresssion.

\*---------------------------------------------------------------------------*/

#ifndef Foam_bitSetOrBoolList_H
#define Foam_bitSetOrBoolList_H

#include "bitSet.H"
#include "boolList.H"

/*---------------------------------------------------------------------------*\
                      Class bitSetOrBoolList Declaration
\*---------------------------------------------------------------------------*/

namespace Foam
{

class bitSetOrBoolList
{
    const bitSet& bits_;
    const boolList& bools_;

public:

    // Constructors

        //- Construct with a bitSet reference
        explicit bitSetOrBoolList(const bitSet& select)
        :
            bits_(select),
            bools_(boolList::null())
        {}

        //- Construct with a boolList reference
        explicit bitSetOrBoolList(const boolList& select)
        :
            bits_(bitSet::null()),
            bools_(select)
        {}


    // Member Functions

        //- Is empty
        bool empty() const noexcept
        {
            return bits_.empty() && bools_.empty();
        }

        //- Size
        label size() const noexcept
        {
            return bits_.size() + bools_.size();
        }

        //- Test function
        bool test(const label i) const
        {
            return bits_.test(i) || bools_.test(i);
        }

        //- Test predicate
        bool operator()(const label i) const
        {
            // Can also use test(i) etc...
            return bits_.test(i) || bools_.test(i);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
