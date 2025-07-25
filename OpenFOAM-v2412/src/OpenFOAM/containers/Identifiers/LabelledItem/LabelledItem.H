/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2021-2022 OpenCFD Ltd.
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
    Foam::LabelledItem

Description
    A container with an integer index that can be attached to any item.
    The index may be useful for sorting or storing additional information.

SeeAlso
    Foam::objectHit
    Foam::PointIndexHit

SourceFiles
    LabelledItemI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_LabelledItem_H
#define Foam_LabelledItem_H

#include "label.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class T> class LabelledItem;
template<class T> Istream& operator>>(Istream&, LabelledItem<T>&);
template<class T> Ostream& operator<<(Ostream&, const LabelledItem<T>&);

/*---------------------------------------------------------------------------*\
                        Class LabelledItem Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class LabelledItem
:
    public T
{
    // Private Data

        //- The object index
        label index_;


public:


    // Constructors

        //- Default construct item, with index = -1
        LabelledItem()
        :
            T(),
            index_(-1)
        {}

        //- Copy construct item, with index = -1
        explicit LabelledItem(const T& item)
        :
            T(item),
            index_(-1)
        {}

        //- Move construct item, with index = -1
        explicit LabelledItem(T&& item)
        :
            T(std::move(item)),
            index_(-1)
        {}

        //- Construct from components
        LabelledItem(const T& item, label idx)
        :
            T(item),
            index_(idx)
        {}

        //- Construct from Istream
        explicit LabelledItem(Istream& is)
        {
            is >> *this;
        }


    // Member Functions

        //- Return the index
        label index() const noexcept
        {
            return index_;
        }

        //- Non-const access to the index
        label& index() noexcept
        {
            return index_;
        }

        //- Set the index
        void setIndex(const label idx) noexcept
        {
            index_ = idx;
        }


    // Member Operators

        //- Test for equality of components
        bool operator==(const LabelledItem<T>& rhs) const
        {
            return
            (
                index_ == rhs.index_
             && static_cast<const T&>(*this) == static_cast<const T&>(rhs)
            );
        }

        //- Test for inequality of components
        bool operator!=(const LabelledItem<T>& rhs) const
        {
            return !(*this == rhs);
        }


    // IOstream Operators

        friend Istream& operator>> <T>(Istream&, LabelledItem<T>&);
        friend Ostream& operator<< <T>(Ostream&, const LabelledItem<T>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class T>
inline Foam::Istream& Foam::operator>>
(
    Istream& is,
    LabelledItem<T>& item
)
{
    is.readBegin("LabelledItem");
    is >> static_cast<T&>(item) >> item.index();
    is.readEnd("LabelledItem");

    is.check(FUNCTION_NAME);
    return is;
}


template<class T>
inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const LabelledItem<T>& item
)
{
    // Output like Tuple2
    os  << token::BEGIN_LIST
        << static_cast<const T&>(item) << token::SPACE
        << item.index()
        << token::END_LIST;

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
