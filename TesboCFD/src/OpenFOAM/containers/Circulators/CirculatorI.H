/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2015 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Container, bool Const>
bool Foam::CirculatorIters<Container, Const>::equal
(
    const CirculatorIters<Container, Const>& rhs
)
{
    return
    (
        begin_ == rhs.begin_
     && end_ == rhs.end_
     && iter_ == rhs.iter_
     && fulcrum_ == rhs.fulcrum_
    );
}


template<class Container, bool Const>
void Foam::CirculatorIters<Container, Const>::increment()
{
    ++iter_;
    if (iter_ == end_)
    {
        iter_ = begin_;
    }
}


template<class Container, bool Const>
void Foam::CirculatorIters<Container, Const>::decrement()
{
    if (iter_ == begin_)
    {
        iter_ = end_;
    }
    --iter_;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Container, bool Const>
Foam::CirculatorIters<Container, Const>::CirculatorIters()
:
    begin_(0),
    end_(0),
    iter_(0),
    fulcrum_(0)
{}


template<class Container, bool Const>
Foam::CirculatorIters<Container, Const>::CirculatorIters
(
    const iterator& begin,
    const iterator& end
)
:
    begin_(begin),
    end_(end),
    iter_(begin_),
    fulcrum_(begin_)
{}


template<class Container, bool Const>
Foam::CirculatorIters<Container, Const>::CirculatorIters
(
    const CirculatorIters<Container, Const>& rhs
)
:
    begin_(rhs.begin_),
    end_(rhs.end_),
    iter_(rhs.iter_),
    fulcrum_(rhs.fulcrum_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Container, bool Const>
bool Foam::CirculatorIters<Container, Const>::empty() const
{
    return (end_ == begin_);
}


template<class Container, bool Const>
typename Foam::CirculatorIters<Container, Const>::size_type
Foam::CirculatorIters<Container, Const>::size() const
{
    return (end_ - begin_);
}


template<class Container, bool Const>
typename Foam::CirculatorIters<Container, Const>::difference_type
Foam::CirculatorIters<Container, Const>::nRotations() const
{
    return (iter_ - fulcrum_);
}


template<class Container, bool Const>
bool Foam::CirculatorIters<Container, Const>::circulate
(
    const CirculatorBase::direction dir
)
{
    if (dir == CirculatorBase::CLOCKWISE)
    {
        increment();
    }
    else if (dir == CirculatorBase::ANTICLOCKWISE)
    {
        decrement();
    }

    return !(iter_ == fulcrum_);
}


template<class Container, bool Const>
void Foam::CirculatorIters<Container, Const>::setFulcrumToIterator()
{
    fulcrum_ = iter_;
}


template<class Container, bool Const>
void Foam::CirculatorIters<Container, Const>::setIteratorToFulcrum()
{
    iter_ = fulcrum_;
}


template<class Container, bool Const>
typename Foam::CirculatorIters<Container, Const>::reference
Foam::CirculatorIters<Container, Const>::curr() const
{
    return *iter_;
}


template<class Container, bool Const>
typename Foam::CirculatorIters<Container, Const>::reference
Foam::CirculatorIters<Container, Const>::next() const
{
    if (iter_ == end_ - 1)
    {
        return *begin_;
    }

    return *(iter_ + 1);
}


template<class Container, bool Const>
typename Foam::CirculatorIters<Container, Const>::reference
Foam::CirculatorIters<Container, Const>::prev() const
{
    if (iter_ == begin_)
    {
        return *(end_ - 1);
    }

    return *(iter_ - 1);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Container, bool Const>
void Foam::CirculatorIters<Container, Const>::operator=
(
    const CirculatorIters<Container, Const>& rhs
)
{
    if (this == &rhs)
    {
        return;  // Self-assignment is a no-op
    }

    begin_ = rhs.begin_;
    end_ = rhs.end_;
    iter_ = rhs.iter_;
    fulcrum_ = rhs.fulcrum_;
}


template<class Container, bool Const>
Foam::CirculatorIters<Container, Const>&
Foam::CirculatorIters<Container, Const>::operator++()
{
    this->increment();
    return *this;
}


template<class Container, bool Const>
Foam::CirculatorIters<Container, Const>
Foam::CirculatorIters<Container, Const>::operator++(int)
{
    auto old(*this);
    this->increment();
    return old;
}


template<class Container, bool Const>
Foam::CirculatorIters<Container, Const>&
Foam::CirculatorIters<Container, Const>::operator--()
{
    this->decrement();
    return *this;
}


template<class Container, bool Const>
Foam::CirculatorIters<Container, Const>
Foam::CirculatorIters<Container, Const>::operator--(int)
{
    auto old(*this);
    this->decrement();
    return old;
}


template<class Container, bool Const>
bool Foam::CirculatorIters<Container, Const>::operator==
(
    const CirculatorIters<Container, Const>& rhs
) const
{
    return this->equal(rhs);
}


template<class Container, bool Const>
bool Foam::CirculatorIters<Container, Const>::operator!=
(
    const CirculatorIters<Container, Const>& rhs
) const
{
    return !this->equal(rhs);
}


template<class Container, bool Const>
typename Foam::CirculatorIters<Container, Const>::reference
Foam::CirculatorIters<Container, Const>::operator*() const
{
    return *iter_;
}


template<class Container, bool Const>
typename Foam::CirculatorIters<Container, Const>::reference
Foam::CirculatorIters<Container, Const>::operator()() const
{
    return *iter_;
}


template<class Container, bool Const>
typename Foam::CirculatorIters<Container, Const>::difference_type
Foam::CirculatorIters<Container, Const>::operator-
(
    const CirculatorIters<Container, Const>& rhs
) const
{
    return (iter_ - rhs.iter_);
}


// ************************************************************************* //
