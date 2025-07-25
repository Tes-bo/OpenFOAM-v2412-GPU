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

Class
    Foam::FIFOStack

Description
    A FIFO stack based on a singly-linked list.

    Stack operations are push(), pop(), top(), bottom().

\*---------------------------------------------------------------------------*/

#ifndef Foam_FIFOStack_H
#define Foam_FIFOStack_H

#include "SLList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class FIFOStack Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class FIFOStack
:
    public SLList<T>
{
public:

    // Constructors

        //- Inherit constructors from SLList
        using SLList<T>::SLList;


    // Member Functions

        //- Const reference to the top element
        const T& top() const
        {
            return this->back();
        }

        //- Const reference to the bottom element
        const T& bottom() const
        {
            return this->front();
        }

        //- Push an element onto the back of the stack
        void push(const T& elem)
        {
            this->push_back(elem);
        }

        //- Move an element onto the back of the stack
        void push(T&& elem)
        {
            this->push_back(std::move(elem));
        }

        //- Pop the bottom element off the stack
        T pop()
        {
            return this->removeHead();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
