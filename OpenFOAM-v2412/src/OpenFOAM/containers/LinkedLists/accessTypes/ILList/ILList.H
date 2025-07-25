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
    Foam::ILList

Description
    Template class for intrusive linked lists.

SourceFiles
    ILList.C
    ILListIO.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_ILList_H
#define Foam_ILList_H

#include "UILList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class LListBase, class T> class ILList;

template<class LListBase, class T> Istream& operator>>
(
    Istream& is,
    ILList<LListBase, T>& lst
);


/*---------------------------------------------------------------------------*\
                           Class ILList Declaration
\*---------------------------------------------------------------------------*/

template<class LListBase, class T>
class ILList
:
    public UILList<LListBase, T>
{
    // Private Member Functions

        //- Read from Istream using given Istream constructor class
        template<class INew>
        void readIstream(Istream& is, const INew& inew);


public:

    // Constructors

        //- Default construct
        ILList() = default;

        //- Construct and add initial item pointer
        explicit ILList(T* item)
        :
            UILList<LListBase, T>(item)
        {}

        //- Construct from Istream
        explicit ILList(Istream& is);

        //- Copy construct using the 'clone()' method for each element
        ILList(const ILList<LListBase, T>& lst);

        //- Move construct
        ILList(ILList<LListBase, T>&& lst);

        //- Copy constructor with additional argument for clone 'clone()'
        template<class CloneArg>
        ILList(const ILList<LListBase, T>& lst, const CloneArg& cloneArg);

        //- Construct from Istream using given Istream constructor class
        template<class INew>
        ILList(Istream& is, const INew& inew);


    //- Destructor. Calls clear()
    ~ILList();


    // Member Functions

        //- Clear the contents of the list
        void clear();

        //- Remove first element(s) from the list (deletes pointers)
        void pop_front(label n = 1);

        //- Remove the specified element from the list and delete it
        bool erase(T* item);

        //- Transfer the contents of the argument into this List
        //- and annul the argument list.
        void transfer(ILList<LListBase, T>& lst);


    // Member Operators

        //- Copy assignment using the 'clone()' method for each element
        void operator=(const ILList<LListBase, T>& lst);

        //- Move assignment
        void operator=(ILList<LListBase, T>&& lst);


    // Istream Operator

        //- Read from Istream, discarding existing contents.
        friend Istream& operator>> <LListBase, T>
        (
            Istream& is,
            ILList<LListBase, T>& list
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ILList.C"
    #include "ILListIO.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
