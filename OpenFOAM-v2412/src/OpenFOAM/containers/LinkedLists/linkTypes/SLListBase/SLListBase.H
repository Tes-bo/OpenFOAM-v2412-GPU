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
    Foam::SLListBase

Description
    Base for singly-linked lists.

    The iterators associated with the list only have a core functionality
    for navigation, with additional functionality to be added by inheriting
    classes. The node iterators always have a node-pointer as the
    first member data, which allows reinterpret_cast from anything else with
    a nullptr as its first data member.
    The nullObject is such an item (with a nullptr data member).

SourceFiles
    SLListBaseI.H
    SLListBase.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_SLListBase_H
#define Foam_SLListBase_H

#include "label.H"
#include "uLabel.H"
#include "stdFoam.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class SLListBase Declaration
\*---------------------------------------------------------------------------*/

class SLListBase
{
public:

    //- The structure for a singly-linked storage node
    struct link
    {
        //- Pointer to next entry in list
        link* next_ = nullptr;

        //- Default construct
        link() noexcept = default;

        //- Deregister the node (after removal)
        void deregister() noexcept { next_ = nullptr; }
    };


private:

    // Private Data

       //- A pointer to the last element.
       //  last_->next_ points to first element, i.e. circular storage
       link* last_ = nullptr;

       //- Number of elements in the list
       label size_ = 0;


protected:

    // Protected Member Functions

        //- Factory method to return an iterator end
        //  Simply reinterprets a NullObject as a SLListBase iterator.
        template<class IteratorType>
        inline static const IteratorType& iterator_end();

        //- Factory method to return an iterator rend
        //  Deleted for SLListBase
        template<class IteratorType>
        static const IteratorType& iterator_rend() = delete;

        //- Return iterator to first item or end-iterator if list is empty
        //  Removes constness which the caller promises to manage.
        template<class IteratorType>
        inline IteratorType iterator_first() const;

        //- Return iterator to last item or end-iterator if list is empty
        //  Removes constness which the caller promises to manage.
        template<class IteratorType>
        inline IteratorType iterator_last() const;

public:

    // Forward Declarations (iterators)

        class iterator;
        friend class iterator;

        class const_iterator;
        friend class const_iterator;


    // Generated Methods

        //- Default construct
        SLListBase() = default;

        //- No copy construct
        SLListBase(const SLListBase&) = delete;

        //- No copy assignment
        void operator=(const SLListBase&) = delete;

        //- Destructor
        ~SLListBase() = default;


    // Member Functions

        //- True if the list is empty
        bool empty() const noexcept { return !size_; }

        //- The number of elements in list
        label size() const noexcept { return size_; }

        //- Return first entry
        inline link* front();

        //- Return const access to first entry
        inline const link* front() const;

        //- Return last entry
        inline link* back();

        //- Return const access to last entry
        inline const link* back() const;


        //- Add at front of list
        void push_front(link* item);

        //- Add at back of list
        void push_back(link* item);

        //- Remove and return first entry
        link* removeHead();

        // Remove and return element
        link* remove(link* item);

        // Remove and return element specified by iterator
        inline link* remove(iterator& iter);

        //- Clear the list
        inline void clear();

        //- Swap the contents of list
        inline void swap(SLListBase& lst);

        //- Transfer the contents of the argument into this list
        //- and annul the argument list.
        inline void transfer(SLListBase& lst);


    // iterator

        //- A primitive non-const node iterator.
        //  Must normally be extended by inheriting classes.
        class iterator
        {
            friend class SLListBase;
            friend class const_iterator;

            //- The selected node.
            //  MUST be the first member for easy comparison between iterators
            //  and for reinterpret_cast from nullObject
            link* node_;

            //- The list being iterated on (as pointer for bitwise copy)
            SLListBase* list_;

            //- Copy of the node next pointer (to use after removal)
            link copy_;

        public:

            //- Copy construct
            iterator(const iterator&) = default;

            //- Construct for a node on the list
            inline iterator(SLListBase* list, link* item);

            //- The storage node
            inline link* get_node() const noexcept;

            //- Pointing at a valid storage node
            inline bool good() const noexcept;

            //- Cannot move backward through list
            inline void prev() = delete;

            //- Move forward through list
            inline void next();

            //- Copy assignment
            inline void operator=(const iterator& iter);

            inline bool operator==(const iterator& iter) const;
            inline bool operator!=(const iterator& iter) const;
        };


    // STL const_iterator

        //- A primitive const node iterator.
        //  Must normally be extended by inheriting classes.
        class const_iterator
        {
            //- The selected node.
            //  MUST be the first member for easy comparison between iterators
            //  and for reinterpret_cast from nullObject
            const link* node_;

            //- The list being iterated on (as pointer for bitwise copy)
            const SLListBase* list_;

        public:

            //- Copy construct
            const_iterator(const const_iterator&) = default;

            //- Construct for a node on the list
            inline const_iterator(const SLListBase* list, const link* item);

            //- Construct from a non-const iterator
            inline const_iterator(const SLListBase::iterator& iter);

            //- The storage node
            inline const link* get_node() const noexcept;

            //- Pointing at a valid storage node
            inline bool good() const noexcept;

            //- Cannot move backward through list
            inline void prev() = delete;

            //- Move forward through list
            inline void next();

            //- Copy assignment
            const_iterator& operator=(const const_iterator&) = default;

            inline bool operator==(const const_iterator& iter) const;
            inline bool operator!=(const const_iterator& iter) const;
        };


        //- Iterator to first item in list with non-const access
        inline iterator begin();

        //- Iterator to first item in list with const access
        inline const_iterator cbegin() const;

        //- No reverse iteration
        const_iterator crbegin() const = delete;

        //- End of list for iterators
        inline const iterator& end();

        //- End of list for iterators
        inline const const_iterator& cend() const;

        //- No reverse iteration
        const const_iterator& crend() const = delete;


    // Housekeeping

        //- Return first entry
        //FOAM_DEPRECATED_FOR(2022-10, "front()")
        link* first() { return front(); }

        //- Return const access to first entry
        //FOAM_DEPRECATED_FOR(2022-10, "front()")
        const link* first() const { return front(); }

        //- Return last entry
        //FOAM_DEPRECATED_FOR(2022-10, "back()")
        link* last() { return back(); }

        //- Return const access to last entry
        //FOAM_DEPRECATED_FOR(2022-10, "back()")
        const link* last() const { return back(); }

        //- Add at front of list
        //FOAM_DEPRECATED_FOR(2022-10, "push_front()")
        void prepend(link* item) { push_front(item); }

        //- Add at back of list
        //FOAM_DEPRECATED_FOR(2022-10, "push_back()")
        void append(link* item) { push_back(item); }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SLListBaseI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
