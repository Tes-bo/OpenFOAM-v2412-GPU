/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::DirLister

Description
    A class for listing directory contents, or obtaining a list of
    directory contents.

    For many situations, the iterator-based access is the most convenient
    since it allows use of a range-for without any intermediate variables
    or allocations.

    List all directories or files:
    \code
    for (const word& item : DirLister("."))
    {
        Info<< "dir or file: " << item << nl;
    }
    \endcode

    List files or directories only:
    \code
    for (const word& item : DirLister(".").files())
    {
        Info<< "file: " << item << nl;
    }
    for (const word& item : DirLister(".").dirs())
    {
        Info<< "dir: " << item << nl;
    }
    \endcode
    These can be combined with a unary predicate to select only specific
    items:
    \code
    {
        wordRes matchProcs
        {
            wordRe("processors"),
            wordRe("processor\\d+", wordRe::REGEX)
        };

        for
        (
            const word& item
          : DirLister::dirs(".").where(matchProcs)
        )
        {
           Info<< "processor dir: " << item << nl;
        }
    }
    \endcode
    The unary predicate can similarly be used with the list or sorted
    methods:
    \code
    {
        wordRes matchProcs
        {
            wordRe("processor[0-9][0-9]*", wordRe::REGEX),
            wordRe("processors")
        };

        fileNameList procDirs
        (
            DirLister::dirs(".").csorted<fileName>(matchProcs)
        );
    }
    \endcode

See Also
    Foam::readDir

SourceFiles
    DirListerI.H
    DirLister.cxx
    DirLister.txx

\*---------------------------------------------------------------------------*/

#ifndef Foam_DirLister_H
#define Foam_DirLister_H

#include "fileNameList.H"
#include "wordList.H"
#include <dirent.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class DirLister Declaration
\*---------------------------------------------------------------------------*/

class DirLister
{
    // Private data

        fileName dirName_;

        fileName::Type requestedType_;

        bool followLink_;

        bool stripgz_;


protected:

    // Protected Member Functions

        //- Placeholder for additional selection filters
        virtual inline bool accept(const word& name) const { return true; }

        //- Traverse to next entry
        word next(DIR* dirPtr) const;


public:

    //- Specialization for filtering with a predicate
    template<class UnaryPredicate> class Filtered;


    // Constructors

        //- Construct from components - list files and directories
        inline DirLister
        (
            const fileName& directory,
            const fileName::Type = fileName::UNDEFINED,
            const bool filtergz = true,
            const bool followLink = true
        );


        //- Convenience constructor for listing directories
        inline static DirLister dirs
        (
            const fileName& directory,
            const bool followLink = true
        );

        //- Convenience constructor for listing files
        inline static DirLister files
        (
            const fileName& directory,
            const bool filtergz = true,
            const bool followLink = true
        );

        //- Convenience for constructing a filtered iterator
        template<class UnaryPredicate>
        inline Filtered<UnaryPredicate> where
        (
            const UnaryPredicate& pred,
            const bool prune = false
        ) const;


    //- Destructor
    virtual ~DirLister() = default;


    // Member Functions

        //- Return the directory name
        const fileName& dirName() const noexcept { return dirName_; }

        //- Return a complete list of names
        template<class StringType=Foam::word>
        List<StringType> list() const;

        //- Return complete list of names
        template<class StringType=Foam::word, class UnaryPredicate>
        List<StringType> list
        (
            const UnaryPredicate& pred,
            const bool prune = false
        ) const;

        //- Return a complete list of names, sorted in natural order
        template<class StringType=Foam::word>
        List<StringType> csorted() const;

        //- Return complete list of names, sorted in natural order
        template<class StringType=Foam::word, class UnaryPredicate>
        List<StringType> csorted
        (
            const UnaryPredicate& pred,
            const bool prune = false
        ) const;


    // Iterators

        //- Const iterator for traversing directory contents
        class const_iterator
        {
            const DirLister* lister_;
            DIR* dirptr_;
            word name_;

            //- Open directory and set dirptr_ (nullptr on failure)
            bool open(const fileName& dir);

            //- Close dirptr_
            void close();

            bool next();

        public:

            //- Construct end iterator
            inline const_iterator();

            //- Construct begin iterator
            inline const_iterator(const DirLister* lister);

            //- Destructor
            inline ~const_iterator();

            inline const_iterator& operator++();
            const word& operator*() const noexcept { return name_; }

            inline bool operator==(const const_iterator& iter) const;
            inline bool operator!=(const const_iterator& iter) const;
        };


        inline const_iterator cbegin() const;
        inline const_iterator cend() const;

        inline const_iterator begin() const;
        inline const_iterator end() const;

};


//- A DirLister specialization for filtering with a predicate
template<class UnaryPredicate>
class DirLister::Filtered
:
    public DirLister
{
    const UnaryPredicate& pred_;
    bool prune_;

protected:

    friend class DirLister;

    virtual inline bool accept(const word& name) const
    {
        return (pred_(name) ? !prune_ : prune_);
    }

    //- Constructor
    inline Filtered
    (
        const DirLister& lister,
        const UnaryPredicate& pred,
        const bool prune = false
    );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DirListerI.H"
#include "DirLister.txx"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
