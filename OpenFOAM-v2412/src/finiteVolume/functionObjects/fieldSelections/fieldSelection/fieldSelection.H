/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2019 OpenCFD Ltd.
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
    Foam::functionObjects::fieldSelection

Description
    Helper class to manage field selections

    The class holds a list of field name filters which are then applied to a
    set of field objects (in derived classes) from which the resulting set is
    available via the selection() function.  This returns a list of
    (fieldName, component) objects, e.g. for U.component(0) this is (U, 0).

SourceFiles
    fieldSelection.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_fieldSelection_H
#define Foam_functionObjects_fieldSelection_H

#include "fieldInfo.H"
#include "DynamicList.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class dictionary;
class objectRegistry;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class fieldSelection Declaration
\*---------------------------------------------------------------------------*/

class fieldSelection
:
    public List<fieldInfo>
{
    // Private Member Functions

        //- Reset the field filters to the given field names
        template<class Container>
        bool resetFieldFiltersImpl(const Container& names);

        //- No copy construct
        fieldSelection(const fieldSelection&) = delete;


protected:

    // Protected Member Data

        //- Reference to the database
        const objectRegistry& obr_;

        //- Flag to indicate whether components are allowed
        const bool includeComponents_;

        //- Current field selection
        List<fieldInfo> selection_;


    // Protected Member Functions

        //- Add registered objects of a given type
        template<class Type>
        void addRegistered(DynamicList<fieldInfo>& set) const;


public:

    //- Construct from object registry
    explicit fieldSelection
    (
        const objectRegistry& obr,
        const bool includeComponents = false
    );


    //- Destructor
    virtual ~fieldSelection() = default;


    // Member Functions

        //- The current field selection
        const List<fieldInfo>& selection() const noexcept
        {
            return selection_;
        }

        //- Return the current filters
        inline HashSet<wordRe> filters() const;

        //- Return the current field selection, in sorted order
        inline wordList selectionNames() const;

        //- Reset the field filters to the given field names
        virtual bool resetFieldFilters(const HashSet<wordRe>& names);

        //- Reset the field filters to the given field name
        virtual bool resetFieldFilters(const wordRe& name);

        //- Read the fieldSelection data from dictionary
        virtual bool read(const dictionary& dict);

        //- Return whether the field names contain a pattern
        virtual bool containsPattern() const;

        //- Clear the current selection
        virtual void clearSelection();

        //- Update the selection
        virtual bool updateSelection();

        //- Check that all requested fielda have been found
        virtual bool checkSelection();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fieldSelectionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fieldSelectionTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
