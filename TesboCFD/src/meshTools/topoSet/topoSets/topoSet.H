/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2024 OpenCFD Ltd.
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
    Foam::topoSet

Description
    General set of labels of mesh quantity (points, cells, faces).

SourceFiles
    topoSet.C
    topoSetTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_topoSet_H
#define Foam_topoSet_H

#include "HashSet.H"
#include "regIOobject.H"
#include "labelList.H"
#include "typeInfo.H"
#include "autoPtr.H"
#include "pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class mapPolyMesh;
class polyMesh;
class primitiveMesh;
class mapDistributePolyMesh;

/*---------------------------------------------------------------------------*\
                           Class topoSet Declaration
\*---------------------------------------------------------------------------*/

class topoSet
:
    public regIOobject,
    public labelHashSet
{
protected:

    // Protected Member Functions

        //- Read into labelHashSet if IOobject flags set. Return true if read.
        bool readIOcontents(const word& wantedType, labelHashSet& contents);

        //- Check limits on addressable range.
        void checkLabels(const labelUList& labels, const label maxSize);

        //- Check limits on addressable range.
        void checkLabels(const labelHashSet& labels, const label maxSize);

        //- Update map from map.
        //  Used to update cell/face labels after morphing
        virtual void updateLabels(const labelUList& map);

        //- Check limits on addressable range.
        virtual void check(const label maxSize);

        //- Write part of contents nicely formatted.
        //  Prints labels only.
        //  \returns number of labels written
        label writeDebug
        (
            Ostream& os,
            const label maxElem,
            labelHashSet::const_iterator& iter
        ) const;

        //- Write part of contents nicely formatted.
        //  Prints label and corresponding coordinate.
        //  \returns number of labels written
        label writeDebug
        (
            Ostream& os,
            const pointField& coords,
            const label maxElem,
            labelHashSet::const_iterator& iter
        ) const;

        //- Write labels and coordinates columnwise to os. Truncate to maxLen.
        void writeDebug
        (
            Ostream& os,
            const pointField& coords,
            const label maxLen
        ) const;

        //- No copy construct
        topoSet(const topoSet&) = delete;


public:

    //- Runtime type information
    TypeName("topoSet");


    // Static

        //- Debug switch to disallow the use of generic sets
        static int disallowGenericSets;

        //- Name of file set will use.
        static fileName localPath(const polyMesh& mesh, const word& name);

        //- Find IOobject in the polyMesh/sets/ (used as constructor helper)
        static IOobject findIOobject
        (
            const polyMesh& mesh,
            const word& name,
            IOobjectOption::readOption rOpt = IOobjectOption::MUST_READ,
            IOobjectOption::writeOption wOpt = IOobjectOption::NO_WRITE,
            IOobjectOption::registerOption reg = IOobjectOption::LEGACY_REGISTER
        );

        //- Find IOobject in the polyMesh/sets/ (used as constructor helper)
        static IOobject findIOobject
        (
            const polyMesh& mesh,
            const word& name,
            IOobjectOption::registerOption reg,
            IOobjectOption::readOption rOpt = IOobjectOption::MUST_READ,
            IOobjectOption::writeOption wOpt = IOobjectOption::NO_WRITE
        )
        {
            return findIOobject(mesh, name, rOpt, wOpt, reg);
        }

        //- Find IOobject in the polyMesh/sets/ (used as constructor helper)
        static IOobject findIOobject
        (
            const Time& runTime,
            const word& name,
            IOobjectOption::readOption rOpt = IOobjectOption::MUST_READ,
            IOobjectOption::writeOption wOpt = IOobjectOption::NO_WRITE,
            IOobjectOption::registerOption reg = IOobjectOption::LEGACY_REGISTER
        );


    // Declare run-time constructor selection table

        // For the direct constructor
        declareRunTimeSelectionTable
        (
            autoPtr,
            topoSet,
            word,
            (
                const polyMesh& mesh,
                const word& name,
                IOobjectOption::readOption rOpt,
                IOobjectOption::writeOption wOpt
            ),
            (mesh, name, rOpt, wOpt)
        );

        // For the constructor from size
        declareRunTimeSelectionTable
        (
            autoPtr,
            topoSet,
            size,
            (
                const polyMesh& mesh,
                const word& name,
                const label size,
                IOobjectOption::writeOption wOpt
            ),
            (mesh, name, size, wOpt)
        );

        // For the constructor as copy
        declareRunTimeSelectionTable
        (
            autoPtr,
            topoSet,
            set,
            (
                const polyMesh& mesh,
                const word& name,
                const topoSet& set,
                IOobjectOption::writeOption wOpt
            ),
            (mesh, name, set, wOpt)
        );


    // Constructors

        //- Construct from IOobject as explicitly passed type.
        //  Can't use typeName info here since subclasses not yet instantiated
        topoSet(const IOobject& io, const word& wantedType);

        //- Construct empty (no-read) with IOobject information
        topoSet(const IOobject& io, const Foam::zero);

        //- Construct empty (no-read) with initial labelHashSet capacity,
        //- with IOobject information
        topoSet(const IOobject& io, const label initialCapacity);

        //- Copy construct (no-read) from labelHashSet,
        //- with IOobject information.
        topoSet(const IOobject& io, const labelHashSet& labels);

        //- Move construct (no-read) from labelHashSet,
        //- with IOobject information.
        topoSet(const IOobject& io, labelHashSet&& labels);

        //- Construct from polyMesh and name.
        //  Searches for a polyMesh/sets directory but not beyond the
        //  mesh.facesInstance().
        topoSet
        (
            const polyMesh& mesh,
            const word& wantedType,
            const word& name,
            IOobjectOption::readOption rOpt = IOobjectOption::MUST_READ,
            IOobjectOption::writeOption wOpt = IOobjectOption::NO_WRITE,
            IOobjectOption::registerOption reg = IOobjectOption::LEGACY_REGISTER
        );

        //- Construct empty (no-read) with initial labelHashSet capacity.
        //  Searches for a polyMesh/sets directory but not beyond the
        //  mesh.facesInstance().
        topoSet
        (
            const polyMesh& mesh,
            const word& name,
            const label initialCapacity,
            IOobjectOption::writeOption wOpt = IOobjectOption::NO_WRITE,
            IOobjectOption::registerOption reg = IOobjectOption::LEGACY_REGISTER
        );

        //- Copy construct (no-read) from labelHashSet.
        //  Searches for a polyMesh/sets directory but not beyond the
        //  mesh.facesInstance().
        topoSet
        (
            const polyMesh& mesh,
            const word& name,
            const labelHashSet& labels,
            IOobjectOption::writeOption wOpt = IOobjectOption::NO_WRITE,
            IOobjectOption::registerOption reg = IOobjectOption::LEGACY_REGISTER
        );

        //- Move construct (no-read) from labelHashSet.
        //  Searches for a polyMesh/sets directory but not beyond the
        //  mesh.facesInstance().
        topoSet
        (
            const polyMesh& mesh,
            const word& name,
            labelHashSet&& labels,
            IOobjectOption::writeOption wOpt = IOobjectOption::NO_WRITE,
            IOobjectOption::registerOption reg = IOobjectOption::LEGACY_REGISTER
        );

        //- Copy construct (no-read) from list of labels.
        //  Searches for a polyMesh/sets directory but not beyond the
        //  mesh.facesInstance().
        topoSet
        (
            const polyMesh& mesh,
            const word& name,
            const labelUList& labels,
            IOobjectOption::writeOption wOpt = IOobjectOption::NO_WRITE,
            IOobjectOption::registerOption reg = IOobjectOption::LEGACY_REGISTER
        );


        //- Clone
        autoPtr<topoSet> clone() const
        {
            NotImplemented;
            return nullptr;
        }


    // Selectors

        //- Return a pointer to a toposet read from file
        static autoPtr<topoSet> New
        (
            const word& setType,
            const polyMesh& mesh,
            const word& name,
            IOobjectOption::readOption rOpt = IOobjectOption::MUST_READ,
            IOobjectOption::writeOption wOpt = IOobjectOption::NO_WRITE
        );

        //- Return a pointer to a new toposet of given size
        static autoPtr<topoSet> New
        (
            const word& setType,
            const polyMesh& mesh,
            const word& name,
            const label size,
            IOobjectOption::writeOption wOpt = IOobjectOption::NO_WRITE
        );

        //- Return a pointer to a new toposet as copy of another toposet
        static autoPtr<topoSet> New
        (
            const word& setType,
            const polyMesh& mesh,
            const word& name,
            const topoSet& set,
            IOobjectOption::writeOption wOpt = IOobjectOption::NO_WRITE
        );


    //- Destructor
    virtual ~topoSet() = default;


    // Member Functions

        //- Has the given index?
        virtual bool contains(const label id) const;

        //- Has the given index?
        virtual bool found(const label id) const;

        //- Set an index
        virtual bool set(const label id);

        //- Unset an index
        virtual bool unset(const label id);

        //- Set multiple indices
        virtual void set(const labelUList& labels);

        //- Unset multiple indices
        virtual void unset(const labelUList& labels);

        //- Invert contents.
        //  Insert all members [0,maxLen) which were not in set.
        virtual void invert(const label maxLen);

        //- Subset contents. Only elements present in both sets remain.
        virtual void subset(const labelUList& elems);

        //- Subset contents. Only elements present in both sets remain.
        virtual void subset(const topoSet& set);

        //- Add given elements to the set
        virtual void addSet(const labelUList& elems);

        //- Add given elements to the set
        virtual void addSet(const topoSet& set);

        //- Subtract given elements from the set
        virtual void subtractSet(const labelUList& elems);

        //- Subtract given elements from the set
        virtual void subtractSet(const topoSet& set);


        //- Sync set across coupled patches.
        virtual void sync(const polyMesh& mesh);

        //- Write labels columnwise to os. Truncate to maxLen.
        virtual void writeDebug(Ostream& os, const label maxLen) const;

        //- Like above but also writes mesh related quantity
        //- (usually coordinate).
        virtual void writeDebug
        (
            Ostream& os,
            const primitiveMesh&,
            const label maxLen
        ) const = 0;

        //- Write contents.
        virtual bool writeData(Ostream&) const;

        //- Update any stored data for new labels. Not implemented.
        virtual void updateMesh(const mapPolyMesh& morphMap);

        //- Update any stored data for mesh redistribution.
        virtual void distribute(const mapDistributePolyMesh&) = 0;

        //- Return max allowable index (+1). Not implemented.
        virtual label maxSize(const polyMesh& mesh) const = 0;

        //- Helper: call updateMesh on all items in container
        //- (and updates instance)
        template<class Container>
        static void updateMesh
        (
            const fileName& instance,
            const mapPolyMesh&,
            Container& items
        );

        //- Helper: set instance on all items in container
        template<class Container>
        static void setInstance
        (
            const fileName& instance,
            Container& items
        );

        //- Helper: remove all sets files from mesh instance
        static void removeFiles(const polyMesh&);


    // Member Operators

        //- Copy labelHashSet part only
        void operator=(const topoSet&);


    // Housekeeping

        //- Deprecated(2018-10) subtract elements present in set.
        //  \deprecated(2018-10) - use subtractSet instead
        FOAM_DEPRECATED_FOR(2018-10, "subtractSet()")
        virtual void deleteSet(const topoSet& set) { this->subtractSet(set); }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "topoSetTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
