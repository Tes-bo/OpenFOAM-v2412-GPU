/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2021 OpenCFD Ltd.
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
    Foam::function1Base

Description
    Base class for template-invariant parts of Function1

SourceFiles
    function1Base.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_function1Base_H
#define Foam_function1Base_H

#include "dictionary.H"
#include "objectRegistry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class Type> class Function1;

/*---------------------------------------------------------------------------*\
                     Class function1Base Declaration
\*---------------------------------------------------------------------------*/

class function1Base
:
    public refCount
{
    // Private Member Functions

        //- The associated registry, the time registry or nullptr
        const objectRegistry* whichDb(const bool useTime) const noexcept;


protected:

    // Protected Data

        //- Name of entry
        const word name_;

        //- Pointer to an object registry
        const objectRegistry* obrPtr_;


    // Protected Member Functions

        //- No copy assignment
        void operator=(const function1Base&) = delete;


public:

    // Constructors

        //- Construct from entry name and optional registry
        explicit function1Base
        (
            const word& entryName,
            const objectRegistry* obrPtr = nullptr
        );

        //- Construct from entry name, dictionary (unused) and optional registry
        function1Base
        (
            const word& entryName,
            const dictionary& dict,
            const objectRegistry* obrPtr = nullptr
        );

        //- Copy construct
        explicit function1Base(const function1Base& rhs);


    //- Destructor
    virtual ~function1Base();


    // Member Functions

    // Access

        //- The name of the entry
        const word& name() const noexcept
        {
            return name_;
        }

        //- Return the associated registry or nullptr.
        const objectRegistry* whichDb() const noexcept
        {
            return obrPtr_;
        }

        //- Reset the associated objectRegistry
        void resetDb(const objectRegistry* obrPtr = nullptr) noexcept;

        //- Reset the associated objectRegistry
        void resetDb(const objectRegistry& db) noexcept;

        //- Return the object registry
        //  FatalError if object registry is not set
        const objectRegistry& obr() const;

        //- Return true if this function was created with the time database
        bool isTime() const noexcept;

        //- Return the time database
        //  FatalError if object registry is not set
        const Time& time() const;

        //- Return the mesh database if this Function1 was created using a mesh
        //  Note: relies on refCast failure if the type is not correct
        template<class MeshType>
        const MeshType& mesh(const word& regionName = word::null) const
        {
            const objectRegistry* ptr = whichDb(!regionName.empty());

            if (!ptr)
            {
                FatalErrorInFunction
                    << "Object registry not set"
                    << abort(FatalError);
            }

            if (regionName.empty())
            {
                return refCast<const MeshType>(*ptr);
            }
            else
            {
                return ptr->lookupObject<MeshType>(regionName);
            }
        }


    // Manipulation

        //- Convert time
        virtual void userTimeToTime(const Time& t);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
