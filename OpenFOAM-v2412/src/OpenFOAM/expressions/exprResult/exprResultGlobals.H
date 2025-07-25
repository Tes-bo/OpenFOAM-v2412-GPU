/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2018 Bernhard Gschaider
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    Foam::expressions::exprResultGlobals

Description
    A globally available registry of expression results.
    These are currently registered on Time (may change in the future).

Note
    The storage mechanism is similar to a MeshObject, but always stores
    on Time and flushes stale values according to the time index.

SourceFiles
    exprResultGlobals.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_expressions_exprResultGlobals_H
#define Foam_expressions_exprResultGlobals_H

#include "exprResult.H"
#include "autoPtr.H"
#include "HashPtrTable.H"
#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace expressions
{

/*---------------------------------------------------------------------------*\
                      Class exprResultGlobals Declaration
\*---------------------------------------------------------------------------*/

class exprResultGlobals
:
    public regIOobject
{
public:

    // Public Classes

        //- The table of results
        class Table
        :
            public HashPtrTable<exprResult>
        {
        public:

            //- Default construct
            Table() = default;

            //- Copy (clone) construct
            Table(const Table& tbl);

            //- Move construct
            Table(Table&& tbl);

            //- Read construct from stream
            explicit Table(Istream& is);
        };


private:

    // Private Data

        //- The scoped table of results
        HashTable<Table> variables_;

        //- The currently (or previously) used time-index
        label timeIndex_;


    // Private Member Functions

        //- Construct on Time for registry
        explicit exprResultGlobals(const objectRegistry& obr);

        //- Reset all variables
        void reset();

        //- Get or create a table for the scope
        Table& getOrCreateScope(const word& scope)
        {
            return variables_(scope);
        }


public:

    //- Runtime type information
    TypeNameNoDebug("exprResultGlobals");


    // Selectors

        //- Static constructor for singleton
        static exprResultGlobals& New(const objectRegistry& obr);

        //- Static destructor for singleton
        static bool Delete(const objectRegistry& obr);


    //- Destructor
    virtual ~exprResultGlobals() = default;


    // Member Functions

        //- Get an existing table for the namespace
        Table& getNamespace(const word& name);


        //- Return a global variable, if it exists, or a exprResult::null
        const exprResult& get
        (
            const word& name,
            const wordUList& scopes
        ) const;

        //- Add named result to specified scope
        exprResult& addValue
        (
            const word& name,
            const word& scope,
            const exprResult& value,
            const bool overwrite = true
        );

        //- Add named result to specified scope
        exprResult& addValue
        (
            const word& name,
            const word& scope,
            autoPtr<exprResult>&& value,
            const bool overwrite = true
        );

        //- Extract result from dictionary and add to the scope
        //
        //  Controlled by dictionary entries:
        //  - globalName (mandatory)
        //  - globalScope (optional)
        //  - resultType (optional)
        exprResult& addValue
        (
            const dictionary& dict,
            const word& scope = "",
            const bool overwrite = true
        );

        //- Remove named result from specified scope
        //  \return true if result was removed
        bool removeValue
        (
            const word& name,
            const word& scope
        );

        //- Write variables
        virtual bool writeData(Ostream& os) const;

        //- Read variables
        virtual bool readData(Istream& os);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace expressions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
