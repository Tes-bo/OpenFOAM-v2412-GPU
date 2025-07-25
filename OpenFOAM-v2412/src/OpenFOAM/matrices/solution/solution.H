/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020-2023 OpenCFD Ltd.
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
    Foam::solution

Description
    Selector class for relaxation factors, solver type and solution.

    The file will contain these types of entries:

    \par solvers { }
    A sub-dictionary listing of linear solver settings. For example,
    \verbatim
    solvers
    {
        p
        {
            solver          PCG;
            preconditioner  DIC;
            tolerance       1e-06;
            relTol          0.05;
        }
    }
    \endverbatim

    \par select
    This is a rarely used feature to select between different groups of
    settings within the file. If unspecified, the file contents are used
    directly.

    There are currently no other specified entries, but the file will
    often also contain other solution specifications. For example,
    \verbatim
    PIMPLE
    {
        nOuterCorrectors 2;
        nCorrectors     1;
        nNonOrthogonalCorrectors 0;
    }
    \endverbatim

Note
    The solution data are treated as \c READ_MODIFIED even if the
    requested \p readOption is nominally MUST_READ or READ_IF_PRESENT.
    This allows run-time modifications to behave as expected. <br>
    The optional fallback dictionary content for constructors is used
    when a file is missing or for a NO_READ, with a null pointer being
    treated like an empty dictionary.

SourceFiles
    solution.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_solution_H
#define Foam_solution_H

#include "IOdictionary.H"
#include "HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class Type> class Function1;

/*---------------------------------------------------------------------------*\
                          Class solution Declaration
\*---------------------------------------------------------------------------*/

class solution
:
    public IOdictionary
{
    // Private Data

        //- Dictionary of temporary fields to cache
        dictionary cache_;

        //- Switch for the caching mechanism
        bool caching_;

        //- Dictionary of relaxation factors for all the fields
        dictionary fieldRelaxDict_;

        //- Cache of Function1s in above dictionary
        mutable HashPtrTable<Function1<scalar>> fieldRelaxCache_;

        //- Dictionary of relaxation factors for all the equations
        dictionary eqnRelaxDict_;

        //- Cache of Function1s in above dictionary
        mutable HashPtrTable<Function1<scalar>> eqnRelaxCache_;

        //- Optional default relaxation factor for all the fields
        autoPtr<Function1<scalar>> fieldRelaxDefault_;

        //- Optional default relaxation factor for all the equations
        autoPtr<Function1<scalar>> eqnRelaxDefault_;

        //- Dictionary of solver parameters for all the fields
        dictionary solvers_;


    // Private Member Functions

        //- Read settings from the dictionary
        void read(const dictionary&);

        //- The entire dictionary or the optional "select" sub-dictionary.
        const dictionary& selectedDict() const;


    // Generated Methods

        //- No copy construct
        solution(const solution&) = delete;

        //- No copy assignment
        void operator=(const solution&) = delete;


public:

    //- Update from older solver controls syntax
    //  Usually verbose, since we want to know about the changes
    //  Returns the number of settings changed
    static label upgradeSolverDict(dictionary& dict, const bool verbose=true);

    //- Debug switch (registered name: "solution")
    static int debug;


    // Constructors

        //- Construct for objectRegistry, readOption, (system) dictionary name.
        solution
        (
            const objectRegistry& obr,
            IOobjectOption::readOption rOpt,
            const fileName& dictName,
            const dictionary* fallback = nullptr
        );

        //- Construct for objectRegistry, (system) dictionary name
        //- using the readOption from the registry.
        solution
        (
            const objectRegistry& obr,
            const fileName& dictName,
            const dictionary* fallback = nullptr
        );


    //- Destructor. Non-default in header (incomplete types)
    virtual ~solution();


    // Member Functions

    // Access

        //- True if the given field should be cached
        bool cache(const word& name) const;

        //- True if the relaxation factor is given for the field
        bool relaxField(const word& name) const;

        //- Get the relaxation factor specified for the field
        //- or the specified "default" entry, if present.
        //- Does not change \p factor if neither direct nor "default"
        //- can be used,
        //  \return True if found
        bool relaxField(const word& name, scalar& factor) const;

        //- True if the relaxation factor is given for the equation
        bool relaxEquation(const word& name) const;

        //- Get the relaxation factor specified for the equation
        //- or the specified "default" entry, if present.
        //- Does not change \p factor if neither direct nor "default"
        //- can be used,
        //  \return True if found
        bool relaxEquation(const word& name, scalar& factor) const;

        //- Get the relaxation factor for the given field.
        //- Fatal if not found.
        scalar fieldRelaxationFactor(const word& name) const;

        //- Get the relaxation factor for the given equation.
        //- Fatal if not found.
        scalar equationRelaxationFactor(const word& name) const;

        //- The entire dictionary or the optional "select" sub-dictionary.
        const dictionary& solutionDict() const;

        //- Return \p name sub-dictionary within the solutionDict().
        //  Same as \c solutionDict().subDict(...)
        const dictionary& solutionDict(const word& name) const;

        //- The solver controls dictionary (all fields)
        const dictionary& solversDict() const;

        //- The solver controls dictionary for the given field.
        //- Same as \c solversDict().subDict(...)
        const dictionary& solverDict(const word& name) const;

        //- The solver controls dictionary for the given field.
        //- Same as solverDict(...)
        const dictionary& solver(const word& name) const;


    // Read

        //- Read the solution dictionary
        bool read();


    // Other

        //- Helper for printing cache message
        template<class FieldType>
        static void cachePrintMessage
        (
            const char* message,
            const word& name,
            const FieldType& fld
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "solutionTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
