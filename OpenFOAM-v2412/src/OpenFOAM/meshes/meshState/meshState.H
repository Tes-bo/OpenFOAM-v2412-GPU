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
    Foam::meshState

Description
    Database for mesh data, solution data, solver performance and other
    reduced data.

    polyMesh has a meshState data member so that
    every GeometricField has access.

SourceFiles
    meshState.C
    meshStateTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_meshState_H
#define Foam_meshState_H

#include "IOdictionary.H"
#include "solverPerformance.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class meshState Declaration
\*---------------------------------------------------------------------------*/

class meshState
:
    public IOdictionary
{
    // Private Data

        //- Previously used time-index, used for reset between iterations
        mutable label prevTimeIndex_;


    // Private Member Functions

        //- True if boolean entry exists and is set
        static bool getBoolEntry(const dictionary& dict, const word& keyword);

        //- Add/remove boolean entry
        static void setBoolEntry
        (
            dictionary& dict,
            const word& keyword,
            bool on
        );

        //- No copy construct
        meshState(const meshState&) = delete;

        //- No copy assignment
        void operator=(const meshState&) = delete;


public:

    // Static Data Members

        //- Name for "controls" sub-dictionary
        static const word controlsDictName;

        //- Name for "mesh" sub-dictionary
        static const word meshDictName;

        //- Name for "solver" (solver performance) sub-dictionary
        static const word solverPerformanceDictName;


    // Constructors

        //- Construct for objectRegistry (registered with specified name)
        //- optionally with initial content
        meshState
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary* content = nullptr
        );

        //- Construct for objectRegistry (registered as "data")
        explicit meshState(const objectRegistry& obr)
        :
            meshState("data", obr)
        {}

        //- Construct for objectRegistry (registered with specified name)
        //- copying initial content
        meshState
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict
        )
        :
            meshState(name, obr, &dict)
        {}

        //- Construct for objectRegistry (registered as "data")
        //- copying initial content
        meshState(const objectRegistry& obr, const dictionary& dict)
        :
            meshState("data", obr, &dict)
        {}


    // Member Functions

        //- Reset the dictionary
        void reset(const meshState& ms);

        //- Return the dictionary of controls
        dictionary& controlsDict();
        const dictionary& controlsDict() const;

        //- Return the dictionary of mesh data, typically populated by
        //- the polyMesh::checkXXX functions, checkMesh etc.
        dictionary& meshDict();
        const dictionary& meshDict() const;

        //- Return the dictionary of solver performance data which
        //- includes initial and final residuals for convergence checking
        dictionary& solverPerformanceDict();
        const dictionary& solverPerformanceDict() const;

        //- True if "firstIteration" entry exists (in controls) and is set
        bool isFirstIteration() const;

        //- True if "finalIteration" entry exists (in controls) and is set
        bool isFinalIteration() const;

        //- Add/remove "firstIteration" entry (from controls)
        void setFirstIteration(bool on);

        //- Add/remove "finalIteration" entry (from controls)
        void setFinalIteration(bool on);

        //- Add/set the solverPerformance entry for the named field
        template<class Type>
        void setSolverPerformance
        (
            const word& name,
            const SolverPerformance<Type>& sp
        ) const;

        //- Add/set the solverPerformance entry, using its fieldName
        template<class Type>
        void setSolverPerformance
        (
            const SolverPerformance<Type>& sp
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "meshStateTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
