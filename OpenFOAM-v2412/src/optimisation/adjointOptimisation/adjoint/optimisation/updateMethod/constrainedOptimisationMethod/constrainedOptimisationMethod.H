/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2021 PCOpt/NTUA
    Copyright (C) 2013-2021 FOSS GP
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::constrainedOptimisationMethod

Description
    Abstract base class for optimisation methods supporting constraints.
    Does not add functionality to updateMethod but act as a means
    to categorize constrained optimisation methods for performing sanity checks

SourceFiles
    constrainedOptimisationMethod.C

\*---------------------------------------------------------------------------*/

#ifndef constrainedOptimisationMethod_H
#define constrainedOptimisationMethod_H

#include "runTimeSelectionTables.H"
#include "fvMesh.H"
#include "dictionary.H"
#include "designVariables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class constrainedOptimisationMethod Declaration
\*---------------------------------------------------------------------------*/

class constrainedOptimisationMethod
{
private:

    // Private Member Functions

        //- No copy construct
        constrainedOptimisationMethod
        (
            const constrainedOptimisationMethod&
        ) = delete;

        //- No copy assignment
        void operator=(const constrainedOptimisationMethod&) = delete;


public:

    //- Runtime type information
    TypeName("constrainedOptimisationMethod");

    // Declare table with available constrainedOptimisation methods.
    // Not used for run-time selection, just for keeping track of the
    // available methods

        declareRunTimeSelectionTable
        (
            autoPtr,
            constrainedOptimisationMethod,
            dictionary,
            (
                const fvMesh& mesh,
                const dictionary& dict,
                autoPtr<designVariables>& designVars,
                const label nConstraints,
                const word& type
            ),
            (mesh, dict, designVars, nConstraints, type)
        );


    // Constructors

        //- Construct from components
        constrainedOptimisationMethod
        (
            const fvMesh& mesh,
            const dictionary& dict,
            autoPtr<designVariables>& designVars,
            const label nConstraints,
            const word& type
        );


    //- Destructor
    virtual ~constrainedOptimisationMethod() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
