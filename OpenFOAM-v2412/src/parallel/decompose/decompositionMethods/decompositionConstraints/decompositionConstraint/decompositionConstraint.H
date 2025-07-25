/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2016 OpenFOAM Foundation
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
    Foam::decompositionConstraint

Description
    Abstract class for handling decomposition constraints.

SourceFiles
    decompositionConstraint.C

\*---------------------------------------------------------------------------*/

#ifndef decompositionConstraint_H
#define decompositionConstraint_H

#include "dictionary.H"
#include "boolList.H"
#include "labelList.H"
#include "labelPair.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class polyMesh;

/*---------------------------------------------------------------------------*\
                   Class decompositionConstraint Declaration
\*---------------------------------------------------------------------------*/

class decompositionConstraint
{
protected:

    // Protected data

        //- Constraint coefficients dictionary
        dictionary coeffDict_;


    // Protected Member Functions

        //- Get minimum label across coupled boundary faces
        void getMinBoundaryValue
        (
            const polyMesh& mesh,
            const labelList& decomposition,
            labelList& destProc
        ) const;


        //- No copy construct
        decompositionConstraint(const decompositionConstraint&) = delete;

        //- No copy assignment
        void operator=(const decompositionConstraint&) = delete;

public:

    //- Runtime type information
    TypeName("decompositionConstraint");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            decompositionConstraint,
            dictionary,
            (
                const dictionary& dict
            ),
            (dict)
        );


    // Constructors

        //- Construct with constraint dictionary
        explicit decompositionConstraint(const dictionary& constraintDict);

        //- Construct with constraint dictionary and model type.
        //  The model type could be used for defining a coefficients
        //  sub-dictionary.
        decompositionConstraint
        (
            const dictionary& constraintDict,
            const word& modelType
        );


    // Selectors

        //- Return a reference to the selected decompositionConstraint
        static autoPtr<decompositionConstraint> New
        (
            const dictionary& constraintDict
        );

        //- Return a reference to the selected decompositionConstraint
        static autoPtr<decompositionConstraint> New
        (
            const dictionary& constraintDict,
            const word& modelType
        );


    //- Destructor
    virtual ~decompositionConstraint() = default;


    // Member Functions

        //- Add this constraint to list of constraints
        virtual void add
        (
            const polyMesh& mesh,
            boolList& blockedFace,
            PtrList<labelList>& specifiedProcessorFaces,
            labelList& specifiedProcessor,
            List<labelPair>& explicitConnections
        ) const = 0;

        //- Apply any additional post-decomposition constraints.
        //  Usually no need to do anything since decomposition method
        //  should have already obeyed the constraints
        virtual void apply
        (
            const polyMesh& mesh,
            const boolList& blockedFace,
            const PtrList<labelList>& specifiedProcessorFaces,
            const labelList& specifiedProcessor,
            const List<labelPair>& explicitConnections,
            labelList& decomposition
        ) const
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
