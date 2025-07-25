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
    Foam::decompositionConstraints::singleProcessorFaceSets

Description
    Constraint to keep all cells connected to face or point of faceSet on
    a single processor.

    \heading Dictionary parameters
    \table
        Property    | Description                       | Required  | Default
        type        | singleProcessorFaceSets           | yes   |
        sets        | List of sets and processers       | yes   |
    \endtable

    Example,
    \verbatim
    constraints
    {
        procs
        {
            type    singleProcessorFaceSets;
            sets
            (
                (faceSet1   0)
                (faceSet5   12)
                (faceSet25 510)
            );
        }
    }
    \endverbatim

SourceFiles
    singleProcessorFaceSetsConstraint.C

\*---------------------------------------------------------------------------*/

#ifndef singleProcessorFaceSetsConstraint_H
#define singleProcessorFaceSetsConstraint_H

#include "decompositionConstraint.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace decompositionConstraints
{

/*---------------------------------------------------------------------------*\
                   Class singleProcessorFaceSets Declaration
\*---------------------------------------------------------------------------*/

class singleProcessorFaceSets
:
    public decompositionConstraint
{
    // Private data

        //- List of faceSet+processor
        List<Tuple2<word, label>> setNameAndProcs_;

        //- Report names of sets and processor numnbers
        void printInfo() const;


public:

    //- Runtime type information
    TypeName("singleProcessorFaceSets");


    // Constructors

        //- Construct with constraint dictionary
        explicit singleProcessorFaceSets(const dictionary& dict);

        //- Copy construct from components
        explicit singleProcessorFaceSets
        (
            const List<Tuple2<word, label>>& setNameAndProcs
        );

        //- Read construct a List of (setName processor)
        explicit singleProcessorFaceSets(Istream& is);


    //- Destructor
    virtual ~singleProcessorFaceSets() = default;


    // Member Functions

        //- Add this constraint to list of constraints
        virtual void add
        (
            const polyMesh& mesh,
            boolList& blockedFace,
            PtrList<labelList>& specifiedProcessorFaces,
            labelList& specifiedProcessor,
            List<labelPair>& explicitConnections
        ) const;

        //- Add this constraint post-decomposition
        virtual void apply
        (
            const polyMesh& mesh,
            const boolList& blockedFace,
            const PtrList<labelList>& specifiedProcessorFaces,
            const labelList& specifiedProcessor,
            const List<labelPair>& explicitConnections,
            labelList& decomposition
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace decompositionConstraints
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
