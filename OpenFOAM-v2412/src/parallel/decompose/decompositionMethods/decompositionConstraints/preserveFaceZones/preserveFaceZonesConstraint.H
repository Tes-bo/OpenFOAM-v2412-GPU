/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2016 OpenFOAM Foundation
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
    Foam::decompositionConstraints::preserveFaceZones

Description
    Constraint to keep/move owner and neighbour of faceZone onto same
    processor.

    \heading Dictionary parameters
    \table
        Property    | Description                       | Required  | Default
        type        | preserveFaceZones                 | yes   |
        zones       | List of face zone names or regexs | yes   |
    \endtable

SourceFiles
    preserveFaceZonesConstraint.C

\*---------------------------------------------------------------------------*/

#ifndef preserveFaceZonesConstraint_H
#define preserveFaceZonesConstraint_H

#include "decompositionConstraint.H"
#include "wordRes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace decompositionConstraints
{

/*---------------------------------------------------------------------------*\
                      Class preserveFaceZones Declaration
\*---------------------------------------------------------------------------*/

class preserveFaceZones
:
    public decompositionConstraint
{
    // Private data

        //- List of zones to keep together
        wordRes zones_;


public:

    //- Runtime type information
    TypeName("preserveFaceZones");


    // Constructors

        //- Construct with constraint dictionary
        explicit preserveFaceZones(const dictionary& dict);

        //- Construct from components
        explicit preserveFaceZones(const UList<wordRe>& zones);


    //- Destructor
    virtual ~preserveFaceZones() = default;


    // Member Functions

        //- Add my constraints to list of constraints
        virtual void add
        (
            const polyMesh& mesh,
            boolList& blockedFace,
            PtrList<labelList>& specifiedProcessorFaces,
            labelList& specifiedProcessor,
            List<labelPair>& explicitConnections
        ) const;

        //- Apply any additional post-decomposition constraints
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
