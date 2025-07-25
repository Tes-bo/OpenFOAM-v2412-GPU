/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2023 PCOpt/NTUA
    Copyright (C) 2020-2023 FOSS GP
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
    Foam::betaMax

Description
    Base class for selecting the betaMax value, i.e. the value multiplying the
    Brinkman penalisation term, for topology optimisation.

SourceFiles
    betaMax.C

\*---------------------------------------------------------------------------*/

#ifndef betaMax_H
#define betaMax_H

#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class betaMax Declaration
\*---------------------------------------------------------------------------*/

class betaMax
{
private:

    // Private Member Functions

        //- No copy construct
        betaMax(const betaMax&) = delete;

        //- No copy assignment
        void operator=(const betaMax&) = delete;


protected:

    // Protected Data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- betaMax value
        scalar value_;


    // Protected Member Functions

        //- Compute the characteristic length
        scalar computeLength(const dictionary& dict) const;

        //- Mark all common inlet - processor edges
        boolList markProcessorEdges(const polyPatch& patch) const;


public:

    //- Runtime type information
    TypeName("betaMax");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            betaMax,
            dictionary,
            (
                const fvMesh& mesh,
                const dictionary& dict
            ),
            (mesh, dict)
        );


    // Constructors

        //- Construct from components
        betaMax
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Selectors

        //- Construct and return the selected betaMax model
        static autoPtr<betaMax> New
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~betaMax() = default;


    // Member Functions

        //- Get value
        virtual scalar value() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
