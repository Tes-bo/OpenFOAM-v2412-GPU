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
    Foam::ATCstandard

Description
    The ATC formualtion resulting by differentiating the
    Non-conservative form of the momentum equations.

SourceFiles
    ATCstandard.C

\*---------------------------------------------------------------------------*/

#ifndef ATCstandard_H
#define ATCstandard_H

#include "ATCModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class ATCstandard Declaration
\*---------------------------------------------------------------------------*/

class ATCstandard
:
    public ATCModel
{
private:

    // Private data

        //- The gradU used in the computation of the standard ATC
        //  Cached to avoid costly recomputation in each adjoint iteration
        volTensorField gradU_;

    // Private Member Functions

        //- No copy construct
        ATCstandard(const ATCstandard&) = delete;

        //- No copy assignment
        void operator=(const ATCstandard&) = delete;


public:

    //- Runtime type information
    TypeName("standard");


    // Constructors

        //- Construct from components
        ATCstandard
        (
            const fvMesh& mesh,
            const incompressibleVars& primalVars,
            const incompressibleAdjointVars& adjointVars,
            const dictionary& dict
        );


    //- Destructor
    virtual ~ATCstandard() = default;


    // Member Functions

        //- Add ATC
        virtual void addATC(fvVectorMatrix& UaEqn);

        //- Get the FI sensitivity derivatives term coming from the ATC
        virtual tmp<volTensorField> getFISensitivityTerm() const;

        //- Update quantities related with the primal fields
        virtual void updatePrimalBasedQuantities();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
