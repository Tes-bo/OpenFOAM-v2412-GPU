/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
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
    Foam::ATCUaGradU

Description
    The ATC formualtion resulting by differentiating the
    Conservative form of the Momentum equations.

SourceFiles
    ATCUaGradU.C

\*---------------------------------------------------------------------------*/

#ifndef ATCUaGradU_H
#define ATCUaGradU_H

#include "ATCModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class ATCUaGradU Declaration
\*---------------------------------------------------------------------------*/

class ATCUaGradU
:
    public ATCModel
{
private:

    // Private Member Functions

        //- No copy construct
        ATCUaGradU(const ATCUaGradU&) = delete;

        //- No copy assignment
        void operator=(const ATCUaGradU&) = delete;


public:

    //- Runtime type information
    TypeName("UaGradU");


    // Constructors

        //- Construct from components
        ATCUaGradU
        (
            const fvMesh& mesh,
            const incompressibleVars& primalVars,
            const incompressibleAdjointVars& adjointVars,
            const dictionary& dict
        );


    //- Destructor
    virtual ~ATCUaGradU() = default;


    // Member Functions

        //- Add ATC
        virtual void addATC(fvVectorMatrix& UaEqn);

        //- Get the FI sensitivity derivatives term coming from the ATC
        virtual tmp<volTensorField> getFISensitivityTerm() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
