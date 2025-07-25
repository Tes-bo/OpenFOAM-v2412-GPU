/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::multiphaseEuler::dragModels::GidaspowSchillerNaumann

Description
    H, Enwald, E. Peirano, A-E Almstedt
    'Eulerian Two-Phase Flow Theory Applied to Fluidization'
    Int. J. Multiphase Flow, Vol. 22, Suppl, pp. 21-66 (1996)
    Eq. 86-87, p. 40

    This is identical to the Wen and Yu, Rowe model Table 3.6 p.56  in
    the Ph.D. thesis of Berend van Wachem
    'Derivation, Implementation and Validation
                    of
          Computer Simulation Models
         for Gas-Solid Fluidized Beds'

SourceFiles
    GidaspowSchillerNaumann.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_multiphaseEuler_dragModels_GidaspowSchillerNaumann_H
#define Foam_multiphaseEuler_dragModels_GidaspowSchillerNaumann_H

#include "dragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace multiphaseEuler
{
namespace dragModels
{

/*---------------------------------------------------------------------------*\
                           Class GidaspowSchillerNaumann Declaration
\*---------------------------------------------------------------------------*/

class GidaspowSchillerNaumann
:
    public dragModel
{

public:

    //- Runtime type information
    TypeName("GidaspowSchillerNaumann");


    // Constructors

        //- Construct from components
        GidaspowSchillerNaumann
        (
            const dictionary& interfaceDict,
            const phaseModel& phase1,
            const phaseModel& phase2
        );


    //- Destructor
    virtual ~GidaspowSchillerNaumann();


    // Member Functions

        tmp<volScalarField> K(const volScalarField& Ur) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dragModels
} // End namespace multiphaseEuler
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
