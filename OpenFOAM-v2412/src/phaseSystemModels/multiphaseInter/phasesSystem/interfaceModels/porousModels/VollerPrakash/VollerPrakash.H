/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::porousModels::VollerPrakash

Description
    Porous model to apply a porous media on a solid phase.

    Based on the references:
    -# V.R. Voller and C. Prakash, A fixed grid numerical modelling
       methodology for convection-diffusion mushy phase-change problems,
       Int. J. Heat Mass Transfer 30(8):17091719, 1987.
    -# C.R. Swaminathan. and V.R. Voller, A general enthalpy model for
       modeling solidification processes, Metallurgical Transactions
       23B:651664, 1992.
Usage
    Example usage:
    \verbatim
    interfacePorous
    (
        (solid and liquid)
        {
            type            VollerPrakash;
            solidPhase      alpha.solid;
            Cu              1e7;
        }
    );
    \endverbatim

    Where:

    \table
        Property    | Description           | Required    | Default value
        solidPhase  | Solid phase name      | yes
        Cu          | Model constant        | yes
    \endtable

SourceFiles
    VollerPrakash.C

\*---------------------------------------------------------------------------*/

#ifndef porousModels_VollerPrakash_H
#define porousModels_VollerPrakash_H

#include "porousModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace porousModels
{

/*---------------------------------------------------------------------------*\
              Class VollerPrakash Declaration
\*---------------------------------------------------------------------------*/

class VollerPrakash
:
    public porousModel
{
    // Private data

        //- Constant surface tension value
        const scalar Cu_;

        //- Solid phase name
        word solidPhase_;


public:

    //- Runtime type information
    TypeName("VollerPrakash");


    // Constructors

        //- Construct from a dictionary and a phase pair
        VollerPrakash
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~VollerPrakash() = default;


    // Member Functions

        //- Aspect ratio
        virtual tmp<volScalarField> S() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace porousModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
