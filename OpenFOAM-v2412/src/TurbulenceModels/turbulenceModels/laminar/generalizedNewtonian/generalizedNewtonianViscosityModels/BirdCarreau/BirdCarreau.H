/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2020 OpenFOAM Foundation
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
    Foam::laminarModels::generalizedNewtonianViscosityModels::BirdCarreau

Description
    An incompressible Bird-Carreau non-Newtonian viscosity model.

    The Bird-Carreau-Yasuda form is also supported if the optional "a"
    coefficient is specified.  "a" defaults to 2 for the Bird-Carreau model.

    Example specification for a polymer:
    \verbatim
         viscosityModel BirdCarreau;

         nuInf       0;
         n           0.5;
    \endverbatim

SourceFiles
    BirdCarreau.C

\*---------------------------------------------------------------------------*/

#ifndef BirdCarreau_H
#define BirdCarreau_H

#include "generalizedNewtonianViscosityModel.H"
#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarModels
{
namespace generalizedNewtonianViscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class BirdCarreau Declaration
\*---------------------------------------------------------------------------*/

class BirdCarreau
:
    public generalizedNewtonianViscosityModel
{
    // Private Data

        dictionary BirdCarreauCoeffs_;

        dimensionedScalar nuInf_;
        dimensionedScalar k_;
        dimensionedScalar n_;
        dimensionedScalar a_;


public:

    //- Runtime type information
    TypeName("BirdCarreau");


    // Constructors

        //- Construct from dictionary (components)
        explicit BirdCarreau(const dictionary& viscosityProperties);


    //- Destructor
    virtual ~BirdCarreau() = default;


    // Member Functions

        //- Read transportProperties dictionary
        virtual bool read(const dictionary& viscosityProperties);

        //- Return the laminar viscosity
        virtual tmp<volScalarField> nu
        (
            const volScalarField& nu0,
            const volScalarField& strainRate
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace generalizedNewtonianViscosityModels
} // End namespace laminarModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
