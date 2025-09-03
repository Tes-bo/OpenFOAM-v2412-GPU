/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2020 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::laminarModels::generalizedNewtonianViscosityModels::strainRateFunction

Description
    Run-time selected strain-rate function generalized Newtonian viscosity model

    Example linear function of strain-rate:
    \verbatim
        generalizedNewtonianModel strainRateFunction;

        function polynomial ((0 0.1) (1 1.3));
    \endverbatim

See also
    Foam::generalizedNewtonianViscosityModel
    Foam::Function1

SourceFiles
    strainRateFunction.C

\*---------------------------------------------------------------------------*/

#ifndef strainRateFunction_H
#define strainRateFunction_H

#include "generalizedNewtonianViscosityModel.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarModels
{
namespace generalizedNewtonianViscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class strainRateFunction Declaration
\*---------------------------------------------------------------------------*/

class strainRateFunction
:
    public generalizedNewtonianViscosityModel
{
    // Private Data

        //- Coefficients dictionary
        dictionary strainRateFunctionCoeffs_;

        //- Strain-rate function
        autoPtr<Function1<scalar>> strainRateFunction_;


public:

    //- Runtime type information
    TypeName("strainRateFunction");


    // Constructors

        //- Construct from dictionary (components)
        explicit strainRateFunction(const dictionary& viscosityProperties);


    //- Destructor
    virtual ~strainRateFunction() = default;


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
