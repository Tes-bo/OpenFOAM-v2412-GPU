/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenCFD Ltd.
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
    Foam::thermoIncompressibleTwoPhaseMixture

Description
    A two-phase incompressible transportModel

SourceFiles
    thermoIncompressibleTwoPhaseMixture.C

\*---------------------------------------------------------------------------*/

#ifndef thermoIncompressibleTwoPhaseMixture_H
#define thermoIncompressibleTwoPhaseMixture_H

#include "incompressibleTwoPhaseMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class thermoIncompressibleTwoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class thermoIncompressibleTwoPhaseMixture
:
    public incompressibleTwoPhaseMixture
{
protected:

    // Protected data


        //- Thermal variables
        dimensionedScalar kappa1_;
        dimensionedScalar kappa2_;

        dimensionedScalar Cp1_;
        dimensionedScalar Cp2_;

        dimensionedScalar Cv1_;
        dimensionedScalar Cv2_;

        dimensionedScalar Hf1_;
        dimensionedScalar Hf2_;


public:

    TypeName("thermoIncompressibleTwoPhaseMixture");


    // Constructors

        //- Construct from U and phi
        thermoIncompressibleTwoPhaseMixture
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~thermoIncompressibleTwoPhaseMixture() = default;


    // Access function

        //- Return const-access to phase2 kappa
        const dimensionedScalar& kappa2() const
        {
            return kappa2_;
        };

        //- Return const-access to phase1 kappa
        const dimensionedScalar& kappa1() const
        {
            return kappa1_;
        };

        //- Return const-access to phase1 Cp
        const dimensionedScalar& Cp1() const
        {
            return Cp1_;
        };

        //- Return const-access to phase1 Cp
        const dimensionedScalar& Cp2() const
        {
            return Cp2_;
        };

        //- Return const-access to phase1 Cv
        const dimensionedScalar& Cv1() const
        {
            return Cv1_;
        };

        //- Return const-access to phase1 Cv
        const dimensionedScalar& Cv2() const
        {
            return Cv2_;
        };

        //- Return latent heat for phase 1
        const dimensionedScalar& Hf1() const
        {
            return Hf1_;
        };

        //- Return latent heat for phase 2
        const dimensionedScalar& Hf2() const
        {
            return Hf2_;
        };

        //- Read base transportProperties dictionary
        virtual bool read();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
