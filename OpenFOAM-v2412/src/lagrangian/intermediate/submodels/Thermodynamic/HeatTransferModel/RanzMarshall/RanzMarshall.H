/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::RanzMarshall

Group
    grpLagrangianIntermediateHeatTransferSubModels

Description
    Nusselt-number model using the empirical Ranz-Marshall correlation
    to be used in modelling of the fluid-particle heat transfer coefficient:

    \f[
        \mathrm{Nu} = a + b \, \mathrm{Re}_p^{m} \, \mathrm{Pr}^{n}
    \f]
    with

    \f[
        \mathrm{Re}_p =
            \frac{\rho_c \, | \mathbf{u}_\mathrm{rel} | \, d_p}{\mu_c}
    \f]

    \f[
        \mathrm{Pr} = \frac{ C_p \, \mu_c }{ \kappa_c }
    \f]
    where

    \vartable
        \mathrm{Nu}    | Nusselt number
        \mathrm{Re}_p  | Particle Reynolds number
        \mathrm{Pr}    | Prandtl number
        d_p            | Particle diameter
        \rho_c         | Density of carrier in the film surrounding particle
        \mu_c    | Dynamic viscosity of carrier in the film surrounding particle
        \mathbf{u}_\mathrm{rel} | Relative velocity between particle and carrier
        a              | Correlation coefficient
        b              | Correlation coefficient
        m              | Correlation exponent of particle Reynolds number
        n              | Correlation exponent of Prandtl number
        C_p            | Specific heat capacity
        \kappa_c       | Thermal conductivity of carrier in the film
    \endvartable

    Reference:
    \verbatim
        Standard model:
            Ranz, W. E., & Marshall, W. R. (1952).
            Evaporation from drops - part 1.
            Chem. Eng. Prog, 48, 22, pp. 141-146.

            Ranz, W. E., & Marshall, W. R. (1952).
            Evaporation from drops - part 2.
            Chem. Eng. Prog, 48, 4, pp. 173-180.

        Expressions (tag:AOB), p. 18:
            Amsden, A. A., O'Rourke, P. J., & Butler, T. D. (1989).
            KIVA-II: A computer program for chemically
            reactive flows with sprays (No. LA-11560-MS).
            Los Alamos National Lab.(LANL), Los Alamos, NM (United States).
            DOI:10.2172/6228444
    \endverbatim

Usage
    Minimal example by using \c constant/\<CloudProperties\>:
    \verbatim
    subModels
    {
        // Mandatory entries
        heatTransferModel    RanzMarshall;

        // Optional entries
        RanzMarshallCoeffs
        {
            a    2.0;
            b    0.6;
            m    0.5;
            n    0.66666;
        }
    }
    \endverbatim

    where the entries mean:
    \table
      Property           | Description              | Type | Reqd | Deflt
      heatTransferModel  | Type name: RanzMarshall  | word | yes  | -
      a                  | Correlation coefficient  | scalar | no | 2.0
      b                  | Correlation coefficient  | scalar | no | 0.6
      m  | Correlation exponent of particle Reynolds number | scalar | no | 0.5
      n  | Correlation exponent of Prandtl number   | scalar | no | 1.0/3.0
    \endtable

SourceFiles
    RanzMarshall.C

\*---------------------------------------------------------------------------*/

#ifndef RanzMarshall_H
#define RanzMarshall_H

#include "HeatTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class RanzMarshall Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class RanzMarshall
:
    public HeatTransferModel<CloudType>
{
    // Private Data

        //- Correlation coefficient
        const scalar a_;

        //- Correlation coefficient
        const scalar b_;

        //- Correlation exponent of particle Reynolds number
        const scalar m_;

        //- Correlation exponent of Prandtl number
        const scalar n_;


public:

    //- Runtime type information
    TypeName("RanzMarshall");


    // Generated Methods

        //- No copy assignment
        void operator=(const RanzMarshall&) = delete;


    // Constructors

        //- Construct from dictionary
        RanzMarshall(const dictionary& dict, CloudType& cloud);

        //- Copy construct
        RanzMarshall(const RanzMarshall<CloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<HeatTransferModel<CloudType>> clone() const
        {
            return autoPtr<HeatTransferModel<CloudType>>
            (
                new RanzMarshall<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~RanzMarshall() = default;


    // Member Functions

        // Evaluation

            //- Return Nusselt number
            virtual scalar Nu
            (
                const scalar Re,
                const scalar Pr
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "RanzMarshall.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
