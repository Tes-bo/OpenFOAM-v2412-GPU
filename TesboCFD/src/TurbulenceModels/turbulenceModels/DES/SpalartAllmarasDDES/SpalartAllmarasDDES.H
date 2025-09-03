/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2022, 2024 Upstream CFD GmbH
    Copyright (C) 2019-2024 OpenCFD Ltd.
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
    Foam::LESModels::SpalartAllmarasDDES

Group
    grpDESTurbulence

Description
    SpalartAllmaras DDES turbulence model for incompressible and compressible
    flows.

    Reference:
    \verbatim
        Spalart, P. R., Deck, S., Shur, M. L., Squires,
        K. D., Strelets, M. K., & Travin, A. (2006).
        A new version of detached-eddy simulation,
        resistant to ambiguous grid densities.
        Theoretical and computational fluid dynamics, 20(3), 181-195.
        DOI:10.1007/s00162-006-0015-0
    \endverbatim

    Reference for enhanced shielding function formulation:
    \verbatim
        Deck, S., Renard, N. (2020).
        Towards an enhanced protection of attached boundary layers in hybrid
        RANS/LES methods.
        Journal of Computational Physics, 400, 108970.
        DOI:10.1016/j.jcp.2019.108970
    \endverbatim

SourceFiles
    SpalartAllmarasDDES.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_SpalartAllmarasDDES_H
#define Foam_SpalartAllmarasDDES_H

#include "SpalartAllmarasBase.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                     Class SpalartAllmarasDDES Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class SpalartAllmarasDDES
:
    public SpalartAllmarasDES<BasicTurbulenceModel>
{
public:

    // Public enumerations

        //- Shielding modes
        enum class shieldingMode
        {
            standard,
            ZDES2020
        };

        //- Shielding mode names
        static const Enum<shieldingMode> shieldingModeNames;


private:

    // Private Member Functions

        //- Return the shielding function
        tmp<volScalarField> fd(const volScalarField& magGradU) const;

        //- No copy construct
        SpalartAllmarasDDES(const SpalartAllmarasDDES&) = delete;

        //- No copy assignment
        void operator=(const SpalartAllmarasDDES&) = delete;


protected:

    // Protected Data

        //- Shielding mode
        shieldingMode shielding_;

        // Model coefficients

            dimensionedScalar Cd1_;
            dimensionedScalar Cd2_;
            dimensionedScalar Cd3_;
            dimensionedScalar Cd4_;
            dimensionedScalar betaZDES_;
            Switch usefP2_;


    // Protected Member Functions

        //- Return the production term
        virtual tmp<volScalarField> Stilda
        (
            const volScalarField& chi,
            const volScalarField& fv1,
            const volTensorField& gradU,
            const volScalarField& dTilda
        ) const;

        //- Return the length scale
        virtual tmp<volScalarField> dTilda
        (
            const volScalarField& chi,
            const volScalarField& fv1,
            const volTensorField& gradU
        ) const;


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("SpalartAllmarasDDES");


    // Constructors

        //- Construct from components
        SpalartAllmarasDDES
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~SpalartAllmarasDDES() = default;


    // Member Functions

        //- Read from dictionary
        virtual bool read();

        //- Return the shielding function
        virtual tmp<volScalarField> fd() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SpalartAllmarasDDES.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
