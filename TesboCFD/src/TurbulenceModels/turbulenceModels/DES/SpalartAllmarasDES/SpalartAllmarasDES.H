/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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
    Foam::LESModels::SpalartAllmarasDES

Group
    grpDESTurbulence

Description
    SpalartAllmarasDES DES turbulence model for incompressible and
    compressible flows.

    Reference:
    \verbatim
        Spalart, P. R., Jou, W. H., Strelets, M., & Allmaras, S. R. (1997).
        Comments on the feasibility of LES for wings, and on a hybrid
        RANS/LES approach.
        Advances in DNS/LES, 1, 4-8.
    \endverbatim

    Including the low Reynolds number correction documented in
    \verbatim
        Spalart, P. R., Deck, S., Shur, M. L., Squires,
        K. D., Strelets, M. K., & Travin, A. (2006).
        A new version of detached-eddy simulation,
        resistant to ambiguous grid densities.
        Theoretical and computational fluid dynamics, 20(3), 181-195.
        DOI:10.1007/s00162-006-0015-0
    \endverbatim

Note
    The default value of the DES constant implemented was calibrated for
    OpenFOAM using decaying isotropic turbulence and agrees with the value
    suggested in the reference publication.

SourceFiles
    SpalartAllmarasDES.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_SpalartAllmarasDES_H
#define Foam_SpalartAllmarasDES_H

#include "SpalartAllmarasBase.H"
#include "DESModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                        Class SpalartAllmarasDES Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class SpalartAllmarasDES
:
    public SpalartAllmarasBase<DESModel<BasicTurbulenceModel>>
{
    // Private Member Functions

        //- No copy construct
        SpalartAllmarasDES(const SpalartAllmarasDES&) = delete;

        //- No copy assignment
        void operator=(const SpalartAllmarasDES&) = delete;


protected:

    // Protected Data

        //- Switch to activate grey-area enhanced sigma-(D)DES
        Switch useSigma_;

        // Model constants

            //- DES coefficient
            dimensionedScalar CDES_;

            //- Flag for low Reynolds number correction
            Switch lowReCorrection_;
            dimensionedScalar fwStar_;


    // Protected Member Functions

        //- Return the low Reynolds number correction function
        virtual tmp<volScalarField> psi
        (
            const volScalarField& chi,
            const volScalarField& fv1
        ) const;

        //- Return the LES length scale
        virtual tmp<volScalarField> lengthScaleLES
        (
            const volScalarField& chi,
            const volScalarField& fv1
        ) const;

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

        virtual void correctNut();

public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("SpalartAllmarasDES");


    // Constructors

        //- Construct from components
        SpalartAllmarasDES
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
    virtual ~SpalartAllmarasDES() = default;


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the LES field indicator
        virtual tmp<volScalarField> LESRegion() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SpalartAllmarasDES.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
