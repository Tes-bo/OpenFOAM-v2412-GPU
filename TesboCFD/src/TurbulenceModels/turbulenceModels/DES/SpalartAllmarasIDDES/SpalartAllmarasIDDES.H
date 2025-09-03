/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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
    Foam::LESModels::SpalartAllmarasIDDES

Group
    grpDESTurbulence

Description
    SpalartAllmaras IDDES turbulence model
    for incompressible and compressible flows.

    Reference:
    \verbatim
        Shur, M. L., Spalart, P. R., Strelets, M. K., & Travin, A. K. (2008).
        A hybrid RANS-LES approach with delayed-DES
        and wall-modelled LES capabilities.
        International journal of heat and fluid flow, 29(6), 1638-1649.
        DOI:10.1016/j.ijheatfluidflow.2008.07.001
    \endverbatim

SourceFiles
    SpalartAllmarasIDDES.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_SpalartAllmarasIDDES_H
#define Foam_SpalartAllmarasIDDES_H

#include "SpalartAllmarasDES.H"
#include "IDDESDelta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                    Class SpalartAllmarasIDDES Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class SpalartAllmarasIDDES
:
    public SpalartAllmarasDES<BasicTurbulenceModel>
{
    // Private Member Functions

        //- Check that the supplied delta is an IDDESDelta
        const IDDESDelta& setDelta() const;

        tmp<volScalarField> alpha() const;

        tmp<volScalarField> ft(const volScalarField& magGradU) const;

        tmp<volScalarField> fl(const volScalarField& magGradU) const;

        //- Delay function
        tmp<volScalarField> fdt(const volScalarField& magGradU) const;

        //- No copy construct
        SpalartAllmarasIDDES(const SpalartAllmarasIDDES&) = delete;

        //- No copy assignment
        void operator=(const SpalartAllmarasIDDES&) = delete;


protected:

    // Protected Data

        // Model coefficients

            dimensionedScalar Cdt1_;
            dimensionedScalar Cdt2_;
            dimensionedScalar Cl_;
            dimensionedScalar Ct_;
            Switch fe_;


        //- IDDES delta
        const IDDESDelta& IDDESDelta_;


    // Protected Member Functions

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
    TypeName("SpalartAllmarasIDDES");


    // Constructors

        //- Construct from components
        SpalartAllmarasIDDES
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
    virtual ~SpalartAllmarasIDDES() = default;


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the shielding function
        virtual tmp<volScalarField> fd() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SpalartAllmarasIDDES.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
