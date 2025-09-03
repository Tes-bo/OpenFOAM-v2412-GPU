/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::resolutionIndexModels::CelikEtaIndex

Description
    Computes a single-mesh resolution index according to Celik et al.'s index
    using Kolmogorov length scale, which is used as a LES/DES quality/post
    verification metric that does not require any experimental or DNS data.

    \f[
        \Gamma_{Celik,\eta}(\mathbf{x}, t) =
            \frac{1}{1 + \alpha_\eta \left(\frac{h}{\eta_{eff}}\right)^m}
    \f]

    with

    \f[
        \eta_{eff} = \left(\frac{\nu^3}{\epsilon}\right)^{1/4}
    \f]

    \f[
        \epsilon = \nu_{eff} \frac{k_{sgs}}{C_k \Delta^2}
    \f]

    \f[
        \nu_{eff} = \nu_{num} + \nu_{sgs} + \nu
    \f]

    \f[
        \nu_{num} = {sgn}(k_{num}) C_\nu \Delta \sqrt{k_{num}}
    \f]

    \f[
        k_{num} = C_n \left(\frac{h}{\Delta}\right)^2 k_{sgs}
    \f]

    where
    \vartable
        \Gamma_{Celik,\eta}(\mathbf{x}, t) | Celik et al.'s index   [-]
        \alpha_\eta | Empirical constant                            [-]
        h          | Characteristic length scale with \f$h = V^{1/3} \f$ [m]
        V          | Cell volume                                    [m^3]
        \eta_{eff} | Kolmogorov length scale                        [m]
        m          | Empirical exponent                             [-]
        \nu        | Kinematic viscosity of fluid                   [m^2/s]
        \epsilon   | Kinetic energy dissipation rate                [m^2/s^3]
        \nu_{eff}  | Effective eddy viscosity                       [m^2/s]
        \nu_{num}  | Numerical eddy viscosity                       [m^2/s]
        \nu_{sgs}  | Subgrid-scale eddy viscosity                   [m^2/s]
        k_{num}    | Numerical turbulent kinetic energy             [m^2/s^2]
        C_\nu      | Empirical constant                             [-]
        \Delta     | Filter length scale                            [m]
        k_{sgs}    | Subgrid-scale turbulent kinetic energy         [m^2/s^2]
        C_n        | Empirical constant                             [-]
        C_k        | Empirical constant                             [-]
    \endvartable

    References:
    \verbatim
        Governing equations (tag:CCY):
            Celik, I. B., Cehreli Z. N., Yavuz I. (2005).
            Index of resolution quality for large eddy simulations.
            Journal of Fluids Engineering. 127:949–958.
            DOI:10.1115/1.1990201

        Governing equations (tag:CKJ):
            Celik, I., Klein, M., & Janicka, J. (2009).
            Assessment measures for engineering LES applications.
            Journal of fluids engineering, 131(3).
            DOI:10.1115/1.3059703
    \endverbatim

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    resolutionIndex1
    {
        // Inherited entries
        ...
        model       CelikEtaIndex;

        // Optional entries
        alphaEta    <scalar>;
        m           <scalar>;
        Cnu         <scalar>;
        Cn          <scalar>;
        Ck          <scalar>;
        k           <word>;
        delta       <word>;
        nu          <word>;
        nut         <word>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property    | Description                       | Type | Reqd | Deflt
      model       | Model name: CelikEtaIndex         | word | yes  | -
      alphaEta    | Empirical constant                | scalar | no | 0.05
      m           | Empirical exponent                | scalar | no | 0.5
      Cnu         | Empirical constant                | scalar | no | 0.1
      Cn          | Empirical constant                | scalar | no | 1.0
      Ck          | Empirical constant                | scalar | no | 0.0376
      k           | Name of subgrid-scale turbulent kinetic energy field <!--
                  -->                                 | word   | no | k
      delta       | Name of filter field              | word   | no | delta
      nu          | Name of kinematic viscosity field | word   | no | nu
      nut         | Name of subgrid-scale eddy viscosity field | word | no | nut
    \endtable

SourceFiles
    CelikEtaIndex.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_resolutionIndexModels_CelikEtaIndex_H
#define Foam_resolutionIndexModels_CelikEtaIndex_H

#include "resolutionIndexModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace resolutionIndexModels
{

/*---------------------------------------------------------------------------*\
                        Class CelikEtaIndex Declaration
\*---------------------------------------------------------------------------*/

class CelikEtaIndex
:
    public resolutionIndexModel
{
    // Private Data

        //- Empirical constant
        scalar alphaEta_;

        //- Empirical exponent
        scalar m_;

        //- Empirical constant
        scalar Cnu_;

        //- Empirical constant
        scalar Cn_;

        //- Empirical constant
        scalar Ck_;

        //- Name of subgrid-scale turbulent kinetic energy field
        word kName_;

        //- Name of filter field
        word deltaName_;

        //- Name of kinematic viscosity field
        word nuName_;

        //- Name of subgrid-scale eddy viscosity field
        word nutName_;


    // Private Member Functions

        //- Return Kolmogorov length scale field
        tmp<volScalarField> eta() const;

        //- Return kinetic energy dissipation rate field
        tmp<volScalarField> epsilon() const;

        //- Return effective eddy viscosity field
        tmp<volScalarField> nuEff() const;

        //- Return numerical eddy viscosity field
        tmp<volScalarField> nuNum() const;

        //- Return numerical turbulent kinetic energy field
        tmp<volScalarField> kNum() const;


public:

    //- Runtime type information
    TypeName("CelikEtaIndex");


    // Constructors

        //- Construct from components
        CelikEtaIndex
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- No copy construct
        CelikEtaIndex(const CelikEtaIndex&) = delete;

        //- No copy assignment
        void operator=(const CelikEtaIndex&) = delete;


    // Destructor
    virtual ~CelikEtaIndex() = default;


    // Member Functions

        //- Read top-level dictionary
        virtual bool read(const dictionary& dict);

        //- Calculate the result field
        virtual bool execute();

        //- Write the result field
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace resolutionIndexModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
