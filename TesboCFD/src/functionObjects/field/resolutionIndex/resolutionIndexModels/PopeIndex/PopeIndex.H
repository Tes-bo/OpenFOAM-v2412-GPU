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
    Foam::resolutionIndexModels::PopeIndex

Description
    Computes a single-mesh resolution index according to Pope's index,
    which is used as a LES/DES quality/post-verification metric that does
    not require any experimental or DNS data.

    \f[
        \Gamma_{Pope}(\mathbf{x}, t) = \frac{k_{res}}{k_{tot}}
    \f]

    with

    \f[
        k_{tot} = k_{res} + k_{sgs} + |k_{num}|
    \f]

    where
    \vartable
        \Gamma_{Pope}(\mathbf{x}, t) | Pope's index                 [-]
        k_{tot} | Total turbulent kinetic energy                    [m^2/s^2]
        k_{res} | Resolved turbulent kinetic energy                 [m^2/s^2]
        k_{sgs} | Subgrid-scale turbulent kinetic energy            [m^2/s^2]
        k_{num} | Numerical turbulent kinetic energy                [m^2/s^2]
    \endvartable

    Inclusion of \f$k_{num}\f$ is optional, and set as \c true by default:

    \f[
        k_{num} = C_n \left(\frac{h}{\Delta}\right)^2 k_{sgs}
    \f]

    where
    \vartable
        C_n     | Empirical constant                                  [-]
        h       | Characteristic length scale with \f$h = V^{1/3} \f$ [m]
        V       | Cell volume                                         [m^3]
        \Delta  | Filter length scale                                 [m]
    \endvartable

    Typical criterion for acceptable-quality resolution:

    \f[
        \Gamma_{Pope}(\mathbf{x}) \geq 0.8
    \f]

    Required fields:
    \verbatim
      U        | Velocity                                   [m/s]
      UMean    | Mean velocity                              [m/s]
      k        | Subgrid-scale turbulent kinetic energy     [m^2/s^2]
      delta    | Filter length                              [m]
    \endverbatim

    References:
    \verbatim
        Governing equation (tag:P):
            Pope, S. B. (2000).
            Turbulent flows.
            Cambridge, UK: Cambridge Univ. Press
            DOI:10.1017/CBO9780511840531

        Governing equations for the denominator kNum term (tag:CKJ):
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
        model       PopeIndex;

        // Optional entries
        U           <word>;
        UMean       <word>;
        k           <word>;
        delta       <word>;
        includeKnum <bool>;

        // Conditional entries
        // when includeKnum = true
        Cn          <scalar>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property    | Description                       | Type | Reqd | Deflt
      model       | Model name: PopeIndex             | word | yes  | -
      U           | Name of velocity field            | word | no   | U
      UMean       | Name of mean velocity field       | word | no   | UMean
      k           | Name of subgrid-scale turbulent kinetic energy field <!--
                  -->                                 | word   | no | k
      delta       | Name of filter field              | word   | no | delta
      includeKnum | Flag to include numerical k field | bool   | no | true
      Cn          | Empirical constant                | choice | no | 1.0
    \endtable

Note
  - Some studies measured \f$\Gamma_{Pope} > 1\f$ with \f$k_{res}\f$ comparisons
    between a LES and a corresponding filtered DNS. Nonphysical results may
    occur.

SourceFiles
    PopeIndex.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_resolutionIndexModels_PopeIndex_H
#define Foam_resolutionIndexModels_PopeIndex_H

#include "resolutionIndexModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace resolutionIndexModels
{

/*---------------------------------------------------------------------------*\
                          Class PopeIndex Declaration
\*---------------------------------------------------------------------------*/

class PopeIndex
:
    public resolutionIndexModel
{
    // Private Data

        //- Flag to include numerical turbulent kinetic energy field
        bool includeKnum_;

        //- Empirical constant in numerical turbulent kinetic energy estimation
        scalar Cn_;

        //- Name of velocity field
        word UName_;

        //- Name of mean velocity field
        word UMeanName_;

        //- Name of subgrid-scale turbulent kinetic energy field
        word kName_;

        //- Name of filter field
        word deltaName_;


    // Private Member Functions

        //- Return numerical turbulent kinetic energy field
        tmp<volScalarField> kNum() const;


public:

    //- Runtime type information
    TypeName("PopeIndex");


    // Constructors

        //- Construct from components
        PopeIndex
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- No copy construct
        PopeIndex(const PopeIndex&) = delete;

        //- No copy assignment
        void operator=(const PopeIndex&) = delete;


    // Destructor
    virtual ~PopeIndex() = default;


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
