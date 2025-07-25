/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenFOAM Foundation
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
    Foam::regionModels::areaSurfaceFilmModels::BrunDrippingInjection

Description
    Film-dripping mass-transfer model.

    If the film thickness exceeds the critical value needed to generate one or
    more drops, the equivalent mass is removed from the film.  The critical film
    thickness is calculated from the Rayleigh-Taylor stability analysis of film
    flow on an inclined plane, proposed by (Brun et. al., 2015).

    Reference:
    \verbatim
        Brun, P. T., Damiano, A., Rieu, P., Balestra, G., & Gallaire, F. (2015).
        Rayleigh-Taylor instability under an inclined plane.
        Physics of Fluids, 27(8), 084107.
        DOI:10.1063/1.4927857
    \endverbatim

    The diameter of the drops formed are obtained from the local
    capillary length multiplied by the \c dCoeff coefficient.

    Reference:
    \verbatim
        Lefebvre, A. (1988).
        Atomization and sprays
        (Vol. 1040, No. 2756). CRC press.
    \endverbatim

Usage
    Minimal example:
    \verbatim
    injectionModels
    (
        BrunDrippingInjection
    );

    BrunDrippingInjectionCoeffs
    {
        // Optional entries
        ubarStar        <scalar>;
        dCoeff          <scalar>;
        deltaStable     <scalar>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                      | Type   | Reqd  | Deflt
      type      | Type name: BrunDrippingInjection | word   | yes   | -
      ubarStar  | Critical non-dimensional interface velocity <!--
                -->                                | scalar | yes   | 1.62208
      dCoeff    | Coefficient relating the diameter of the drops <!--
                --> formed to the capillary length | scalar | yes   | 3.3
      deltaStable  | Stable film thickness         | scalar | yes   | 0
    \endtable

    The inherited entries are elaborated in:
      - \link injectionModel.H \endlink

SourceFiles
    BrunDrippingInjection.C

\*---------------------------------------------------------------------------*/

#ifndef areaSurfaceFilmModels_BrunDrippingInjection_H
#define areaSurfaceFilmModels_BrunDrippingInjection_H

#include "injectionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                    Class BrunDrippingInjection Declaration
\*---------------------------------------------------------------------------*/

class BrunDrippingInjection
:
    public injectionModel
{
    // Private Member Functions

        //- No copy construct
        BrunDrippingInjection(const BrunDrippingInjection&) = delete;

        //- No copy assignment
        void operator=(const BrunDrippingInjection&) = delete;


protected:

    // Protected Data

        //- Critical non-dimensional interface velocity
        //  Coefficient in the film angle stability function.
        scalar ubarStar_;

        //- Coefficient relating the diameter of the drops formed to
        //- the capillary length.
        scalar dCoeff_;

        //- Stable film thickness - drips only formed if thickness
        //- exceeds this threshold value
        scalar deltaStable_;

        //- Diameters of particles to inject into the dripping
        scalarList diameter_;


public:

    //- Runtime type information
    TypeName("BrunDrippingInjection");


    // Constructors

        //- Construct from surface film model
        BrunDrippingInjection
        (
            liquidFilmBase& film,
            const dictionary& dict
        );


    //- Destructor
    virtual ~BrunDrippingInjection() = default;


    // Member Functions

        //- Correct
        virtual void correct
        (
            scalarField& availableMass,
            scalarField& massToInject,
            scalarField& diameterToInject
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
