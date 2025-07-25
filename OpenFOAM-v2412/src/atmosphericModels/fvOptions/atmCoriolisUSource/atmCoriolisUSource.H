/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 CENER
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
    Foam::fv::atmCoriolisUSource

Group
    grpAtmFvOptions

Description
    Applies corrections to incorporate the horizontal and vertical components
    of the Coriolis force for which the rotating frame is Earth.

    The Coriolis force is an inertial or fictitious force that acts on
    objects that are in motion within a frame of reference that rotates with
    respect to an inertial frame.

    In the atmospheric boundary layer context, for the "Coriolis effect",
    the rotating reference frame implied is almost always Earth.
    Because Earth spins, Earth-bound observers need to account for the
    Coriolis force to correctly analyze the motion of objects. Earth
    completes one rotation per day, so for motions of everyday objects the
    Coriolis force is usually quite small compared with other forces; its
    effects generally become noticeable only for motions occurring over large
    distances and long periods of time, such as large-scale movement of air in
    the atmosphere or water in the ocean. Such motions are constrained by the
    surface of Earth, so only the horizontal component of the Coriolis
    force is generally important.

    Corrections applied on:
    \verbatim
      U    | Velocity                  [m/s]
    \endverbatim

    Required fields:
    \verbatim
      U    | Velocity                  [m/s]
    \endverbatim

    References:
    \verbatim
        Coriolis force. (n.d.).
        In Wikipedia. Retrieved Feb 26, 2020, from https://w.wiki/JE5
    \endverbatim

Usage
    Example by using \c constant/fvOptions:
    \verbatim
    atmCoriolisUSource1
    {
        // Mandatory entries (unmodifiable)
        type             atmCoriolisUSource;

        atmCoriolisUSourceCoeffs
        {
            // Mandatory (inherited) entries (unmodifiable)
            selectionMode    all;

            // Conditional mandatory entries (unmodifiable)
            // Select either of the below

            // Option-1: to directly input rotation vector
            Omega        (0 0 5.65156e-5);

            // Option-2: to indirectly input rotation vector
            // by a latitude-period pair
            latitude                   51.971;
            planetaryRotationPeriod    23.9344694;
        }

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
        Property    | Description                        | Type | Req'd  | Dflt
        type        | Type name: atmCoriolisUSource      | word | yes    | -
        latitude    | Geographic coordinate specifying the north–south  <!--
                  --> position of a point on the surface of a planetary <!--
                  --> body [degree]         | scalar | conditional  | 0.0
        planetaryRotationPeriod | Rotation period of the planetary body <!--
                              --> [h]       | scalar | conditional  | 23.9344694
        Omega    | Rotation vector of the rotating reference frame      <!--
               --> relative to the inertial frame [rad/s]               <!--
               -->                          | vector | conditional  | (0 0 0)
    \endtable

    The inherited entries are elaborated in:
      - \link fvOption.H \endlink
      - \link cellSetOption.H \endlink

Note
    - Dimensional consistencies are hard-coded; therefore, no internal
    check is performed for potential dimension inconsistencies.
    - The magnitude of the \c latitude is limited to [0, 90], yet its value
    is allowed to be negative for the southern hemisphere.
    - The Coriolis force for cell whose volume is less than \c VSMALL is
    equated to zero.

SourceFiles
    atmCoriolisUSource.C

\*---------------------------------------------------------------------------*/

#ifndef atmCoriolisUSource_H
#define atmCoriolisUSource_H

#include "cellSetOption.H"
#include "dimensionedVector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                  Class atmCoriolisUSource Declaration
\*---------------------------------------------------------------------------*/

class atmCoriolisUSource
:
    public fv::cellSetOption
{
    // Private Data

        //- Latitude on the planetary body
        const scalar latitude_;

        //- Rotation period of the planetary body
        const scalar planetaryRotationPeriod_;

        //- Planetary rotation vector
        const dimensionedVector Omega_;


    // Private Member Functions

        //- Rotation vector of the planetary body
        vector planetaryRotationVector() const;


public:

    //- Runtime type information
    TypeName("atmCoriolisUSource");


    // Constructors

        //- Construct from explicit source name and mesh
        atmCoriolisUSource
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        atmCoriolisUSource(const atmCoriolisUSource&) = delete;

        //- No copy assignment
        void operator=(const atmCoriolisUSource&) = delete;


    //- Destructor
    virtual ~atmCoriolisUSource() = default;


    // Member Functions

        //- Add explicit contribution to incompressible momentum equation
        virtual void addSup
        (
            fvMatrix<vector>& eqn,
            const label fieldi
        );

        //- Add explicit contribution to compressible momentum equation
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<vector>& eqn,
            const label fieldi
        );

        //- Add explicit contribution to phase momentum equation
        virtual void addSup
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            fvMatrix<vector>& eqn,
            const label fieldi
        );

        //- Read source dictionary
        virtual bool read(const dictionary& dict)
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
