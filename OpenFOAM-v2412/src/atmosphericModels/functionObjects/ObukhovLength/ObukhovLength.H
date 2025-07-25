/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 ENERCON GmbH
    Copyright (C) 2020 OpenCFD Ltd
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
    Foam::functionObjects::ObukhovLength

Group
    grpFieldFunctionObjects

Description
    Computes the Obukhov length field and associated friction velocity field.

    When scaled by the ground-normal height, i.e. \c z, the Obukhov length
    becomes a dimensionless stability parameter, i.e. \c z/L, for atmospheric
    boundary layer modelling, expressing the relative roles of buoyancy and
    shear in the production and dissipation of turbulent kinetic energy.

    \f[
        u^* = \sqrt{\max (\nu_t \sqrt{2 |\symm{\grad{u}}|^2}, VSMALL)}
    \f]

    \f[
        L_o = - \frac{(u^*)^3}{sign(B) \kappa \max (|B|, VSMALL)}
    \f]

    with
    \f[
        B = \alpha_t \beta \frac{\grad{T} \cdot \vec{g}}{\rho}
    \f]

    where
    \vartable
      u^*        | Friction velocity                        [m/s]
      \nu_t      | Turbulent viscosity                      [m2/s]
      u          | Velocity                                 [m/s]
      L_o        | Obukhov length                           [m]
      B          | Buoyancy production term                 [m2/s3]
      \alpha_t   | Kinematic turbulent thermal conductivity [m2/s]/[kg/m/s]
      \rho       | Density of fluid                         [-]/[kg/m3]
      \beta      | Thermal expansion coefficient            [1/K]
      T          | Temperature                              [K]
      g          | Gravitational acceleration               [m/s2]
      VSMALL     | A very small number to avoid floating point exceptions
    \endvartable

    References:
    \verbatim
    Definition of the Obukhov length (tag:O):
        Obukhov, A. M. (1946).
        Turbulentnost'v temperaturnoj-neodnorodnoj atmosfere.
        Trudy Inst. Theor. Geofiz. AN SSSR, 1, 95-115.

        Obukhov, A. M. (1971).
        Turbulence in an atmosphere with a non-uniform temperature.
        Boundary-layer meteorology, 2(1), 7-29.
        DOI:10.1007/BF00718085

        Obukhov length. (2016).
        In American Meteorological Society - Glossary of meteorology.
        Retrieved Apr 1, 2020, from glossary.ametsoc.org/wiki/Obukhov_length
    \endverbatim

    Required fields:
    \verbatim
      U           | Velocity                                 [m/s]
      T           | Temperature                              [K]
      nut         | Turbulent viscosity                      [m2/s]
      alphat      | Kinematic turbulent thermal conductivity [m2/s]/[kg/m/s]
      g           | Gravitational acceleration               [m/s2]
    \endverbatim

    \table
      Operand        | Type           | Location
      input          | -              | -
      output file    | -              | -
      output field 1 | volScalarField | $FOAM_CASE/\<time\>/\<ObukhovLength\>
      output field 2 | volScalarField | $FOAM_CASE/\<time\>/\<Ustar\>
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    ObukhovLength1
    {
        // Mandatory entries (unmodifiable)
        type            ObukhovLength;
        libs            (fieldFunctionObjects);

        // Optional entries (runtime modifiable)
        U               U;
        result1         ObukhovLength;
        result2         Ustar;
        rhoRef          1.0;
        kappa           0.4;
        beta            3e-3;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                        | Type | Req'd | Dflt
      type      | Type name: ObukhovLength           | word |  yes  | -
      libs      | Library name: fieldFunctionObjects | word |  yes  | -
      U         | Name of the velocity field         | word |  no   | U
      result1   | Name of the output field for ObukhovLength | word | no <!--
                                                 --> | ObukhovLength
      result2   | Name of the output field for Ustar | word | no    | Ustar
      rhoRef    | Reference density (to convert from kinematic to static <!--
                                       --> pressure) | scalar | no  | 1.0
      kappa     | von Kármán constant                | scalar | no  | 0.40
      beta      | Thermal expansion coefficient      | scalar | no | 3e-3
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink

    Minimal example by using the \c postProcess utility:
    \verbatim
        postProcess -func "ObukhovLength(<UField>)"
    \endverbatim

See also
    - Foam::functionObjects::fvMeshFunctionObject
    - ExtendedCodeGuide::atmosphericModels::functionObjects::ObukhovLength

SourceFiles
    ObukhovLength.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_ObukhovLength_H
#define functionObjects_ObukhovLength_H

#include "fvMeshFunctionObject.H"
#include "turbulentTransportModel.H"
#include "gravityMeshObject.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class ObukhovLength Declaration
\*---------------------------------------------------------------------------*/

class ObukhovLength
:
    public fvMeshFunctionObject
{
protected:

    // Protected Data

        //- Name of velocity field
        word UName_;

        //- Name of the output field for ObukhovLength
        word resultName1_;

        //- Name of the output field for Ustar
        word resultName2_;

        //- Reference density (to convert from kinematic to static pressure)
        scalar rhoRef_;

        //- von Kármán constant [-]
        scalar kappa_;

        //- Thermal expansion coefficient [1/K]
        dimensionedScalar beta_;

        //- Gravitational acceleration vector [m/s2]
        const dimensionedVector g_;


    // Protected Member Functions

        //- Hard-coded Obukhov length field and friction velocity
        //  \return true if fields did not previously exist
        bool calcOL();


public:

    //- Run-time type information
    TypeName("ObukhovLength");


    // Constructors

        //- Construct from Time and dictionary
        ObukhovLength
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        ObukhovLength(const ObukhovLength&) = delete;

        //- No copy assignment
        void operator=(const ObukhovLength&) = delete;


    //- Destructor
    virtual ~ObukhovLength() = default;


    // Member Functions

        //- Remove (checkOut) the output fields from the object registry
        void removeObukhovLength();

        //- Read the data
        virtual bool read(const dictionary& dict);

        //- Calculate the output fields
        virtual bool execute();

        //- Write the output fields
        virtual bool write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh& mpm);

        //- Update for mesh point-motion
        virtual void movePoints(const polyMesh& m);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
