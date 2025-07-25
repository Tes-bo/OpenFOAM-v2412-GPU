/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2018 OpenFOAM Foundation
    Copyright (C) 2015-2021 OpenCFD Ltd.
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
    Foam::atmBoundaryLayer

Group
    grpRASBoundaryConditions grpInletBoundaryConditions

Description
    Base class to set log-law type ground-normal inlet boundary conditions for
    wind velocity and turbulence quantities for homogeneous, two-dimensional,
    dry-air, equilibrium and neutral atmospheric boundary layer (ABL) modelling.

    The ground-normal profile expressions are due to \c YGCJ
    (refer to references below) whereat \c RH expressions were generalised:

    \f[
        u = \frac{u^*}{\kappa} \ln \left( \frac{z - d + z_0}{z_0} \right)
    \f]

    \f[
        v = w = 0
    \f]

    \f[
        k = \frac{(u^*)^2}{\sqrt{C_\mu}}
        \sqrt{C_1 \ln \left( \frac{z - d + z_0}{z_0} \right) + C_2}
    \f]

    \f[
        \epsilon = \frac{(u^*)^3}{\kappa (z - d + z_0)}
        \sqrt{C_1 \ln \left( \frac{z - d + z_0}{z_0} \right) + C_2}
    \f]

    \f[
        \omega = \frac{u^*}{\kappa \sqrt{C_\mu}} \frac{1}{z - d + z_0}
    \f]

    \f[
        u^* =
            \frac{u_{ref} \kappa}{\ln\left(\frac{z_{ref} + z_0}{z_0}\right)}
    \f]

    where
    \vartable
      u        | Ground-normal streamwise flow speed profile          [m/s]
      v        | Spanwise flow speed                                  [m/s]
      w        | Ground-normal flow speed                             [m/s]
      k        | Ground-normal turbulent kinetic energy (TKE) profile [m^2/s^2]
      \epsilon | Ground-normal TKE dissipation rate profile           [m^2/s^3]
      \omega   | Ground-normal specific dissipation rate profile      [m^2/s^3]
      u^*      | Friction velocity                                    [m/s]
      \kappa   | von Kármán constant                                  [-]
      C_\mu    | Empirical model constant                             [-]
      z        | Ground-normal coordinate component                   [m]
      d        | Ground-normal displacement height                    [m]
      z_0      | Aerodynamic roughness length                         [m]
      u_{ref}  | Reference mean streamwise wind speed at \f$z_{ref}\f$ [m/s]
      z_{ref}  | Reference height being used in \f$u^*\f$ estimations [m]
      C_1      | Curve-fitting coefficient for \c YGCJ profiles       [-]
      C_2      | Curve-fitting coefficient for \c YGCJ profiles       [-]
    \endvartable

    Reference:
    \verbatim
        The ground-normal profile expressions (tag:RH):
            Richards, P. J., & Hoxey, R. P. (1993).
            Appropriate boundary conditions for computational wind
            engineering models using the k-ε turbulence model.
            In Computational Wind Engineering 1 (pp. 145-153).
            DOI:10.1016/B978-0-444-81688-7.50018-8

        Modifications to preserve the profiles downstream (tag:HW):
            Hargreaves, D. M., & Wright, N. G. (2007).
            On the use of the k–ε model in commercial CFD software
            to model the neutral atmospheric boundary layer.
            Journal of wind engineering and
            industrial aerodynamics, 95(5), 355-369.
            DOI:10.1016/j.jweia.2006.08.002

        Expression generalisations to allow height
        variation for turbulence quantities (tag:YGCJ):
            Yang, Y., Gu, M., Chen, S., & Jin, X. (2009).
            New inflow boundary conditions for modelling the neutral equilibrium
            atmospheric boundary layer in computational wind engineering.
            J. of Wind Engineering and Industrial Aerodynamics, 97(2), 88-95.
            DOI:10.1016/j.jweia.2008.12.001

        The generalised ground-normal profile expression for omega (tag:YGJ):
            Yang, Y., Gu, M., & Jin, X., (2009).
            New inflow boundary conditions for modelling the
            neutral equilibrium atmospheric boundary layer in SST k-ω model.
            In: The Seventh Asia-Pacific Conference on Wind Engineering,
            November 8-12, Taipei, Taiwan.

        Theoretical remarks (tag:E):
            Emeis, S. (2013).
            Wind Energy Meteorology: Atmospheric
            Physics for Wind Power Generation.
            Springer-Verlag Berlin Heidelberg.
            DOI:10.1007/978-3-642-30523-8
    \endverbatim

Usage
    Example of the entries provided for the inherited boundary conditions:
    \verbatim
    inlet
    {
        // Mandatory and other optional entries
        ...

        // Mandatory (inherited) entries (runtime modifiable)
        flowDir         (1 0 0);
        zDir            (0 0 1);
        Uref            10.0;
        Zref            0.0;
        z0              uniform 0.1;
        d               uniform 0.0;

        // Optional (inherited) entries (unmodifiable)
        kappa           0.41;
        Cmu             0.09;
        initABL         true;
        phi             phi;
        C1              0.0;
        C2              1.0;

        // Conditional mandatory (inherited) entries (runtime modifiable)
        value           uniform 0;    // when initABL=false
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description              | Type          | Reqd  | Deflt
      flowDir   | Flow direction           | Function1<vector> | yes | -
      zDir      | Ground-normal direction  | Function1<vector> | yes | -
      Uref      | Reference mean streamwise flow speed being used in <!--
             --> \f$u^*\f$ estimations [m/s] | Function1<scalar> | yes | -
      Zref      | Reference height being used in \f$u^*\f$ estimations [m]  <!--
                                       --> | Function1<scalar> | yes | -
      z0        | Surface roughness length [m] <!--
                                       --> | PatchFunction1<scalar> | yes | -
      d         | Displacement height [m] - see Notes <!--
                                       --> | PatchFunction1<scalar> | yes | -
      kappa     | von Kármán constant      | scalar | no  | 0.41
      Cmu       | Empirical model constant | scalar | no  | 0.09
      initABL   | Flag to initialise profiles with the theoretical <!--
            --> ABL expressions, otherwise use "value" list        <!--
                                       --> | bool   | no  | true
      value     | ABL profile content when initABL=false <!--
                                       --> | scalarList | conditional | -
      phi       | Name of the flux field   | word         | no        | phi
      C1        | Curve-fitting coefficient \c YGCJ profiles | scalar | no | 0.0
      C2        | Curve-fitting coefficient \c YGCJ profiles | scalar | no | 1.0
    \endtable

Note
    - The \c RH expressions are special cases of those in \c YGCJ when \c C1=0
    and \c C2=1. Both \c C1 and \c C2 can be determined by nonlinear fitting
    of (\c YGCJ:Eqs. 19-20) with an experimental dataset for \c k. By default,
    \c atmBoundaryLayerInlet boundary conditions compute \c RH expressions.
    - \c z is the ground-normal height relative to the global minimum height
    of the inlet patch; therefore, the minimum of \c z is always zero
    irrespective of the absolute z-coordinate of the computational patch.
    - The derived ABL expressions automatically satisfy the simplified transport
    equation for \c k. Yet the same expressions only satisfy the simplified
    transport equation for \c epsilon when the model constants \c sigmaEpsilon
    is 1.11 with \c kappa=0.4 (\c HW:p. 358).
    - \c atmBoundaryLayerInlet boundary conditions inherit \c inletOutlet
    traits, so that a given inlet condition can be supplied from all sides of
    the domain, e.g. a ground-normal cylinder domain having a single
    inlet/outlet boundary where the changes between inlet and outlet depend
    on the wind direction and patch normals, so that any change in inflow
    orientation can be handled with the same mesh.
    - \c d is the displacement height, and "is relevant for flows over forests
    and cities" (E:p. 28). "The displacement height gives the vertical
    displacement of the entire flow regime over areas which are densely covered
    with obstacles such as trees or buildings" (E:p. 28).

See also
    - Foam::atmBoundaryLayerInletVelocityFvPatchVectorField
    - Foam::atmBoundaryLayerInletKFvPatchScalarField
    - Foam::atmBoundaryLayerInletEpsilonFvPatchScalarField
    - Foam::atmBoundaryLayerInletOmegaFvPatchScalarField

SourceFiles
    atmBoundaryLayer.C

\*---------------------------------------------------------------------------*/

#ifndef atmBoundaryLayer_H
#define atmBoundaryLayer_H

#include "fvPatchFields.H"
#include "Function1.H"
#include "PatchFunction1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class atmBoundaryLayer Declaration
\*---------------------------------------------------------------------------*/

class atmBoundaryLayer
{
protected:

    // Protected Data

        //- Flag to initialise profiles with the theoretical ABL expressions,
        //- otherwise initialises by using "value" entry content
        bool initABL_;


private:

    // Private Data

        //- von Kármán constant
        const scalar kappa_;

        //- Empirical model constant
        const scalar Cmu_;

        //- Curve-fitting coefficient
        const scalar C1_;

        //- Curve-fitting coefficient
        const scalar C2_;

        //- Minimum coordinate vector of this patch
        const vector ppMin_;

        //- Reference to the time database
        const Time& time_;

        //- Reference to the patch
        const polyPatch& patch_;

        //- Streamwise flow direction
        autoPtr<Function1<vector>> flowDir_;

        //- Direction of the ground-normal coordinate
        autoPtr<Function1<vector>> zDir_;

        //- Reference mean streamwise flow speed being used in Ustar estimations
        autoPtr<Function1<scalar>> Uref_;

        //- Reference height being used in Ustar estimations
        autoPtr<Function1<scalar>> Zref_;

        //- Surface roughness length
        autoPtr<PatchFunction1<scalar>> z0_;

        //- Displacement height
        autoPtr<PatchFunction1<scalar>> d_;


public:

    // Constructors

        //- Construct null from time
        atmBoundaryLayer(const Time& time, const polyPatch& pp);

        //- Construct from the time database and dictionary
        atmBoundaryLayer
        (
            const Time& time,
            const polyPatch& pp,
            const dictionary& dict
        );

        //- Construct by mapping given atmBoundaryLayer onto a new patch
        atmBoundaryLayer
        (
            const atmBoundaryLayer& abl,
            const fvPatch& patch,
            const fvPatchFieldMapper& mapper
        );

        //- Construct as copy
        atmBoundaryLayer(const atmBoundaryLayer&);


    // Member Functions

        // Access

            //- Return flow direction
            vector flowDir() const;

            //- Return the ground-normal direction
            vector zDir() const;

            //- Return friction velocity
            tmp<scalarField> Ustar(const scalarField& z0) const;


        // Mapping

            //- Map (and resize as needed) from self given a mapping object
            void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            void rmap(const atmBoundaryLayer&, const labelList&);


        // Evaluation

            //- Return the velocity distribution for the ATM
            tmp<vectorField> U(const vectorField& pCf) const;

            //- Return the turbulent kinetic energy distribution for the ATM
            tmp<scalarField> k(const vectorField& pCf) const;

            //- Return the turbulent dissipation rate distribution for the ATM
            tmp<scalarField> epsilon(const vectorField& pCf) const;

            //- Return the specific dissipation rate distribution for the ATM
            tmp<scalarField> omega(const vectorField& pCf) const;


        //- Write
        void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
