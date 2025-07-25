/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2023 OpenCFD Ltd.
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
    Foam::fv::jouleHeatingSource

Group
    grpFvOptionsSources

Description
    Evolves an electrical potential equation

    \f[
        \grad \left( \sigma \grad V \right)
    \f]

    where \f$ V \f$ is electrical potential
    and \f$\sigma\f$ is the electrical current.

    To provide a Joule heating contribution according to:

    Differential form of Joule heating - power per unit volume:

    \f[
        \frac{d(P)}{d(V)} = J \cdot E
    \f]

    where \f$ J \f$ is the current density and \f$ E \f$ the electric field.
    If no magnetic field is present:

    \f[
        J = \sigma E
    \f]

    The electric field given by

    \f[
        E = \grad V
    \f]

    Therefore:

    \f[
        \frac{d(P)}{d(V)} = J \cdot E
                          = (sigma E) \cdot E
                          = (sigma \grad V) \cdot \grad V
    \f]

Usage
    Minimal example by using \c constant/fvOptions:
    \verbatim
    jouleHeatingSource1
    {
        // Mandatory entries (unmodifiable)
        type                 jouleHeatingSource;

        // Mandatory entries (runtime modifiable)
        anisotropicElectricalConductivity true;

        // Optional entries (runtime modifiable)
        T       <Tname>;

        // Conditional mandatory entries (runtime modifiable)

            // when anisotropicElectricalConductivity=true
            coordinateSystem
            {
                origin  (0 0 0);
                e1      (1 0 0);
                e3      (0 0 1);
            }

        // Conditional optional entries (runtime modifiable)

            // when anisotropicElectricalConductivity=false
            // Specify the conductivity as a function of temperature
            // If not supplied, this will be read from the time directory
            sigma           table
            (
                (273        1e5)
                (1000       1e5)
            );

            // when anisotropicElectricalConductivity=true
            sigma           (31900 63800 127600);
            //sigma           table
            //(
            //    (0      (0 0 0))
            //    (1000   (127600 127600 127600))
            //);

        // Mandatory/Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                       | Type  | Reqd | Dflt
      type      | Type name: jouleHeatingSource     | word  | yes  | -
      anisotropicElectricalConductivity | Flag to indicate that <!--
                --> if the electrical conductivity is anisotropic <!--
                -->                                 | bool  | yes  | -
      T         | Name of operand temperature field | word  | no   | T
      sigma     | Electrical conductivity as a function of temperature <!--
                -->                                 | table | no   | -
      coordinateSystem | User-specified coordinate system | coordSystem | no | -
    \endtable

    The inherited entries are elaborated in:
     - \link fvOption.H \endlink

Note
  - \c anisotropicElectricalConductivity=true enables
  anisotropic (vectorial) electrical conductivity.
  - \c anisotropicElectricalConductivity=false enables
  isotropic (scalar) electrical conductivity.
  - The electrical conductivity can be specified using either:
    - If the \c sigma entry is present the electrical conductivity is specified
    as a function of temperature using a \c Function1 type.
    - If not present the \c sigma field will be read from file.
    - If the \c anisotropicElectricalConductivity flag is set to \c true,
    \c sigma should be specified as a vector quantity.

See also
  - Foam::Function1
  - Foam::coordSystem

SourceFiles
    jouleHeatingSource.C
    jouleHeatingSourceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef fv_jouleHeatingSource_H
#define fv_jouleHeatingSource_H

#include "fvOption.H"
#include "Function1.H"
#include "coordinateSystem.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class jouleHeatingSource Declaration
\*---------------------------------------------------------------------------*/

class jouleHeatingSource
:
    public fv::option
{
    // Private Data

        //- Name of temperature field
        word TName_;

        //- Electrical potential field / [V]
        volScalarField V_;

        //- Flag to indicate that the electrical conductivity is anisotropic
        bool anisotropicElectricalConductivity_;

        //- Electrical conductivity as a scalar function of temperature
        autoPtr<Function1<scalar>> scalarSigmaVsTPtr_;

        //- Electrical conductivity as a vector function of temperature
        autoPtr<Function1<vector>> vectorSigmaVsTPtr_;

        //- Coordinate system - used for vectorial electrical conductivity
        autoPtr<coordinateSystem> csysPtr_;

        //- Current time index (used for updating)
        label curTimeIndex_;


    // Private Member Functions

        //- Transform the anisotropic electrical conductivity into global system
        tmp<volSymmTensorField> transformSigma
        (
            const volVectorField& sigmaLocal
        ) const;

        //- Initialise the electrical conductivity field
        template<class Type>
        void initialiseSigma
        (
            const dictionary& dict,
            autoPtr<Function1<Type>>& sigmaVsTPtr
        );

        //- Update the electrical conductivity field
        template<class Type>
        const GeometricField<Type, fvPatchField, volMesh>&
        updateSigma(const autoPtr<Function1<Type>>& sigmaVsTPtr) const;


public:

    //- Runtime type information
    TypeName("jouleHeatingSource");


    // Constructors

        //- Construct from explicit source name and mesh
        jouleHeatingSource
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        jouleHeatingSource(const jouleHeatingSource&) = delete;

        //- No copy assignment
        void operator=(const jouleHeatingSource&) = delete;


    //- Destructor
    virtual ~jouleHeatingSource() = default;


    // Member Functions

        // Evaluation

            //- Add explicit contribution to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const label fieldi
            );


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "jouleHeatingSourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
