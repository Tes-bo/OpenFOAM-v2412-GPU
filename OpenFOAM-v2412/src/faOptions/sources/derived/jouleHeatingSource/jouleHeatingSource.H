/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::fa::jouleHeatingSource

Group
    grpFaOptionsSources

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
    Minimal example by using \c constant/faOptions:
    \verbatim
    jouleHeatingSource1
    {
        // Mandatory entries (unmodifiable)
        type                 jouleHeatingSource;

        // Mandatory entries (runtime modifiable)
        anisotropicElectricalConductivity true;

        // Optional entries (runtime modifiable)
        T                   <Tname>;
        nIter               -1;

        // Conditional mandatory entries (runtime modifiable)

            // when the entry "sigma" is present
            sigma           <Function1>;

            // when when the entry "sigma" is not present
            // read "sigma" from file

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
      nIter     | Number of iterations for electrical potential equation <!--
                --> solution                        | label | no   | -1
    \endtable

    The inherited entries are elaborated in:
     - \link faOption.H \endlink
     - \link faceSetOption.H \endlink

Note
  - \c anisotropicElectricalConductivity=true enables
  anisotropic (tensorial) electrical conductivity.
  - \c anisotropicElectricalConductivity=false enables
  isotropic (scalar) electrical conductivity.
  - The electrical conductivity can be specified using either:
    - If the \c sigma entry is present the electrical conductivity is specified
    as a function of temperature using a \c Function1 type.
    - If not present the \c sigma field will be read from file.
    - If the \c anisotropicElectricalConductivity flag is set to \c true,
    \c sigma should be specified as a tensor quantity.

See also
  - Foam::Function1

SourceFiles
    jouleHeatingSource.C
    jouleHeatingSourceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_fa_jouleHeatingSource_H
#define Foam_fa_jouleHeatingSource_H

#include "faOption.H"
#include "Function1.H"
#include "areaFields.H"
#include "faceSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fa
{

/*---------------------------------------------------------------------------*\
                        Class jouleHeatingSource Declaration
\*---------------------------------------------------------------------------*/

class jouleHeatingSource
:
    public fa::faceSetOption
{
    // Private Data

        //- Name of temperature field
        word TName_;

        //- Electrical potential field / [V]
        areaScalarField V_;

        //- Electrical conductivity as a scalar function of temperature
        autoPtr<Function1<scalar>> scalarSigmaVsTPtr_;

        //- Electrical conductivity as a tensor function of temperature
        autoPtr<Function1<tensor>> tensorSigmaVsTPtr_;

        //- Current time index (used for updating)
        label curTimeIndex_;

        //- Number of iterations for electrical potential equation solution
        label nIter_;

        //- Flag to indicate that the electrical conductivity is anisotropic
        bool anisotropicElectricalConductivity_;


    // Private Member Functions

        //- Initialise the electrical conductivity field
        template<class Type>
        void initialiseSigma
        (
            const dictionary& dict,
            autoPtr<Function1<Type>>& sigmaVsTPtr
        );

        //- Update the electrical conductivity field
        template<class Type>
        const GeometricField<Type, faPatchField, areaMesh>&
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
                const areaScalarField& h,
                const areaScalarField& rho,
                faMatrix<scalar>& eqn,
                const label fieldi
            );


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "jouleHeatingSourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
