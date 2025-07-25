/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
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
    Foam::functionObjects::turbulenceFields

Group
    grpFieldFunctionObjects

Description
    Computes various turbulence-related quantities that are not typically
    output during calculations, and stores/writes them on the mesh database
    for further manipulation.

    Fields are stored as copies of the original with a user-defined prefix
    e.g. a prefix of \c turbulenceModel yields the following for field \c R:

    \verbatim
      turbulenceModel:R
    \endverbatim

    Operands:
    \table
      Operand           | Type             | Location
      input             | -                | -
      output file       | -                | -
      output field      | vol\<Type\>Field | $FOAM_CASE/\<time\>/\<outField\>
    \endtable

    where \c \<Type\>=Scalar/SymmTensor.

    References:
    \verbatim
        Estimation expressions for L (tag:P), Eq. 10.37:
            Pope, S. B. (2000).
            Turbulent flows.
            Cambridge, UK: Cambridge Univ. Press
            DOI:10.1017/CBO9780511840531
    \endverbatim

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    turbulenceFields1
    {
        // Mandatory entries (unmodifiable)
        type        turbulenceFields;
        libs        (fieldFunctionObjects);

        // Mandatory entries (runtime modifiable)

            // Either of the below
            // Option-1
            fields      (R devRhoReff);

            // Option-2
            field       R;

        // Optional entries (runtime modifiable)
        prefix          <word>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                           | Type | Reqd | Deflt
      type         | Type name: turbulenceFields           | word | yes  | -
      libs         | Library name: fieldFunctionObjects    | word | yes  | -
      fields       | Names of fields to store (see below)  | wordList | yes | -
      field        | Name of a field to store (see below)  | word | yes  | -
      prefix  | Name of output-field prefix | word | no | turbulenceProperties
    \endtable

    where \c fields can include:
    \verbatim
      k           | turbulent kinetic energy
      epsilon     | turbulent kinetic energy dissipation rate
      omega       | specific dissipation rate
      nuTilda     | modified turbulent viscosity
      nut         | turbulent viscosity (incompressible)
      nuEff       | effective turbulent viscosity (incompressible)
      mut         | turbulent viscosity (compressible)
      muEff       | effective turbulent viscosity (compressible)
      alphat      | turbulence thermal diffusivity (compressible)
      alphaEff    | effective turbulence thermal diffusivity (compressible)
      R           | Reynolds stress tensor
      devReff     | deviatoric part of the effective Reynolds stress
      devRhoReff  | divergence of the Reynolds stress
      L           | integral-length scale / mixing-length scale
      I           | turbulence intensity
    \endverbatim

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink

    Minimal example by using the \c postProcess utility:
    \verbatim
        <solver> -postProcess -func turbulenceFields
    \endverbatim

Note
  - Multiphase applications are not supported.
  - The governing expression of \c nuTilda is
    an approximation based on a dimensional analysis.

See also
  - Foam::functionObject
  - Foam::functionObjects::fvMeshFunctionObject
  - ExtendedCodeGuide::functionObjects::field::turbulenceFields

SourceFiles
    turbulenceFields.C
    turbulenceFieldsTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_turbulenceFields_H
#define functionObjects_turbulenceFields_H

#include "fvMeshFunctionObject.H"
#include "HashSet.H"
#include "Enum.H"
#include "volFieldsFwd.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                      Class turbulenceFields Declaration
\*---------------------------------------------------------------------------*/

class turbulenceFields
:
    public fvMeshFunctionObject
{
public:

    // Public Enumerations

        //- Options for the turbulence fields (compressible)
        enum compressibleField
        {
            cfK,            //!< "Turbulent kinetic energy"
            cfEpsilon,      //!< "Turbulent kinetic energy dissipation rate"
            cfOmega,        //!< "Specific dissipation rate"
            cfNuTilda,      //!< "Modified turbulent viscosity"
            cfMut,          //!< "Turbulent dynamic viscosity"
            cfMuEff,        //!< "Effective turbulent dynamic viscosity"
            cfAlphat,       //!< "Turbulence thermal diffusivity"
            cfAlphaEff,     //!< "Effective turbulence thermal diffusivity"
            cfR,            //!< "Reynolds stress tensor"
            cfDevRhoReff,   //!< "Divergence of the Reynolds stress"
            cfL,            //!< "Integral-length/Mixing-length scale"
            cfI,            //!< "Turbulence intensity"
            cfLESRegion,    //!< "DES model LES region indicator field"
            cffd            //!< "DES model shielding function"
        };

        //- Names for compressibleField turbulence fields
        static const Enum<compressibleField> compressibleFieldNames_;

        //- Options for the turbulence fields (incompressible)
        enum incompressibleField
        {
            ifK,            //!< "Turbulent kinetic energy"
            ifEpsilon,      //!< "Turbulent kinetic energy dissipation rate"
            ifOmega,        //!< "Specific dissipation rate"
            ifNuTilda,      //!< "Modified turbulent viscosity"
            ifNut,          //!< "Turbulent viscosity"
            ifNuEff,        //!< "Effective turbulent viscosity"
            ifR,            //!< "Reynolds stress tensor"
            ifDevReff,   //!< "Deviatoric part of the effective Reynolds stress"
            ifL,            //!< "Integral-length/Mixing-length scale"
            ifI,            //!< "Turbulence intensity"
            ifLESRegion,    //!< "DES model LES region indicator field"
            iffd            //!< "DES model shielding function"
        };

        //- Names for incompressibleField turbulence fields
        static const Enum<incompressibleField> incompressibleFieldNames_;

        //- Name of the turbulence properties dictionary
        static const word modelName_;


protected:

    // Protected Data

        //- Flag to track initialisation
        bool initialised_;

        //- Name of output-field prefix
        word prefix_;

        //- Fields to load
        wordHashSet fieldSet_;


    // Protected Member Functions

        //- Unset duplicate fields already registered by other function objects
        void initialise();

        //- Return true if compressible turbulence model is identified
        bool compressible();

        //- Process the turbulence field
        template<class Type>
        void processField
        (
            const word& fieldName,
            const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvalue
        );

        //- Return nuTilda calculated from k and omega
        template<class Model>
        tmp<volScalarField> nuTilda(const Model& model) const;

        //- Return integral length scale, L, calculated from k and epsilon
        template<class Model>
        tmp<volScalarField> L(const Model& model) const;

        //- Return turbulence intensity, I, calculated from k and U
        template<class Model>
        tmp<volScalarField> I(const Model& model) const;


public:

    //- Runtime type information
    TypeName("turbulenceFields");


    // Constructors

        //- Construct from Time and dictionary
        turbulenceFields
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        turbulenceFields(const turbulenceFields&) = delete;

        //- No copy assignment
        void operator=(const turbulenceFields&) = delete;


    //- Destructor
    virtual ~turbulenceFields() = default;


    // Member Functions

        //- Read the controls
        virtual bool read(const dictionary&);

        //- Calculate turbulence fields
        virtual bool execute();

        //- Do nothing.
        //  The turbulence fields are registered and written automatically
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "turbulenceFieldsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
