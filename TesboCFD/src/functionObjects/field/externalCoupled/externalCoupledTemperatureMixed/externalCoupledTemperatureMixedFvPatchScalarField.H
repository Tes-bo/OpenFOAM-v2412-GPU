/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
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
    Foam::externalCoupledTemperatureMixedFvPatchScalarField

Group
    grpCoupledBoundaryConditions

Description
    This boundary condition provides a temperatue interface to an external
    application.

    Values are transferred as plain text files, where OpenFOAM data are
    written as follows:

    \verbatim
        # Patch: <patch name>
        <magSf1> <value1> <qDot1> <htc1>
        <magSf2> <value2> <qDot2> <htc2>
        <magSf3> <value3> <qDot3> <htc2>
        ...
        <magSfN> <valueN> <qDotN> <htcN>
    \endverbatim

    and received as the constituent pieces of the `mixed' condition, i.e.

    \verbatim
        # Patch: <patch name>
        <value1> <gradient1> <valueFracion1>
        <value2> <gradient2> <valueFracion2>
        <value3> <gradient3> <valueFracion3>
        ...
        <valueN> <gradientN> <valueFracionN>
    \endverbatim

    Data is sent/received as a single file for all patches from the directory

    \verbatim
        $FOAM_CASE/<commsDir>
    \endverbatim

    At start-up, the boundary creates a lock file, i.e..

    \verbatim
        OpenFOAM.lock
    \endverbatim

    ... to signal the external source to wait.  During the boundary condition
    update, boundary values are written to file, e.g.

    \verbatim
        <fileName>.out
    \endverbatim

    The lock file is then removed, instructing the external source to take
    control of the program execution.  When ready, the external program
    should create the return values, e.g. to file

    \verbatim
        <fileName>.in
    \endverbatim

    ... and then reinstate the lock file.  The boundary condition will then
    read the return values, and pass program execution back to OpenFOAM.

    To be used in combination with the functionObjects::externalCoupled
    functionObject.

Usage
    \table
        Property | Description                              | Required | Default
        outputTemperature | Output temperature: fluid/wall  | yes |
        htcRefTemperature | Fluid temperature for htc: cell/user  | no | cell
        Tref     | Reference temperature [K] for htc        | conditional |
    \endtable

    The user-specified reference temperature Tref is specified as
    a Foam::Function1 of time but spatially uniform.

SeeAlso
    externalCoupledFunctionObject
    mixedFvPatchField
    externalCoupledMixedFvPatchField

SourceFiles
    externalCoupledTemperatureMixedFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef externalCoupledTemperatureMixedFvPatchScalarField_H
#define externalCoupledTemperatureMixedFvPatchScalarField_H

#include "externalCoupledMixedFvPatchFields.H"
#include "Function1.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class externalCoupledTemperatureMixedFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class externalCoupledTemperatureMixedFvPatchScalarField
:
    public externalCoupledMixedFvPatchField<scalar>
{
    // Data Types

        //- Location for the output temperature
        enum outputTemperatureType
        {
            FLUID,      //!< Use fluid (cell) temperature
            WALL        //!< Use wall (patch) temperature
        };

        //- Names for outputTemperatureType
        static const Enum<outputTemperatureType> outputTemperatureNames;

        //- Reference temperature type for HTC calculation
        enum refTemperatureType
        {
            CELL,       //!< Use cell temperature (default)
            USER        //!< User-specified reference temperature
        };

        //- Names for refTemperatureType
        static const Enum<refTemperatureType> refTemperatureNames;


    // Private Data

        //- Location for the output temperature
        enum outputTemperatureType outTempType_;

        //- Reference temperature type for HTC calculation
        enum refTemperatureType refTempType_;

        //- User-specified reference temperature for HTC calculation
        autoPtr<Function1<scalar>> Tref_;


public:

    //- Runtime type information
    TypeName("externalCoupledTemperature");


    // Constructors

        //- Construct from patch and internal field
        externalCoupledTemperatureMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        externalCoupledTemperatureMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  externalCoupledTemperatureMixedFvPatchScalarField onto a new patch
        externalCoupledTemperatureMixedFvPatchScalarField
        (
            const externalCoupledTemperatureMixedFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        externalCoupledTemperatureMixedFvPatchScalarField
        (
            const externalCoupledTemperatureMixedFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        externalCoupledTemperatureMixedFvPatchScalarField
        (
            const externalCoupledTemperatureMixedFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return fvPatchField<scalar>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return fvPatchField<scalar>::Clone(*this, iF);
        }


    //- Destructor
    virtual ~externalCoupledTemperatureMixedFvPatchScalarField() = default;


    // Member Functions

        //- Write header
        virtual void writeHeader(Ostream& os) const;

        //- Write data
        virtual void writeData(Ostream& os) const;

        //- Read data
        virtual void readData(Istream& is);

        //- Write
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
