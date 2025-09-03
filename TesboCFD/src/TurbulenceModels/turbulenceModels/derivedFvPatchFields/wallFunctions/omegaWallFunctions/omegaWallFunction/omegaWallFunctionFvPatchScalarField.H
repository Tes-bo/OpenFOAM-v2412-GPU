/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016, 2019 OpenFOAM Foundation
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
    Foam::omegaWallFunctionFvPatchScalarField

Group
    grpWallFunctions

Description
    This boundary condition provides a wall function for the specific
    dissipation rate (i.e. \c omega) and the turbulent kinetic energy
    production contribution (i.e. \c G) for low- and high-Reynolds number
    applications.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        type            omegaWallFunction;

        // Optional entries
        beta1           0.075;

        // Inherited entries
        ...
    }
    \endverbatim

    \table
      Property  | Description                     | Type   | Reqd | Deflt
      type      | Type name: omegaWallFunction    | word   | yes  | -
      beta1     | Model coefficient               | scalar | no   | 0.075
    \endtable

    The inherited entries are elaborated in:
      - \link wallFunctionCoefficients.H \endlink
      - \link wallFunctionBlenders.H \endlink

SourceFiles
    omegaWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef omegaWallFunctionFvPatchScalarField_H
#define omegaWallFunctionFvPatchScalarField_H

#include "fixedValueFvPatchField.H"
#include "wallFunctionCoefficients.H"
#include "wallFunctionBlenders.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class turbulenceModel;

/*---------------------------------------------------------------------------*\
             Class omegaWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class omegaWallFunctionFvPatchScalarField
:
    public fixedValueFvPatchField<scalar>,
    private wallFunctionBlenders
{
protected:

    // Protected Data

        //- Tolerance used in weighted calculations
        static scalar tolerance_;

        //- Initialised flag
        bool initialised_;

        //- Master patch ID
        label master_;

        //- beta1 coefficient
        scalar beta1_;

        //- Wall-function coefficients
        wallFunctionCoefficients wallCoeffs_;

        //- Local copy of turbulence G field
        scalarField G_;

        //- Local copy of turbulence omega field
        scalarField omega_;

        //- List of averaging corner weights
        List<List<scalar>> cornerWeights_;


    // Protected Member Functions

        //- Set the master patch - master is responsible for updating all
        //- wall function patches
        virtual void setMaster();

        //- Create the averaging weights for cells which are bounded by
        //- multiple wall function faces
        virtual void createAveragingWeights();

        //- Helper function to return non-const access to an omega patch
        virtual omegaWallFunctionFvPatchScalarField& omegaPatch
        (
            const label patchi
        );

        //- Main driver to calculate the turbulence fields
        virtual void calculateTurbulenceFields
        (
            const turbulenceModel& turbulence,
            scalarField& G0,
            scalarField& omega0
        );

        //- Calculate the omega and G
        virtual void calculate
        (
            const turbulenceModel& turbulence,
            const List<scalar>& cornerWeights,
            const fvPatch& patch,
            scalarField& G,
            scalarField& omega
        );

        //- Return non-const access to the master patch ID
        virtual label& master()
        {
            return master_;
        }

        //- Write local wall function variables
        void writeLocalEntries(Ostream&) const;


public:

    //- Runtime type information
    TypeName("omegaWallFunction");


    // Constructors

        //- Construct from patch and internal field
        omegaWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        omegaWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- omegaWallFunctionFvPatchScalarField
        //- onto a new patch
        omegaWallFunctionFvPatchScalarField
        (
            const omegaWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        omegaWallFunctionFvPatchScalarField
        (
            const omegaWallFunctionFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        omegaWallFunctionFvPatchScalarField
        (
            const omegaWallFunctionFvPatchScalarField&,
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
    virtual ~omegaWallFunctionFvPatchScalarField() = default;


    // Member Functions

        // Access

            //- Return non-const access to the master's G field
            scalarField& G(bool init = false);

            //- Return non-const access to the master's omega field
            scalarField& omega(bool init = false);


        // Evaluation

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Update the coefficients associated with the patch field
            virtual void updateWeightedCoeffs(const scalarField& weights);

            //- Manipulate matrix
            virtual void manipulateMatrix(fvMatrix<scalar>& matrix);

            //- Manipulate matrix with given weights
            virtual void manipulateMatrix
            (
                fvMatrix<scalar>& matrix,
                const scalarField& weights
            );


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
