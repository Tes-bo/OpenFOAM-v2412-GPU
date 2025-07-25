/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2020 PCOpt/NTUA
    Copyright (C) 2013-2020 FOSS GP
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::adjointWallVelocityFvPatchVectorField

Description
    Adjoint wall velocity boundary condition.  If nutUSpaldingWallFunction is
    employed in the flow solution, the corresponding adjoint wall function is
    used. Otherwise, the typical low-Re boundary condition is applied

    Reference:
    \verbatim
        For both the low- and high-Re variants

            Papoutsis-Kiachagias, E. M., & Giannakoglou, K. C. (2014).
            Continuous Adjoint Methods for Turbulent Flows, Applied to Shape
            and Topology Optimization: Industrial Applications.
            Archives of Computational Methods in Engineering, 23(2), 255-299.
            http://doi.org/10.1007/s11831-014-9141-9
    \endverbatim

SourceFiles
    adjointWallVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef adjointWallVelocityFvPatchVectorField_H
#define adjointWallVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "adjointBoundaryConditions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class adjointWallVelocity Declaration
\*---------------------------------------------------------------------------*/

class adjointWallVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField,
    public adjointVectorBoundaryCondition
{
private:

    // Private Data

        scalar kappa_;
        scalar E_;


public:

    //- Runtime type information
    TypeName("adjointWallVelocity");


    // Constructors

        //- Construct from patch and internal field
        adjointWallVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        adjointWallVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given adjointWallVelocityFvPatchVectorField
        //- onto a new patch
        adjointWallVelocityFvPatchVectorField
        (
            const adjointWallVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        adjointWallVelocityFvPatchVectorField
        (
            const adjointWallVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<vector>> clone() const
        {
            return fvPatchField<vector>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<vector>> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return fvPatchField<vector>::Clone(*this, iF);
        }


    // Member functions

        //- In case of High-Re runs based on the nutUSpaldingWallFunction
        //- add source terms in the first cell centre off the wall
        virtual void manipulateMatrix(fvMatrix<vector>& matrix);

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
