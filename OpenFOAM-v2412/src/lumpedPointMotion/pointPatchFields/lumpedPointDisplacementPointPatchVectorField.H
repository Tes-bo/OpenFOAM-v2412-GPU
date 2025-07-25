/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenCFD Ltd.
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
    Foam::lumpedPointDisplacementPointPatchVectorField

Description
    This is the point-patch responsible for managing the force
    integration on a 'lumped-point' basis, waiting for the external
    application, reading back the response from the external program
    and updating the locations of the associated patch points
    accordingly.

    The internal patch type name is 'lumpedPointDisplacement'.

    \heading Patch usage
    Example:
    \verbatim
    walls
    {
        type    lumpedPointDisplacement;
        value   uniform (0 0 0);
        controllers ( controllerName1 controllerName2 );
    }
    \endverbatim

SourceFiles
    lumpedPointDisplacementPointPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef lumpedPointDisplacementPointPatchVectorField_H
#define lumpedPointDisplacementPointPatchVectorField_H

#include "fixedValuePointPatchField.H"
#include "lumpedPointMovement.H"
#include "lumpedPointIOMovement.H"
#include "pointField.H"
#include "pointFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class interpolationWeights;

/*---------------------------------------------------------------------------*\
        Class lumpedPointDisplacementPointPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class lumpedPointDisplacementPointPatchVectorField
:
    public fixedValuePointPatchField<vector>
{
    // Private Data

        //- Names of the movement controller(s) in use
        wordList controllers_;

        //- Times when data (forces) were written
        //  With first=time, second=prevTime
        Tuple2<scalar, scalar> dataWritten_;

        //- Backup method for getting "points0" without a motion solver
        mutable autoPtr<pointIOField> points0Ptr_;


    // Convenience typedefs

        typedef lumpedPointDisplacementPointPatchVectorField patchType;
        typedef DimensionedField<vector, pointMesh> fieldType;


protected:

    // Protected Member Functions

        //- The starting locations (obtained from the motionSolver).
        const pointField& points0() const;

        //- The auto-vivifying singleton for movement.
        const lumpedPointMovement& movement() const;


public:

    //- Runtime type information
    TypeName("lumpedPointDisplacement");


    // Constructors

        //- Construct from patch and internal field
        lumpedPointDisplacementPointPatchVectorField
        (
            const pointPatch& p,
            const DimensionedField<vector, pointMesh>& iF
        );

        //- Construct from patch, internal field and dictionary
        lumpedPointDisplacementPointPatchVectorField
        (
            const pointPatch& p,
            const DimensionedField<vector, pointMesh>& iF,
            const dictionary& dict
        );

        //- Construct by mapping given patch field onto a new patch
        lumpedPointDisplacementPointPatchVectorField
        (
            const lumpedPointDisplacementPointPatchVectorField& rhs,
            const pointPatch& p,
            const DimensionedField<vector, pointMesh>& iF,
            const pointPatchFieldMapper& mapper
        );

        //- Construct as copy setting internal field reference
        lumpedPointDisplacementPointPatchVectorField
        (
            const lumpedPointDisplacementPointPatchVectorField& rhs,
            const DimensionedField<vector, pointMesh>& iF
        );

        //- Return a clone
        virtual autoPtr<pointPatchField<vector>> clone() const
        {
            return pointPatchField<vector>::Clone(*this);
        }

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<vector>> clone
        (
            const DimensionedField<vector, pointMesh>& iF
        ) const
        {
            return pointPatchField<vector>::Clone(*this, iF);
        }


    //- Destructor. De-register movement if in use and managed by this patch
    virtual ~lumpedPointDisplacementPointPatchVectorField();


    // Member Functions

        //- The controller names for this patch
        const wordList& controllers() const
        {
            return controllers_;
        }

        //- Set all patch controls for patches of this type
        static label setPatchControls
        (
            const pointVectorField& pvf,
            const pointField& points0
        );

        //- Set all patch controls for patches of this type
        static label setInterpolators
        (
            const pointVectorField& pvf,
            const pointField& points0
        );

        //- The ids for all patches of this type
        static labelList patchIds(const pointVectorField& pvf);

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
