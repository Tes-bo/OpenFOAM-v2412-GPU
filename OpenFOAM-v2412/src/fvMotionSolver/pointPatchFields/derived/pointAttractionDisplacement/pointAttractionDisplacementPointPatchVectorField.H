/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenCFD Ltd.
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
    Foam::pointAttractionDisplacementPointPatchVectorField

Description
    Displacement by attraction to nearest point. Use in a
    displacementMotionSolver as a bc on the pointDisplacement field.

    Needs:
    - file : name of edgeMesh. Either:
        - .extendedFeatureEdgeMesh (read from  constant/extendedFeatureEdgeMesh)
        - or .eMesh in or a known format (read from constant/triSurface)
    - frozenPointsZone : empty or name of pointZone containing points
                         that do not move
    - velocity : maximum displacement velocity


SourceFiles
    pointAttractionDisplacementPointPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef pointAttractionDisplacementPointPatchVectorField_H
#define pointAttractionDisplacementPointPatchVectorField_H

#include "pointPatchFields.H"
#include "extendedFeatureEdgeMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class pointAttractionDisplacementPointPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class pointAttractionDisplacementPointPatchVectorField
:
    public pointPatchVectorField
{
private:

    // Private data

        //- Maximum velocity
        const vector velocity_;

        //- file
        const fileName featFileName_;

        //- pointZone with frozen points
        const word frozenPointsZone_;

        //- Point searching
        mutable autoPtr<indexedOctree<treeDataPoint>> pointTreePtr_;


    // Private Member Functions

        //- Calculate displacement (w.r.t. points0()) to project onto surface
        void calcProjection(vectorField& displacement) const;

        //- No copy assignment
        void operator=
        (
            const pointAttractionDisplacementPointPatchVectorField&
        ) = delete;


protected:

        const indexedOctree<treeDataPoint>& pointTree() const;


public:

    //- Runtime type information
    TypeName("pointAttraction");


    // Constructors

        //- Construct from patch and internal field
        pointAttractionDisplacementPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        pointAttractionDisplacementPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patch field onto a new patch
        pointAttractionDisplacementPointPatchVectorField
        (
            const pointAttractionDisplacementPointPatchVectorField&,
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct as copy
        pointAttractionDisplacementPointPatchVectorField
        (
            const pointAttractionDisplacementPointPatchVectorField&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchVectorField> clone() const
        {
            return autoPtr<pointPatchVectorField>
            (
                new pointAttractionDisplacementPointPatchVectorField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        pointAttractionDisplacementPointPatchVectorField
        (
            const pointAttractionDisplacementPointPatchVectorField&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchVectorField> clone
        (
            const DimensionedField<vector, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchVectorField>
            (
                new pointAttractionDisplacementPointPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }

    // Member Functions

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
