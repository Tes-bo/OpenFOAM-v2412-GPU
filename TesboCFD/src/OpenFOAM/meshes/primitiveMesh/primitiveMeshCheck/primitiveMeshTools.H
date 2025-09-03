/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd.
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
    Foam::primitiveMeshTools

Description
    Collection of static functions operating on primitiveMesh (mainly checks).

SourceFiles
    primitiveMeshTools.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_primitiveMeshTools_H
#define Foam_primitiveMeshTools_H

#include "bitSet.H"
#include "primitiveFields.H"
#include "pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class cell;
class face;
class primitiveMesh;

/*---------------------------------------------------------------------------*\
                     Class primitiveMeshTools Declaration
\*---------------------------------------------------------------------------*/

class primitiveMeshTools
{
public:

    //- Update face centres and areas for the faces in the set faceIDs
    static void updateFaceCentresAndAreas
    (
        const primitiveMesh& mesh,
        const UList<label>& faceIDs,
        const pointField& p,
        vectorField& fCtrs,
        vectorField& fAreas
    );

    //- Update cell centres and volumes for the cells in the set cellIDs
    static void updateCellCentresAndVols
    (
        const primitiveMesh& mesh,
        const vectorField& fCtrs,
        const vectorField& fAreas,
        const List<label>& cellIDs,
        const List<cell>& cells,
        vectorField& cellCtrs_s,
        scalarField& cellVols_s
    );

    //- Calculate face centres and areas for specified faces.
    //  Adjusts the lengths of centres and area normals if required.
    static void makeFaceCentresAndAreas
    (
        const UList<face>& faces,   //!< The faces to consider
        const pointField& p,        //!< Support points for faces
        vectorField& fCtrs,         //!< [out] The face centres
        vectorField& fAreas         //!< [out] The face area normals
    );

    //- Calculate face centres and areas for all mesh faces.
    //  Adjusts the lengths of centres and area normals if required.
    static void makeFaceCentresAndAreas
    (
        const primitiveMesh& mesh,  //!< The mesh faces
        const pointField& p,        //!< Support points for faces
        vectorField& fCtrs,         //!< [out] The face centres
        vectorField& fAreas         //!< [out] The face area normals
    );

    //- Calculate cell centres and volumes from face properties
    static void makeCellCentresAndVols
    (
        const primitiveMesh& mesh,
        const vectorField& fCtrs,
        const vectorField& fAreas,
        vectorField& cellCtrs,
        scalarField& cellVols
    );

    //- Generate non-orthogonality field (internal faces only)
    static tmp<scalarField> faceOrthogonality
    (
        const primitiveMesh& mesh,
        const vectorField& fAreas,
        const vectorField& cellCtrs
    );

    //- Generate face pyramid volume fields
    static void facePyramidVolume
    (
        const primitiveMesh& mesh,
        const pointField& points,
        const vectorField& cellCtrs,
        scalarField& ownPyrVol,
        scalarField& neiPyrVol
    );

    //- Generate skewness field
    static tmp<scalarField> faceSkewness
    (
        const primitiveMesh& mesh,
        const pointField& points,
        const vectorField& fCtrs,
        const vectorField& fAreas,
        const vectorField& cellCtrs
    );

    //- Generate cell openness and cell aspect ratio field
    static void cellClosedness
    (
        const primitiveMesh& mesh,
        const Vector<label>& meshD,
        const vectorField& areas,
        const scalarField& vols,
        scalarField& openness,
        scalarField& aratio
    );

    //- Generate face concavity field. Returns per face the (sin of the)
    //- most concave angle between two consecutive edges
    static tmp<scalarField> faceConcavity
    (
        const scalar maxSin,
        const primitiveMesh& mesh,
        const pointField& p,
        const vectorField& faceAreas
    );

    //- Generate face flatness field. Compares the individual triangles'
    //- normals against the face average normal. Between 0 (fully warped)
    //- and 1 (fully flat)
    static tmp<scalarField> faceFlatness
    (
        const primitiveMesh& mesh,
        const pointField& p,
        const vectorField& fCtrs,
        const vectorField& faceAreas
    );

    //- Generate edge alignment field. Is per face the minimum aligned edge
    //- (does not use edge addressing)
    static tmp<scalarField> edgeAlignment
    (
        const primitiveMesh& mesh,
        const Vector<label>& directions,
        const pointField& p
    );

    //- Generate cell determinant field. Normalised to 1 for an internal cube.
    static tmp<scalarField> cellDeterminant
    (
        const primitiveMesh& mesh,
        const Vector<label>& directions,
        const vectorField& faceAreas,
        const bitSet& internalOrCoupledFace
    );


    // Helpers: single face check

        //- Skewness of single face
        static scalar faceSkewness
        (
            const UList<face>& faces,   //!< The faces to consider
            const pointField& p,        //!< Support points for faces
            const vectorField& fCtrs,
            const vectorField& fAreas,

            const label facei,
            const point& ownCc,
            const point& neiCc
        );

        //- Skewness of single boundary face
        static scalar boundaryFaceSkewness
        (
            const UList<face>& faces,   //!< The faces to consider
            const pointField& p,        //!< Support points for faces
            const vectorField& fCtrs,
            const vectorField& fAreas,

            const label facei,
            const point& ownCc
        );

        //- Skewness of single face
        static scalar faceSkewness
        (
            const primitiveMesh& mesh,  //!< The mesh faces
            const pointField& p,        //!< Support points for faces
            const vectorField& fCtrs,
            const vectorField& fAreas,

            const label facei,
            const point& ownCc,
            const point& neiCc
        );

        //- Skewness of single boundary face
        static scalar boundaryFaceSkewness
        (
            const primitiveMesh& mesh,  //!< The mesh faces
            const pointField& p,        //!< Support points for faces
            const vectorField& fCtrs,
            const vectorField& fAreas,

            const label facei,
            const point& ownCc
        );

        //- Orthogonality of single face
        static scalar faceOrthogonality
        (
            const point& ownCc,
            const point& neiCc,
            const vector& s
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
