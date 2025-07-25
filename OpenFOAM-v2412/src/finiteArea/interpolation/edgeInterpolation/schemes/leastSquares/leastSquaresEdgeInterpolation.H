/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::leastSquaresEdgeInterpolation

Description
    Least-squares edge interpolation scheme class, from
    face centers to points then from points to edges.

    References:
    \verbatim
        Governing equations (tag:P):
            Pesci, C. (2019).
            Computational analysis of fluid interfaces
            influenced by soluble surfactant.
            Darmstadt, Technische Universität. PhD thesis.
            URI:https://tuprints.ulb.tu-darmstadt.de/id/eprint/9303
    \endverbatim

SourceFiles
    leastSquaresEdgeInterpolationMake.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_leastSquaresEdgeInterpolation_H
#define Foam_leastSquaresEdgeInterpolation_H

#include "edgeInterpolationScheme.H"
#include "areaFields.H"
#include "edgeFields.H"
#include "leastSquaresFaGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class leastSquaresEdgeInterpolation Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class leastSquaresEdgeInterpolation
:
    virtual public edgeInterpolationScheme<Type>
{
public:

    //- Runtime type information
    TypeName("leastSquares");


    // Generated Methods

        //- No copy construct
        leastSquaresEdgeInterpolation(const leastSquaresEdgeInterpolation&)
            = delete;

        //- No copy assignment
        void operator=(const leastSquaresEdgeInterpolation&) = delete;


    // Constructors

        //- Construct from mesh
        leastSquaresEdgeInterpolation(const faMesh& mesh)
        :
            edgeInterpolationScheme<Type>(mesh)
        {}

        //- Construct from Istream
        leastSquaresEdgeInterpolation(const faMesh& mesh, Istream&)
        :
            edgeInterpolationScheme<Type>(mesh)
        {}

        //- Construct from faceFlux and Istream
        leastSquaresEdgeInterpolation
        (
            const faMesh& mesh,
            const edgeScalarField&,
            Istream&
        )
        :
            edgeInterpolationScheme<Type>(mesh)
        {}


    // Member Functions

        //- Return the edge-interpolate of the given face field
        tmp<GeometricField<Type, faePatchField, edgeMesh>>
        interpolate
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        ) const;

        //- Return the interpolation weighting factors
        tmp<edgeScalarField> weights
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        ) const
        {
            return this->mesh().edgeInterpolation::weights();
        }
};


template<class Type>
tmp<GeometricField<Type, faePatchField, edgeMesh>>
leastSquaresEdgeInterpolation<Type>::interpolate
(
    const GeometricField<Type, faPatchField, areaMesh>& af
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;
    typedef GeometricField<GradType, faPatchField, areaMesh> GradFieldType;
    typedef GeometricField<Type, faePatchField, edgeMesh> EdgeFieldType;
    typedef Field<Type> FieldType;

    const faMesh& mesh = af.mesh();

    const labelList& meshPts = mesh.patch().meshPoints();
    const labelListList& pointFaceAddr = mesh.patch().pointFaces();
    const areaVectorField& C = mesh.areaCentres();

    tmp<GradFieldType> tgradAf= fa::leastSquaresFaGrad<Type>(mesh).grad(af);
    const GradFieldType& gradAf = tgradAf.cref();

    // Field at surface points, Pi (P:p. 29-30)
    auto tPi = tmp<FieldType>::New(meshPts.size(), Zero);
    FieldType& Pi = tPi.ref();

    forAll(meshPts, i)
    {
        const label nFaces = pointFaceAddr[i].size();

        FieldType Pij(nFaces);
        vectorField dPC(nFaces);

        for (label facei = 0; facei < nFaces; ++facei)
        {
            const label j = pointFaceAddr[i][facei];

            // Euclidean distance between a point and face centres (P:Fig. 3.3)
            dPC[facei] = mesh.points()[i] - C[j];

            // (P:Eq. 3.14)
            Pij[facei] = af[j] + (dPC[facei] & gradAf[j]);
        }

        const scalarField magdPC(max(mag(dPC), SMALL));

        // (P:Eq. 3.15)
        Pi[i] = sum(Pij/magdPC)/sum(scalar(1)/magdPC);
    }

    tgradAf.clear();


    auto tinterp = tmp<EdgeFieldType>::New
    (
        IOobject
        (
            "interpolate(" + af.name() + ')',
            af.instance(),
            af.db()
        ),
        af.mesh(),
        af.dimensions()
    );
    EdgeFieldType& interp = tinterp.ref();
    FieldType& interpEdges = interp.primitiveFieldRef();


    const edgeList& faEdges = mesh.edges();
    const label nInternalEdges = mesh.nInternalEdges();

    // Internal field
    for (label edgei = 0; edgei < nInternalEdges; ++edgei)
    {
        // (P:Eq. 3.16)
        interpEdges[edgei] =
            0.5*
            (
                Pi[faEdges[edgei].start()]
              + Pi[faEdges[edgei].end()]
            );
    }

    // Boundary field from boundary condition [fixedValue]
    // Coupled boundary fields are not handled
    forAll(interp.boundaryField(), patchi)
    {
        interp.boundaryFieldRef()[patchi] = af.boundaryField()[patchi];
    }

    return tinterp;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
