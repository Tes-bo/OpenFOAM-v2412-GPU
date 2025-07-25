/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "edgeInterpolationScheme.H"
#include "areaFields.H"
#include "edgeFields.H"
#include "coupledFaPatchField.H"

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::edgeInterpolationScheme<Type>>
Foam::edgeInterpolationScheme<Type>::New
(
    const faMesh& mesh,
    Istream& schemeData
)
{
    if (edgeInterpolation::debug)
    {
        InfoInFunction
            << "constructing edgeInterpolationScheme<Type>"
            << endl;
    }

    if (schemeData.eof())
    {
        FatalIOErrorInFunction(schemeData)
            << "Discretisation scheme not specified" << nl << nl
            << "Valid schemes are :" << nl
            << MeshConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    const word schemeName(schemeData);

    auto* ctorPtr = MeshConstructorTable(schemeName);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            schemeData,
            "discretisation",
            schemeName,
            *MeshConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return ctorPtr(mesh, schemeData);
}


template<class Type>
Foam::tmp<Foam::edgeInterpolationScheme<Type>>
Foam::edgeInterpolationScheme<Type>::New
(
    const faMesh& mesh,
    const edgeScalarField& faceFlux,
    Istream& schemeData
)
{
    if (edgeInterpolation::debug)
    {
        InfoInFunction
            << "constructing edgeInterpolationScheme<Type>"
            << endl;
    }

    if (schemeData.eof())
    {
        FatalIOErrorInFunction(schemeData)
            << "Discretisation scheme not specified"
            << endl << endl
            << "Valid schemes are :" << endl
            << MeshConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    const word schemeName(schemeData);

    auto* ctorPtr = MeshFluxConstructorTable(schemeName);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            schemeData,
            "discretisation",
            schemeName,
            *MeshFluxConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return ctorPtr(mesh, faceFlux, schemeData);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::edgeInterpolationScheme<Type>::~edgeInterpolationScheme()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::faePatchField, Foam::edgeMesh>>
Foam::edgeInterpolationScheme<Type>::interpolate
(
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const tmp<edgeScalarField>& tlambdas,
    const tmp<edgeScalarField>& tys
)
{
    if (edgeInterpolation::debug)
    {
        InfoInFunction
            << "interpolating "
            << vf.type() << " "
            << vf.name()
            << " from areas to edges "
               "without explicit correction"
            << endl;
    }

    const edgeScalarField& lambdas = tlambdas();
    const edgeScalarField& ys = tys();

    const Field<Type>& vfi = vf.internalField();
    const scalarField& lambda = lambdas.internalField();
    const scalarField& y = ys.internalField();

    const faMesh& mesh = vf.mesh();
    const labelUList& P = mesh.owner();
    const labelUList& N = mesh.neighbour();

    tmp<GeometricField<Type, faePatchField, edgeMesh>> tsf
    (
        new GeometricField<Type, faePatchField, edgeMesh>
        (
            IOobject
            (
                "interpolate("+vf.name()+')',
                vf.instance(),
                vf.db()
            ),
            mesh,
            vf.dimensions()
        )
    );
    GeometricField<Type, faePatchField, edgeMesh>& sf = tsf.ref();

    Field<Type>& sfi = sf.primitiveFieldRef();

    for (label fi=0; fi<P.size(); ++fi)
    {
        const tensorField& curT = mesh.edgeTransformTensors()[fi];

        const tensor& Te = curT[0];
        const tensor& TP = curT[1];
        const tensor& TN = curT[2];

        sfi[fi] =
            transform
            (
                Te.T(),
                lambda[fi]*transform(TP, vfi[P[fi]])
              + y[fi]*transform(TN, vfi[N[fi]])
            );
    }


    // Interpolate across coupled patches using given lambdas and ys

    forAll(lambdas.boundaryField(), pi)
    {
        const faePatchScalarField& pLambda = lambdas.boundaryField()[pi];
        const faePatchScalarField& pY = ys.boundaryField()[pi];

        if (vf.boundaryField()[pi].coupled())
        {
            label size = vf.boundaryField()[pi].patch().size();
            label start = vf.boundaryField()[pi].patch().start();

            Field<Type> pOwnVf(vf.boundaryField()[pi].patchInternalField());
            Field<Type> pNgbVf(vf.boundaryField()[pi].patchNeighbourField());

            Field<Type>& pSf = sf.boundaryFieldRef()[pi];

            for (label i=0; i<size; ++i)
            {
                const tensorField& curT =
                    mesh.edgeTransformTensors()[start + i];

                const tensor& Te = curT[0];
                const tensor& TP = curT[1];
                const tensor& TN = curT[2];

                pSf[i] =
                    transform
                    (
                        Te.T(),
                        pLambda[i]*transform(TP, pOwnVf[i])
                      + pY[i]*transform(TN, pNgbVf[i])
                    );
            }

//             sf.boundaryFieldRef()[pi] =
//                 pLambda*vf.boundaryField()[pi].patchInternalField()
//               + pY*vf.boundaryField()[pi].patchNeighbourField();
        }
        else
        {
            sf.boundaryFieldRef()[pi] = vf.boundaryField()[pi];
        }
    }

    tlambdas.clear();
    tys.clear();

    return tsf;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::faePatchField, Foam::edgeMesh>>
Foam::edgeInterpolationScheme<Type>::interpolate
(
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const tmp<edgeScalarField>& tlambdas
)
{
    if (edgeInterpolation::debug)
    {
        InfoInFunction
            << "interpolating "
            << vf.type() << " "
            << vf.name()
            << " from area to edges "
               "without explicit correction"
            << endl;
    }

    const edgeScalarField& lambdas = tlambdas();

    const Field<Type>& vfi = vf.internalField();
    const scalarField& lambda = lambdas.internalField();

    const faMesh& mesh = vf.mesh();
    const labelUList& P = mesh.owner();
    const labelUList& N = mesh.neighbour();

    tmp<GeometricField<Type, faePatchField, edgeMesh>> tsf
    (
        new GeometricField<Type, faePatchField, edgeMesh>
        (
            IOobject
            (
                "interpolate("+vf.name()+')',
                vf.instance(),
                vf.db()
            ),
            mesh,
            vf.dimensions()
        )
    );
    GeometricField<Type, faePatchField, edgeMesh>& sf = tsf.ref();

    Field<Type>& sfi = sf.primitiveFieldRef();

    for (label eI = 0; eI < P.size(); ++eI)
    {
        const tensorField& curT = mesh.edgeTransformTensors()[eI];

        const tensor& Te = curT[0];
        const tensor& TP = curT[1];
        const tensor& TN = curT[2];

        sfi[eI] =
            transform
            (
                Te.T(),
                lambda[eI]*transform(TP, vfi[P[eI]])
              + (1 - lambda[eI])*transform(TN, vfi[N[eI]])
            );
    }


    // Interpolate across coupled patches using given lambdas

    const auto& vfb = vf.boundaryField();

    forAll(lambdas.boundaryField(), pi)
    {
        const faePatchScalarField& pLambda = lambdas.boundaryField()[pi];

        if (vf.boundaryField()[pi].coupled())
        {
            label size = vfb[pi].patch().size();
            label start = vfb[pi].patch().start();

            Field<Type> pOwnVf(vfb[pi].patchInternalField());
            Field<Type> pNgbVf(vfb[pi].patchNeighbourField());

            Field<Type>& pSf = sf.boundaryFieldRef()[pi];

            for (label i=0; i<size; ++i)
            {
                const tensorField& curT =
                    mesh.edgeTransformTensors()[start + i];

                const tensor& Te = curT[0];
                const tensor& TP = curT[1];
                const tensor& TN = curT[2];

                pSf[i] =
                    transform
                    (
                        Te.T(),
                        pLambda[i]*transform(TP, pOwnVf[i])
                      + (1 - pLambda[i])*transform(TN, pNgbVf[i])
                    );
            }

//             tsf().boundaryFieldRef()[pi] =
//                 pLambda*vf.boundaryField()[pi].patchInternalField()
//              + (1 - pLambda)*vf.boundaryField()[pi].patchNeighbourField();
        }
        else
        {
            sf.boundaryFieldRef()[pi] = vf.boundaryField()[pi];
        }
    }

    tlambdas.clear();

    return tsf;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::faePatchField, Foam::edgeMesh>>
Foam::edgeInterpolationScheme<Type>::euclidianInterpolate
(
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const tmp<edgeScalarField>& tlambdas
)
{
    if (edgeInterpolation::debug)
    {
        InfoInFunction
            << "interpolating "
            << vf.type() << " "
            << vf.name()
            << " from area to edges "
               "without explicit correction"
            << endl;
    }

    const edgeScalarField& lambdas = tlambdas();

    const Field<Type>& vfi = vf.internalField();
    const scalarField& lambda = lambdas.internalField();

    const faMesh& mesh = vf.mesh();
    const labelUList& P = mesh.owner();
    const labelUList& N = mesh.neighbour();

    tmp<GeometricField<Type, faePatchField, edgeMesh>> tsf
    (
        new GeometricField<Type, faePatchField, edgeMesh>
        (
            IOobject
            (
                "interpolate("+vf.name()+')',
                vf.instance(),
                vf.db()
            ),
            mesh,
            vf.dimensions()
        )
    );
    GeometricField<Type, faePatchField, edgeMesh>& sf = tsf.ref();

    Field<Type>& sfi = sf.primitiveFieldRef();

    for (label eI = 0; eI < P.size(); ++eI)
    {
        sfi[eI] = lambda[eI]*vfi[P[eI]] + (1 - lambda[eI])*vfi[N[eI]];
    }


    // Interpolate across coupled patches using given lambdas

    forAll(lambdas.boundaryField(), pi)
    {
        const faePatchScalarField& pLambda = lambdas.boundaryField()[pi];

        if (vf.boundaryField()[pi].coupled())
        {
            tsf.ref().boundaryFieldRef()[pi] =
                pLambda*vf.boundaryField()[pi].patchInternalField()
             + (1.0 - pLambda)*vf.boundaryField()[pi].patchNeighbourField();
        }
        else
        {
            sf.boundaryFieldRef()[pi] = vf.boundaryField()[pi];
        }
    }

    tlambdas.clear();

    return tsf;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::faePatchField, Foam::edgeMesh>>
Foam::edgeInterpolationScheme<Type>::interpolate
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
) const
{
    if (edgeInterpolation::debug)
    {
        InfoInFunction
            << "interpolating "
            << vf.type() << " "
            << vf.name()
            << " from areas to edges"
            << endl;
    }

    tmp<GeometricField<Type, faePatchField, edgeMesh>> tsf =
        interpolate(vf, weights(vf));

    if (corrected())
    {
        tsf.ref() += correction(vf);
    }

    return tsf;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::faePatchField, Foam::edgeMesh>>
Foam::edgeInterpolationScheme<Type>::euclidianInterpolate
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
) const
{
    if (edgeInterpolation::debug)
    {
        InfoInFunction
            << "interpolating "
            << vf.type() << " "
            << vf.name()
            << " from area to edges "
            << endl;
    }

    tmp<GeometricField<Type, faePatchField, edgeMesh>> tsf =
        euclidianInterpolate(vf, weights(vf));

    return tsf;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::faePatchField, Foam::edgeMesh>>
Foam::edgeInterpolationScheme<Type>::interpolate
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tvf
) const
{
    tmp<GeometricField<Type, faePatchField, edgeMesh>> tinterpVf =
        interpolate(tvf());
    tvf.clear();
    return tinterpVf;
}


// ************************************************************************* //
