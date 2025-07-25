/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

#include "AveragingMethod.H"
#include "runTimeSelectionTables.H"
#include "pointMesh.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::AveragingMethod<Type>::updateGrad()
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::AveragingMethod<Type>::AveragingMethod
(
    const IOobject& io,
    const dictionary& dict,
    const fvMesh& mesh,
    const labelList& size
)
:
    regIOobject(io),
    FieldField<Field, Type>(),
    dict_(dict),
    mesh_(mesh)
{
    forAll(size, i)
    {
        FieldField<Field, Type>::append
        (
            new Field<Type>(size[i], Zero)
        );
    }
}


template<class Type>
Foam::AveragingMethod<Type>::AveragingMethod
(
    const AveragingMethod<Type>& am
)
:
    regIOobject(am),
    FieldField<Field, Type>(am),
    dict_(am.dict_),
    mesh_(am.mesh_)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::AveragingMethod<Type>>
Foam::AveragingMethod<Type>::New
(
    const IOobject& io,
    const dictionary& dict,
    const fvMesh& mesh
)
{
    const word modelType
    (
        dict.template getOrDefault<word>(typeName, "basic")
    );

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "averaging limiter",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << abort(FatalIOError);
    }

    return autoPtr<AveragingMethod<Type>>(ctorPtr(io, dict, mesh));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::AveragingMethod<Type>::average()
{
    updateGrad();
}


template<class Type>
void Foam::AveragingMethod<Type>::average
(
    const AveragingMethod<scalar>& weight
)
{
    updateGrad();

    *this /= max(weight, SMALL);
}


template<class Type>
bool Foam::AveragingMethod<Type>::writeData(Ostream& os) const
{
    return os.good();
}


template<class Type>
bool Foam::AveragingMethod<Type>::write(const bool writeOnProc) const
{
    const pointMesh pointMesh_(mesh_);

    // point volumes
    Field<scalar> pointVolume(mesh_.nPoints(), Zero);

    // output fields
    auto tcellValue = GeometricField<Type, fvPatchField, volMesh>::New
    (
        IOobject::scopedName(this->name(), "cellValue"),
        IOobject::NO_REGISTER,
        mesh_,
        dimensioned<Type>(dimless, Zero)
    );
    auto& cellValue = tcellValue.ref();

    auto tcellGrad = GeometricField<TypeGrad, fvPatchField, volMesh>::New
    (
        IOobject::scopedName(this->name(), "cellGrad"),
        IOobject::NO_REGISTER,
        mesh_,
        dimensioned<TypeGrad>(dimless, Zero)
    );
    auto& cellGrad = tcellGrad.ref();

    auto tpointValue = GeometricField<Type, pointPatchField, pointMesh>::New
    (
        IOobject::scopedName(this->name(), "pointValue"),
        IOobject::NO_REGISTER,
        pointMesh_,
        dimensioned<Type>(dimless, Zero)
    );
    auto& pointValue = tpointValue.ref();

    auto tpointGrad = GeometricField<TypeGrad, pointPatchField, pointMesh>::New
    (
        IOobject::scopedName(this->name(), "pointGrad"),
        IOobject::NO_REGISTER,
        pointMesh_,
        dimensioned<TypeGrad>(dimless, Zero)
    );
    auto& pointGrad = tpointGrad.ref();

    // Barycentric coordinates of the tet vertices
    const FixedList<barycentric, 4>
        tetCrds
        ({
            barycentric(1, 0, 0, 0),
            barycentric(0, 1, 0, 0),
            barycentric(0, 0, 1, 0),
            barycentric(0, 0, 0, 1)
        });

    // tet-volume weighted sums
    forAll(mesh_.C(), celli)
    {
        const List<tetIndices> cellTets =
            polyMeshTetDecomposition::cellTetIndices(mesh_, celli);

        forAll(cellTets, tetI)
        {
            const tetIndices& tetIs = cellTets[tetI];
            const triFace triIs = tetIs.faceTriIs(mesh_);
            const scalar v = tetIs.tet(mesh_).mag();

            cellValue[celli] += v*interpolate(tetCrds[0], tetIs);
            cellGrad[celli] += v*interpolateGrad(tetCrds[0], tetIs);

            forAll(triIs, vertexI)
            {
                const label pointi = triIs[vertexI];

                pointVolume[pointi] += v;
                pointValue[pointi] += v*interpolate(tetCrds[vertexI], tetIs);
                pointGrad[pointi] += v*interpolateGrad(tetCrds[vertexI], tetIs);
            }
        }
    }

    // average
    cellValue.primitiveFieldRef() /= mesh_.V();
    cellGrad.primitiveFieldRef() /= mesh_.V();
    pointValue.primitiveFieldRef() /= pointVolume;
    pointGrad.primitiveFieldRef() /= pointVolume;

    // write
    if (!cellValue.write(writeOnProc)) return false;
    if (!cellGrad.write(writeOnProc)) return false;
    if (!pointValue.write(writeOnProc)) return false;
    if (!pointGrad.write(writeOnProc)) return false;

    return true;
}


// ************************************************************************* //
