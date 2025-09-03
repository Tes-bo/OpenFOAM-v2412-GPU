/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 Zeljko Tukovic, FSB Zagreb.
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "skewCorrectedSnGrad.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "linear.H"
#include "fvcGrad.H"
#include "gaussGrad.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::fv::skewCorrectedSnGrad<Type>::fullGradCorrection
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    auto tssf = tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>::New
    (
        IOobject
        (
            "snGradCorr("+vf.name()+')',
            vf.instance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        vf.dimensions()*mesh.nonOrthDeltaCoeffs().dimensions()
    );
    auto& ssf = tssf.ref();

    ssf.setOriented();
    ssf = dimensioned<Type>(ssf.dimensions(), Zero);


    typedef typename
        outerProduct<vector, typename pTraits<Type>::cmptType>::type
        CmptGradType;

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    const vectorField& Sf = mesh.Sf().internalField();
    const scalarField& magSf = mesh.magSf().internalField();

    const vectorField& Cf = mesh.Cf().internalField();
    const vectorField& C = mesh.C().internalField();

    const scalarField& deltaCoeffs =
        mesh.deltaCoeffs().internalField();

    surfaceVectorField kP
    (
        IOobject
        (
            "kP",
            vf.instance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector(dimLength, Zero)
    );
    vectorField& kPI = kP.primitiveFieldRef();

    surfaceVectorField kN
    (
        IOobject
        (
            "kN",
            vf.instance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector(dimLength, Zero)
    );
    vectorField& kNI = kN.primitiveFieldRef();

    kPI = Cf - vectorField(C, owner);
    kPI -= Sf*(Sf&kPI)/sqr(magSf);

    kNI = Cf - vectorField(C, neighbour);
    kNI -= Sf*(Sf&kNI)/sqr(magSf);

    forAll(kP.boundaryField(), patchI)
    {
        if (kP.boundaryField()[patchI].coupled())
        {
            kP.boundaryFieldRef()[patchI] =
                mesh.boundary()[patchI].Cf()
              - mesh.boundary()[patchI].Cn();

            kP.boundaryFieldRef()[patchI] -=
                mesh.boundary()[patchI].Sf()
               *(
                    mesh.boundary()[patchI].Sf()
                  & kP.boundaryField()[patchI]
                )
               /sqr(mesh.boundary()[patchI].magSf());

            kN.boundaryFieldRef()[patchI] =
                mesh.Cf().boundaryField()[patchI]
              - (
                    mesh.boundary()[patchI].Cn()
                  + mesh.boundary()[patchI].delta()
                );

            kN.boundaryFieldRef()[patchI] -=
                mesh.boundary()[patchI].Sf()
               *(
                    mesh.boundary()[patchI].Sf()
                  & kN.boundaryField()[patchI]
                )
               /sqr(mesh.boundary()[patchI].magSf());
        }
    }

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; ++cmpt)
    {
        GeometricField<CmptGradType, fvPatchField, volMesh> cmptGrad
        (
            gradScheme<typename pTraits<Type>::cmptType>::New
            (
                mesh,
                mesh.gradScheme("grad(" + vf.name() + ')')
            )()
           .grad(vf.component(cmpt))
        );

        const Field<CmptGradType>& cmptGradI = cmptGrad.internalField();

        // Skewness and non-rothogonal correction
        {
            ssf.primitiveFieldRef().replace
            (
                cmpt,
                (
                    (kNI&Field<CmptGradType>(cmptGradI, neighbour))
                  - (kPI&Field<CmptGradType>(cmptGradI, owner))
                )
               *deltaCoeffs
            );
        }

        forAll(ssf.boundaryField(), patchI)
        {
            if (ssf.boundaryField()[patchI].coupled())
            {
                ssf.boundaryFieldRef()[patchI].replace
                (
                    cmpt,
                    (
                        (
                            kN.boundaryField()[patchI]
                          & cmptGrad.boundaryField()[patchI]
                           .patchNeighbourField()
                        )
                      - (
                            kP.boundaryField()[patchI]
                          & cmptGrad.boundaryField()[patchI]
                           .patchInternalField()
                        )
                    )
                   *mesh.deltaCoeffs().boundaryField()[patchI]
                );
            }
        }
    }

    // // construct GeometricField<Type, fvsPatchField, surfaceMesh>
    // tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tssf =
    //     linear<typename outerProduct<vector, Type>::type>(mesh).dotInterpolate
    //     (
    //         mesh.nonOrthCorrectionVectors(),
    //         gradScheme<Type>::New
    //         (
    //             mesh,
    //             mesh.gradScheme("grad(" + vf.name() + ')')
    //         )().grad(vf, "grad(" + vf.name() + ')')
    //     );
    //
    // tssf.ref().rename("snGradCorr(" + vf.name() + ')');

    return tssf;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::fv::skewCorrectedSnGrad<Type>::correction
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    auto tssf = tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>::New
    (
        IOobject
        (
            "snGradCorr("+vf.name()+')',
            vf.instance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        vf.dimensions()*mesh.nonOrthDeltaCoeffs().dimensions()
    );
    auto& ssf = tssf.ref();
    ssf.setOriented();

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; ++cmpt)
    {
        ssf.replace
        (
            cmpt,
            skewCorrectedSnGrad<typename pTraits<Type>::cmptType>(mesh)
           .fullGradCorrection(vf.component(cmpt))
        );
    }

    return tssf;
}


// ************************************************************************* //
