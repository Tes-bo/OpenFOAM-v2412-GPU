/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2021 PCOpt/NTUA
    Copyright (C) 2013-2021 FOSS GP
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "adjointBoundaryCondition.H"
#include "emptyFvPatch.H"
#include "adjointSolverManager.H"
#include "HashTable.H"
#include "surfaceInterpolationScheme.H"
#include "ATCUaGradU.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
template<class Type2>
tmp
<
    Field<typename Foam::outerProduct<Foam::vector, Type2>::type>
>
adjointBoundaryCondition<Type>::computePatchGrad(word name)
{
    // Return field
    typedef typename outerProduct<vector, Type2>::type GradType;
    auto tresGrad = tmp<Field<GradType>>::New(patch_.size(), Zero);
    auto& resGrad = tresGrad.ref();

    const labelList& faceCells = patch_.faceCells();
    const fvMesh& mesh = patch_.boundaryMesh().mesh();
    const cellList& cells = mesh.cells();

    // Go through the surfaceInterpolation scheme defined in gradSchemes for
    // consistency
    const GeometricField<Type2, fvPatchField, volMesh>& field =
        mesh.lookupObject<volVectorField>(name);

    // Gives problems when grad(AdjointVar) is computed using a limited scheme,
    // since it is not possible to know a priori how many words to expect in the
    // stream.
    // Interpolation scheme is now read through interpolation schemes.
    /*
    word gradSchemeName("grad(" + name + ')');
    ITstream& is = mesh.gradScheme(gradSchemeName);
    word schemeData(is);
    */

    // If there are many patches calling this function, the computation of
    // the surface field might be a significant computational
    // burden. Cache the interpolated field and fetch from the registry when
    // possible.
    const word surfFieldName("interpolated" + name + "ForBoundaryGrad");
    typedef GeometricField<Type2, fvsPatchField, surfaceMesh> surfFieldType;
    surfFieldType* surfFieldPtr =
        mesh.objectRegistry::template getObjectPtr<surfFieldType>(surfFieldName);

    if (!surfFieldPtr)
    {
        solution::cachePrintMessage("Calculating and caching", name, field);

        surfFieldPtr =
            surfaceInterpolationScheme<Type2>::New
            (
                mesh, mesh.interpolationScheme("interpolate(" + name + ")")
            )().interpolate(field).ptr();
        surfFieldPtr->rename(surfFieldName);
        regIOobject::store(surfFieldPtr);
    }
    else
    {
        if (surfFieldPtr->upToDate(field))
        {
            solution::cachePrintMessage("Reusing", name, field);
        }
        else
        {
            solution::cachePrintMessage("Updating", name, field);
            delete surfFieldPtr;

            surfFieldPtr =
                surfaceInterpolationScheme<Type2>::New
                (
                    mesh, mesh.interpolationScheme("interpolate(" + name + ")")
                )().interpolate(field).ptr();
            surfFieldPtr->rename(surfFieldName);
            regIOobject::store(surfFieldPtr);
        }
    }
    surfFieldType& surfField = *surfFieldPtr;

    // Auxiliary fields
    const surfaceVectorField& Sf = mesh.Sf();
    tmp<vectorField> tnf = patch_.nf();
    const vectorField& nf = tnf();
    const scalarField& V = mesh.V();
    const labelUList& owner = mesh.owner();

    // Compute grad value of cell adjacent to the boundary
    forAll(faceCells, fI)
    {
        const label cI = faceCells[fI];
        const cell& cellI = cells[cI];
        for (const label faceI : cellI) // global face numbering
        {
            label patchID = mesh.boundaryMesh().whichPatch(faceI);
            if (patchID == -1) //face is internal
            {
                const label own = owner[faceI];
                tensor flux = Sf[faceI]*surfField[faceI];
                if (cI == own)
                {
                    resGrad[fI] += flux;
                }
                else
                {
                    resGrad[fI] -= flux;
                }
            }
            else  // Face is boundary. Covers coupled patches as well
            {
                if (!isA<emptyFvPatch>(mesh.boundary()[patchID]))
                {
                    const fvPatch& patchForFlux = mesh.boundary()[patchID];
                    const label boundaryFaceI = faceI - patchForFlux.start();
                    const vectorField& Sfb = Sf.boundaryField()[patchID];
                    resGrad[fI] +=
                        Sfb[boundaryFaceI]
                       *surfField.boundaryField()[patchID][boundaryFaceI];
                }
            }
        }
        resGrad[fI] /= V[cI];
    }

    // This has concluded the computation of the grad at the cell next to the
    // boundary. We now need to compute the grad at the boundary face
    const fvPatchField<Type2>& bField = field.boundaryField()[patch_.index()];
    resGrad = nf*bField.snGrad() + (resGrad - nf*(nf & resGrad));

    return tresGrad;
}


template<class Type>
bool adjointBoundaryCondition<Type>::addATCUaGradUTerm()
{
    if (!addATCUaGradUTerm_.good())
    {
        addATCUaGradUTerm_ = bool(isA<ATCUaGradU>(getATC()));
    }
    return addATCUaGradUTerm_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
adjointBoundaryCondition<Type>::adjointBoundaryCondition
(
    const adjointBoundaryCondition<Type>& adjointBC
)
:
    patch_(adjointBC.patch_),
    managerName_(adjointBC.managerName_),
    adjointSolverName_(adjointBC.adjointSolverName_),
    simulationType_(adjointBC.simulationType_),
    boundaryContrPtr_
    (
        boundaryAdjointContribution::New
        (
            adjointBC.managerName_,
            adjointBC.adjointSolverName_,
            adjointBC.simulationType_,
            adjointBC.patch_
        )
    ),
    addATCUaGradUTerm_(Switch::INVALID)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
adjointBoundaryCondition<Type>::adjointBoundaryCondition
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const word& solverName
)
:
    patch_(p),
    managerName_("objectiveManager" + solverName),
    adjointSolverName_(solverName),
    simulationType_("incompressible"),
    boundaryContrPtr_(nullptr),
    addATCUaGradUTerm_(Switch::INVALID)
{
    // Set the boundaryContribution pointer
    setBoundaryContributionPtr();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const word& adjointBoundaryCondition<Type>::objectiveManagerName() const
{
    return managerName_;
}


template<class Type>
const word& adjointBoundaryCondition<Type>::adjointSolverName() const
{
    return adjointSolverName_;
}


template<class Type>
const word& adjointBoundaryCondition<Type>::simulationType() const
{
    return simulationType_;
}


template<class Type>
void adjointBoundaryCondition<Type>::setBoundaryContributionPtr()
{
    // Note:
    // Check whether there is an objectiveFunctionManager object in the registry
    // Necessary for decomposePar if the libadjoint is loaded
    // through controlDict. A nicer way should be found
    const fvMesh& meshRef = patch_.boundaryMesh().mesh();
    if (meshRef.foundObject<regIOobject>(managerName_))
    {
        boundaryContrPtr_.reset
        (
            boundaryAdjointContribution::New
            (
                managerName_,
                adjointSolverName_,
                simulationType_,
                patch_
            ).ptr()
        );
    }
    else
    {
        WarningInFunction
            << "No objectiveManager " << managerName_ << " available." << nl
            << "Setting boundaryAdjointContributionPtr to nullptr. " << nl
            << "OK for decomposePar."
            << endl;
    }
}


template<class Type>
boundaryAdjointContribution&
adjointBoundaryCondition<Type>::getBoundaryAdjContribution()
{
    return boundaryContrPtr_();
}


template<class Type>
const ATCModel& adjointBoundaryCondition<Type>::getATC() const
{
    return
        patch_.boundaryMesh().mesh().template
            lookupObject<ATCModel>("ATCModel" + adjointSolverName_);
}


template<class Type>
void adjointBoundaryCondition<Type>::updatePrimalBasedQuantities()
{
    // Does nothing in base
}


template<class Type>
tmp
<
    Field<typename Foam::outerProduct<Foam::vector, Type>::type>
>
adjointBoundaryCondition<Type>::dxdbMult() const
{
    return nullptr;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
