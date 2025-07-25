/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
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

#include "shapeDesignVariables.H"
#include "cellQuality.H"
#include "createZeroField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFieldsFwd.H"
#include "adjointEikonalSolver.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(shapeDesignVariables, 0);
    defineRunTimeSelectionTable(shapeDesignVariables, dictionary);
    addToRunTimeSelectionTable
    (
        designVariables,
        shapeDesignVariables,
        designVariables
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::label Foam::shapeDesignVariables::sensSize() const
{
    return size();
}


const Foam::labelList& Foam::shapeDesignVariables::activeSensitivities() const
{
    return activeDesignVariables_;
}


Foam::tmp<Foam::volVectorField>
Foam::shapeDesignVariables::solveMeshMovementEqn
(
    const label patchI,
    const label varID
) const
{
    const dictionary dxdbDict = dict_.subOrEmptyDict("dxdbSolver");
    const label iters = dxdbDict.getOrDefault<label>("iters", 1000);
    const scalar tolerance =
        dxdbDict.getOrDefault<scalar>("tolerance", 1.e-07);
    tmp<volVectorField> tm
    (
        tmp<volVectorField>::New
        (
            variablesSet::autoCreateMeshMovementField
            (
                mesh_,
                "m",
                dimensionSet(dimLength)
            )
        )
    );
    volVectorField& m = tm.ref();

    // Solve for dxdb
    //~~~~~~~~~~~~~~~~
    m.boundaryFieldRef()[patchI] == dxdbFace(patchI, varID);

    // Iterate the direct differentiation of the grid displacement  equation
    for (label iter = 0; iter < iters; ++iter)
    {
        Info<< "Mesh Movement Propagation for varID" << varID
            << ", Iteration : "<< iter << endl;

        fvVectorMatrix mEqn
        (
            fvm::laplacian(m)
        );

        scalar residual = mag(mEqn.solve().initialResidual());

        DebugInfo
            << "Max dxdb " << gMax(mag(m)()) << endl;

        mesh_.time().printExecutionTime(Info);

        // Check convergence
        if (residual < tolerance)
        {
            Info<< "\n***Reached dxdb convergence limit, iteration " << iter
                << "***\n\n";
            break;
        }
    }

    return tm;
}


void Foam::shapeDesignVariables::allocateSensFields()
{
    if (dxdbVolSens_.empty())
    {
        dxdbVolSens_.setSize(sensSize(), Zero);
        dxdbSurfSens_.setSize(sensSize(), Zero);
        dSdbSens_.setSize(sensSize(), Zero);
        dndbSens_.setSize(sensSize(), Zero);
        dxdbDirectSens_.setSize(sensSize(), Zero);
        dVdbSens_.setSize(sensSize(), Zero);
        distanceSens_.setSize(sensSize(), Zero);
        optionsSens_.setSize(sensSize(), Zero);
        bcSens_.setSize(sensSize(), Zero);
    }
}


void Foam::shapeDesignVariables::zeroSensFields()
{
    dxdbVolSens_ = Zero;
    dxdbSurfSens_ = Zero;
    dSdbSens_ = Zero;
    dndbSens_ = Zero;
    dxdbDirectSens_ = Zero;
    dVdbSens_ = Zero;
    distanceSens_ = Zero;
    optionsSens_ = Zero;
    bcSens_ = Zero;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::shapeDesignVariables::shapeDesignVariables
(
    fvMesh& mesh,
    const dictionary& dict
)
:
    designVariables(mesh, dict),
    parametertisedPatches_
    (
        mesh_.boundaryMesh().patchSet(dict.get<wordRes>("patches"))
    ),
    displMethodPtr_
    (
        displacementMethod::New(mesh_, parametertisedPatches_.toc())
    ),
    pointsInit_(nullptr),
    writeEachMesh_(dict.getOrDefault<bool>("writeEachMesh", true)),
    dxdbVolSens_(),
    dxdbSurfSens_(),
    dSdbSens_(),
    dndbSens_(),
    dxdbDirectSens_(),
    dVdbSens_(),
    distanceSens_(),
    optionsSens_(),
    bcSens_(),
    derivativesFolder_
    (
        word("optimisation")/word("derivatives")
       /word(mesh.name() == polyMesh::defaultRegion ? word() : mesh.name())
    )
{
    if (!parametertisedPatches_.size())
    {
        FatalErrorInFunction
            << "None of the provided parameterised patches is valid"
            << endl
            << exit(FatalError);
    }
    mkDir(derivativesFolder_);
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::shapeDesignVariables> Foam::shapeDesignVariables::New
(
    fvMesh& mesh,
    const dictionary& dict
)
{
    const word modelType(dict.get<word>("shapeType"));

    Info<< "shapeDesignVariables type : " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalErrorInLookup
        (
            "shapeType",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalError);
    }

    return autoPtr<shapeDesignVariables>(ctorPtr(mesh, dict));
}


// * * * * * * * * * * * * * * *  Member Functions   * * * * * * * * * * * * //

bool Foam::shapeDesignVariables::readDict(const dictionary& dict)
{
    if (designVariables::readDict(dict))
    {
        parametertisedPatches_ =
            mesh_.boundaryMesh().patchSet(dict.get<wordRes>("patches"));
        displMethodPtr_->setPatchIDs(parametertisedPatches_.toc());

        writeEachMesh_ =
            dict.getOrDefault<bool>("writeEachMesh", true);

        return true;
    }

    return false;
}


void Foam::shapeDesignVariables::storeDesignVariables()
{
    designVariables::storeDesignVariables();

    if (!pointsInit_)
    {
        pointsInit_.reset(new pointField(mesh_.nPoints(), Zero));
    }
    pointsInit_() = mesh_.points();
}


void Foam::shapeDesignVariables::resetDesignVariables()
{
    designVariables::resetDesignVariables();
    mesh_.movePoints(pointsInit_());
}


void Foam::shapeDesignVariables::moveMesh()
{
    // Move mesh
    displMethodPtr_->update();

    if (writeEachMesh_)
    {
        Info<< "  Writing new mesh points for mesh region "
            << mesh_.name() << endl;
        pointIOField points
        (
            IOobject
            (
               "points",
                mesh_.pointsInstance(),
                polyMesh::meshSubDir,
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh_.points()
        );
        points.write();
    }

    // Check mesh quality
    mesh_.checkMesh(true);
}


Foam::tmp<Foam::scalarField> Foam::shapeDesignVariables::assembleSensitivities
(
    adjointSensitivity& adjointSens
)
{
    // Return field
    tmp<scalarField> tsens(tmp<scalarField>::New(sensSize(), Zero));
    scalarField& sens = tsens.ref();

    // Reset sensitivity components to zero
    allocateSensFields();
    zeroSensFields();

    // Grab multipliers from the adjoint sensitivities
    const autoPtr<volTensorField>& gradDxDbMult = adjointSens.gradDxDbMult();
    const autoPtr<scalarField>& divDxDbMult = adjointSens.divDxDbMult();
    const autoPtr<boundaryVectorField>& dxdbMult = adjointSens.dxdbMult();
    const autoPtr<boundaryVectorField>& dSdbMult = adjointSens.dSfdbMult();
    const autoPtr<boundaryVectorField>& dndbMult = adjointSens.dnfdbMult();
    const autoPtr<boundaryVectorField>& dxdbDirectMult =
        adjointSens.dxdbDirectMult();
    const autoPtr<boundaryVectorField>& bcDxDbmult = adjointSens.bcDxDbMult();
    const autoPtr<vectorField>& optionsDxDbMult = adjointSens.optionsDxDbMult();
    const volScalarField::Internal& V = mesh_.V();
    autoPtr<adjointEikonalSolver>& eikonalSolver =
        adjointSens.getAdjointEikonalSolver();

    autoPtr<volTensorField> distanceSens(nullptr);
    if (adjointSens.includeDistance())
    {
        distanceSens.reset
        (
            new volTensorField(eikonalSolver->getFISensitivityTerm())
        );
    }

    // Loop over active design variables only
    for (const label varI : activeSensitivities())
    {
        // FI approach, integrate terms including variations of the grid
        // sensitivities
        if (adjointSens.computeDxDbInternalField())
        {
            // Parameterization info
            tmp<volVectorField> tvolDxDbI = dCdb(varI);
            const volVectorField& volDxDbI = tvolDxDbI();
            tmp<volTensorField> gradDxDb = fvc::grad(volDxDbI);

            // Contributions from the adjoint-related part
            dxdbVolSens_[varI] = gSum((gradDxDbMult() && gradDxDb())*V);

            // Contributions from the distance related part
            if (adjointSens.includeDistance())
            {
                distanceSens_[varI] = gSum((distanceSens() && gradDxDb)*V);
            }
            // Contributions from the multiplier of divDxDb
            if (divDxDbMult)
            {
                dVdbSens_[varI] +=
                    gSum(divDxDbMult()*fvc::div(volDxDbI)().primitiveField()*V);
            }

            // Contributions from fvOptions
            optionsSens_[varI] +=
                gSum((optionsDxDbMult() & volDxDbI.primitiveField())*V);
        }

        // Contribution from boundary terms
        // Most of them (with the expection of dxdbMult) exist in both the
        // FI and E-SI approaches
        for (const label patchI : parametertisedPatches_)
        {
            if (dSdbMult)
            {
                tmp<vectorField> pdSdb = dSdb(patchI, varI);
                dSdbSens_[varI] += gSum(dSdbMult()[patchI] & pdSdb);
            }

            if (dndbMult)
            {
                tmp<vectorField> pdndb = dndb(patchI, varI);
                dndbSens_[varI] += gSum((dndbMult()[patchI] & pdndb));
            }

            tmp<vectorField> pdxdb = dxdbFace(patchI, varI);
            // Main contribution in the E-SI approach
            if (dxdbMult)
            {
                dxdbSurfSens_[varI] += gSum(dxdbMult()[patchI] & pdxdb());
            }
            if (dxdbDirectMult)
            {
                dxdbDirectSens_[varI] +=
                    gSum((dxdbDirectMult()[patchI] & pdxdb()));
            }
            if (bcDxDbmult)
            {
                bcSens_[varI] += gSum((bcDxDbmult()[patchI] & pdxdb()));
            }
        }
    }

    sens =
        dxdbVolSens_ + dxdbSurfSens_ + dSdbSens_ + dndbSens_ + dxdbDirectSens_
      + dVdbSens_ + distanceSens_ + optionsSens_ + bcSens_;

    writeSensitivities(sens, adjointSens);

    return tsens;
}


void Foam::shapeDesignVariables::writeSensitivities
(
    const scalarField& sens,
    const adjointSensitivity& adjointSens
)
{
    OFstream derivFile
    (
        derivativesFolder_/
            type() + adjointSens.getAdjointSolver().solverName()
          + adjointSens.getSuffix() + mesh_.time().timeName()
    );
    unsigned int widthDV = max(int(name(dxdbVolSens_.size()).size()), int(6));
    unsigned int width = IOstream::defaultPrecision() + 7;
    derivFile
        << setw(widthDV) << "#varID" << " "
        << setw(width) << "total"<< " "
        << setw(width) << "dxdbVol" << " "
        << setw(width) << "dxdbSurf" << " "
        << setw(width) << "dSdb" << " "
        << setw(width) << "dndb" << " "
        << setw(width) << "dxdbDirect" << " "
        << setw(width) << "dVdb" << " "
        << setw(width) << "distance" << " "
        << setw(width) << "options" << " "
        << setw(width) << "dvdb" << endl;

    for (const label varI : activeSensitivities())
    {
        derivFile
           << setw(widthDV) << varI << " "
           << setw(width) << sens[varI] << " "
           << setw(width) << dxdbVolSens_[varI] << " "
           << setw(width) << dxdbSurfSens_[varI] << " "
           << setw(width) << dSdbSens_[varI] << " "
           << setw(width) << dndbSens_[varI] << " "
           << setw(width) << dxdbDirectSens_[varI] << " "
           << setw(width) << dVdbSens_[varI] << " "
           << setw(width) << distanceSens_[varI] << " "
           << setw(width) << optionsSens_[varI] << " "
           << setw(width) << bcSens_[varI] << " "
           << endl;
    }
}


Foam::tmp<Foam::vectorField> Foam::shapeDesignVariables::dxdbVol
(
    const label varID
) const
{
    // Deliberately returning a zero-sized field
    return tmp<vectorField>::New(0);
}


Foam::tmp<Foam::vectorField> Foam::shapeDesignVariables::dxdbFace
(
    const label patchI,
    const label varID
) const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::vectorField> Foam::shapeDesignVariables::dndb
(
    const label patchI,
    const label varID
) const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::vectorField> Foam::shapeDesignVariables::dSdb
(
    const label patchI,
    const label varID
) const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::volVectorField>
Foam::shapeDesignVariables::dCdb(const label varID) const
{
    NotImplemented;
    return nullptr;
}


// ************************************************************************* //
