/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2023 OpenCFD Ltd.
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

#include "rotorDiskSource.H"
#include "addToRunTimeSelectionTable.H"
#include "trimModel.H"
#include "fvMatrices.H"
#include "geometricOneField.H"
#include "syncTools.H"
#include "unitConversion.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(rotorDiskSource, 0);
        addToRunTimeSelectionTable(option, rotorDiskSource, dictionary);
    }
}


const Foam::Enum
<
    Foam::fv::rotorDiskSource::geometryModeType
>
Foam::fv::rotorDiskSource::geometryModeTypeNames_
({
    { geometryModeType::gmAuto, "auto" },
    { geometryModeType::gmSpecified, "specified" },
});


const Foam::Enum
<
    Foam::fv::rotorDiskSource::inletFlowType
>
Foam::fv::rotorDiskSource::inletFlowTypeNames_
({
    { inletFlowType::ifFixed, "fixed" },
    { inletFlowType::ifSurfaceNormal, "surfaceNormal" },
    { inletFlowType::ifLocal, "local" },
});


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fv::rotorDiskSource::checkData()
{
    // Set inflow type
    switch (selectionMode())
    {
        case smCellSet:
        case smCellZone:
        case smAll:
        {
            // Set the profile ID for each blade section
            profiles_.connectBlades(blade_.profileName(), blade_.profileID());
            switch (inletFlow_)
            {
                case ifFixed:
                {
                    coeffs_.readEntry("inletVelocity", inletVelocity_);
                    break;
                }
                case ifSurfaceNormal:
                {
                    scalar UIn(coeffs_.get<scalar>("inletNormalVelocity"));
                    inletVelocity_ = -coordSys_.e3()*UIn;
                    break;
                }
                case ifLocal:
                {
                    break;
                }
                default:
                {
                    FatalErrorInFunction
                        << "Unknown inlet velocity type" << abort(FatalError);
                }
            }

            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Source cannot be used with '"
                << selectionModeTypeNames_[selectionMode()]
                << "' mode.  Please use one of: " << nl
                << selectionModeTypeNames_[smCellSet] << nl
                << selectionModeTypeNames_[smCellZone] << nl
                << selectionModeTypeNames_[smAll]
                << exit(FatalError);
        }
    }
}


void Foam::fv::rotorDiskSource::setFaceArea(vector& axis, const bool correct)
{
    area_ = 0.0;

    static const scalar tol = 0.8;

    const label nInternalFaces = mesh_.nInternalFaces();
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    const vectorField& Sf = mesh_.Sf();
    const scalarField& magSf = mesh_.magSf();

    vector n = Zero;

    // Calculate cell addressing for selected cells
    labelList cellAddr(mesh_.nCells(), -1);
    labelUIndList(cellAddr, cells_) = identity(cells_.size());
    labelList nbrFaceCellAddr(mesh_.nBoundaryFaces(), -1);
    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];

        if (pp.coupled())
        {
            forAll(pp, i)
            {
                label facei = pp.start() + i;
                label nbrFacei = facei - nInternalFaces;
                label own = mesh_.faceOwner()[facei];
                nbrFaceCellAddr[nbrFacei] = cellAddr[own];
            }
        }
    }

    // Correct for parallel running
    syncTools::swapBoundaryFaceList(mesh_, nbrFaceCellAddr);

    // Add internal field contributions
    for (label facei = 0; facei < nInternalFaces; facei++)
    {
        const label own = cellAddr[mesh_.faceOwner()[facei]];
        const label nbr = cellAddr[mesh_.faceNeighbour()[facei]];

        if ((own != -1) && (nbr == -1))
        {
            vector nf = Sf[facei]/magSf[facei];

            if ((nf & axis) > tol)
            {
                area_[own] += magSf[facei];
                n += Sf[facei];
            }
        }
        else if ((own == -1) && (nbr != -1))
        {
            vector nf = Sf[facei]/magSf[facei];

            if ((-nf & axis) > tol)
            {
                area_[nbr] += magSf[facei];
                n -= Sf[facei];
            }
        }
    }


    // Add boundary contributions
    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];
        const vectorField& Sfp = mesh_.Sf().boundaryField()[patchi];
        const scalarField& magSfp = mesh_.magSf().boundaryField()[patchi];

        if (pp.coupled())
        {
            forAll(pp, j)
            {
                const label facei = pp.start() + j;
                const label own = cellAddr[mesh_.faceOwner()[facei]];
                const label nbr = nbrFaceCellAddr[facei - nInternalFaces];
                const vector nf = Sfp[j]/magSfp[j];

                if ((own != -1) && (nbr == -1) && ((nf & axis) > tol))
                {
                    area_[own] += magSfp[j];
                    n += Sfp[j];
                }
            }
        }
        else
        {
            forAll(pp, j)
            {
                const label facei = pp.start() + j;
                const label own = cellAddr[mesh_.faceOwner()[facei]];
                const vector nf = Sfp[j]/magSfp[j];

                if ((own != -1) && ((nf & axis) > tol))
                {
                    area_[own] += magSfp[j];
                    n += Sfp[j];
                }
            }
        }
    }

    if (correct)
    {
        reduce(n, sumOp<vector>());
        axis = n/mag(n);
    }

    if (debug)
    {
        auto tarea = volScalarField::New
        (
            IOobject::scopedName(name_, "area"),
            IOobject::NO_REGISTER,
            mesh_,
            dimensionedScalar(dimArea, Zero)
        );
        auto& area = tarea.ref();

        UIndirectList<scalar>(area.primitiveFieldRef(), cells_) = area_;

        Info<< type() << ": " << name_
            << " writing field " << area.name() << endl;

        area.write();
    }
}


void Foam::fv::rotorDiskSource::createCoordinateSystem()
{
    // Construct the local rotor coordinate system
    vector origin(Zero);
    vector axis(Zero);
    vector refDir(Zero);

    geometryModeType gm =
        geometryModeTypeNames_.get("geometryMode", coeffs_);

    switch (gm)
    {
        case gmAuto:
        {
            // Determine rotation origin (cell volume weighted)
            scalar sumV = 0.0;
            const scalarField& V = mesh_.V();
            const vectorField& C = mesh_.C();
            forAll(cells_, i)
            {
                const label celli = cells_[i];
                sumV += V[celli];
                origin += V[celli]*C[celli];
            }
            reduce(origin, sumOp<vector>());
            reduce(sumV, sumOp<scalar>());
            origin /= sumV;

            // Determine first radial vector
            vector dx1(Zero);
            scalar magR = -GREAT;
            forAll(cells_, i)
            {
                const label celli = cells_[i];
                vector test = C[celli] - origin;
                if (mag(test) > magR)
                {
                    dx1 = test;
                    magR = mag(test);
                }
            }
            reduce(dx1, maxMagSqrOp<vector>());
            magR = mag(dx1);

            // Determine second radial vector and cross to determine axis
            forAll(cells_, i)
            {
                const label celli = cells_[i];
                vector dx2 = C[celli] - origin;
                if (mag(dx2) > 0.5*magR)
                {
                    axis = dx1 ^ dx2;
                    if (mag(axis) > SMALL)
                    {
                        break;
                    }
                }
            }
            reduce(axis, maxMagSqrOp<vector>());
            axis.normalise();

            // Correct the axis direction using a point above the rotor
            {
                vector pointAbove(coeffs_.get<vector>("pointAbove"));
                vector dir = pointAbove - origin;
                dir.normalise();
                if ((dir & axis) < 0)
                {
                    axis *= -1.0;
                }
            }

            coeffs_.readEntry("refDirection", refDir);

            // Set the face areas and apply correction to calculated axis
            // e.g. if cellZone is more than a single layer in thickness
            setFaceArea(axis, true);

            break;
        }
        case gmSpecified:
        {
            coeffs_.readEntry("origin", origin);
            coeffs_.readEntry("axis", axis);
            coeffs_.readEntry("refDirection", refDir);

            setFaceArea(axis, false);

            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unknown geometryMode " << geometryModeTypeNames_[gm]
                << ". Available geometry modes include "
                << geometryModeTypeNames_
                << exit(FatalError);
        }
    }

    coordSys_ = coordSystem::cylindrical(origin, axis, refDir);

    const scalar sumArea = gSum(area_);
    const scalar diameter = Foam::sqrt(4.0*sumArea/mathematical::pi);
    Info<< "    Rotor gometry:" << nl
        << "    - disk diameter = " << diameter << nl
        << "    - disk area     = " << sumArea << nl
        << "    - origin        = " << coordSys_.origin() << nl
        << "    - r-axis        = " << coordSys_.e1() << nl
        << "    - psi-axis      = " << coordSys_.e2() << nl
        << "    - z-axis        = " << coordSys_.e3() << endl;
}


void Foam::fv::rotorDiskSource::constructGeometry()
{
    const pointUIndList cc(mesh_.C(), cells_);

    // Optional: for later transform(), invTransform()
    /// Rcyl_.reset(coordSys_.R(cc).ptr());

    forAll(cells_, i)
    {
        if (area_[i] > ROOTVSMALL)
        {
            // Position in (planar) rotor coordinate system
            x_[i] = coordSys_.localPosition(cc[i]);

            // Cache max radius
            rMax_ = max(rMax_, x_[i].x());

            // Swept angle relative to rDir axis [radians] in range 0 -> 2*pi
            scalar psi = x_[i].y();

            // Blade flap angle [radians]
            scalar beta =
                flap_.beta0 - flap_.beta1c*cos(psi) - flap_.beta2s*sin(psi);

            // Determine rotation tensor to convert from planar system into the
            // rotor cone system
            scalar c = cos(beta);
            scalar s = sin(beta);
            Rcone_[i] = tensor(c, 0, -s, 0, 1, 0, s, 0, c);
        }
    }
}


Foam::tmp<Foam::vectorField> Foam::fv::rotorDiskSource::inflowVelocity
(
    const volVectorField& U
) const
{
    switch (inletFlow_)
    {
        case ifFixed:
        case ifSurfaceNormal:
        {
            return tmp<vectorField>::New(mesh_.nCells(), inletVelocity_);

            break;
        }
        case ifLocal:
        {
            return U.primitiveField();

            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unknown inlet flow specification" << abort(FatalError);
        }
    }

    return tmp<vectorField>::New(mesh_.nCells(), Zero);
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::fv::rotorDiskSource::rotorDiskSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh

)
:
    fv::cellSetOption(name, modelType, dict, mesh),
    rhoRef_(1.0),
    omega_(0.0),
    nBlades_(0),
    inletFlow_(ifLocal),
    inletVelocity_(Zero),
    tipEffect_(1.0),
    flap_(),
    x_(cells_.size(), Zero),
    Rcone_(cells_.size(), I),
    area_(cells_.size(), Zero),
    coordSys_(),
    rMax_(0.0),
    trim_(trimModel::New(*this, coeffs_)),
    blade_(coeffs_.subDict("blade")),
    profiles_(coeffs_.subDict("profiles"))
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::rotorDiskSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    auto tforce = volVectorField::New
    (
        IOobject::scopedName(name_, "rotorForce"),
        IOobject::NO_REGISTER,
        mesh_,
        dimensionedVector(eqn.dimensions()/dimVolume, Zero)
    );
    auto& force = tforce.ref();

    // Read the reference density for incompressible flow
    coeffs_.readEntry("rhoRef", rhoRef_);

    const vectorField Uin(inflowVelocity(eqn.psi()));
    trim_->correct(Uin, force);
    calculate(geometricOneField(), Uin, trim_->thetag(), force);

    // Add source to rhs of eqn
    eqn -= force;

    if (mesh_.time().writeTime())
    {
        force.write();
    }
}


void Foam::fv::rotorDiskSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    auto tforce = volVectorField::New
    (
        IOobject::scopedName(name_, "rotorForce"),
        IOobject::NO_REGISTER,
        mesh_,
        dimensionedVector(eqn.dimensions()/dimVolume, Zero)
    );
    auto& force = tforce.ref();

    const vectorField Uin(inflowVelocity(eqn.psi()));
    trim_->correct(rho, Uin, force);
    calculate(rho, Uin, trim_->thetag(), force);

    // Add source to rhs of eqn
    eqn -= force;

    if (mesh_.time().writeTime())
    {
        force.write();
    }
}


bool Foam::fv::rotorDiskSource::read(const dictionary& dict)
{
    if (fv::cellSetOption::read(dict))
    {
        coeffs_.readEntry("fields", fieldNames_);
        fv::option::resetApplied();

        // Read coordinate system/geometry invariant properties
        omega_ = rpmToRads(coeffs_.get<scalar>("rpm"));

        coeffs_.readEntry("nBlades", nBlades_);

        inletFlowTypeNames_.readEntry("inletFlowType", coeffs_, inletFlow_);

        coeffs_.readEntry("tipEffect", tipEffect_);

        const dictionary& flapCoeffs(coeffs_.subDict("flapCoeffs"));
        flap_.beta0 = degToRad(flapCoeffs.get<scalar>("beta0"));
        flap_.beta1c = degToRad(flapCoeffs.get<scalar>("beta1c"));
        flap_.beta2s = degToRad(flapCoeffs.get<scalar>("beta2s"));


        // Create coordinate system
        createCoordinateSystem();

        // Read coordinate system dependent properties
        checkData();

        constructGeometry();

        trim_->read(coeffs_);

        if (debug)
        {
            writeField("thetag", trim_->thetag()(), true);
            writeField("faceArea", area_, true);
        }

        return true;
    }

    return false;
}


// ************************************************************************* //
