/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
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

#include "LBFGS.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(LBFGS, 1);
    addToRunTimeSelectionTable
    (
        updateMethod,
        LBFGS,
        dictionary
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::LBFGS::allocateVectors()
{
    label nVars(activeDesignVars_.size());
    for (label i = 0; i < nPrevSteps_; ++i)
    {
        if (!y_.get(i))
        {
            y_.set(i, new scalarField(nVars, Zero));
        }
        if (!s_.get(i))
        {
            s_.set(i, new scalarField(nVars, Zero));
        }
        if (found("y" + Foam::name(i)))
        {
            y_[i] = scalarField("y" + Foam::name(i), *this, nVars);
        }
        if (found("s" + Foam::name(i)))
        {
            s_[i] = scalarField("s" + Foam::name(i), *this, nVars);
        }
    }
}


void Foam::LBFGS::pivotFields(PtrList<scalarField>& list, const scalarField& f)
{
    if (counter_ > nPrevSteps_)
    {
        // Reorder list by moving pointers down the line
        labelList newOrder(nPrevSteps_, -1);
        newOrder[0] = nPrevSteps_ - 1;
        for (label i = 1; i < nPrevSteps_; ++i)
        {
            newOrder[i] = i - 1;
        }
        list.reorder(newOrder);

        // Fill in last element with the provided field
        list[nPrevSteps_ - 1] = f;
    }
    else
    {
        list[counter_ - 1] = f;
    }
}


void Foam::LBFGS::updateVectors
(
    const scalarField& derivatives,
    const scalarField& derivativesOld
)
{
    // Sanity checks
    if
    (
        (derivatives.size() != derivativesOld.size())
     || (derivatives.size() != designVars_().getVars().size())
    )
    {
        FatalErrorInFunction
            << "Sizes of input derivatives and design variables do not match"
            << exit(FatalError);
    }

    // Update list of y. Can only be done here since derivatives
    // were not known at the end of the previous cycle
    scalarField yRecent(derivatives - derivativesOld, activeDesignVars_);
    // Update list of s.
    // correction_ holds the previous correction
    scalarField sActive(correctionOld_, activeDesignVars_);
    applyDamping(yRecent, sActive);

    pivotFields(y_, yRecent);
    pivotFields(s_, sActive);
}


void Foam::LBFGS::applyDamping(scalarField& y, scalarField& s)
{
    const scalar sy(globalSum(s*y));
    if (useSDamping_)
    {
        const scalarField Hy(invHessianVectorProduct(y, counter_ - 1));
        const scalar yHy(globalSum(y*Hy));
        scalar theta(1);
        if (sy < 0.2*yHy)
        {
            WarningInFunction
                << "y*s is below threshold. Using damped form" << nl
                << "sy, yHy " << sy << " " << yHy << endl;

            theta = 0.8*yHy/(yHy - sy);
        }
        s = theta*s + (1 - theta)*Hy;
    }
    else if (useYDamping_)
    {
        const scalarField Bs(HessianVectorProduct(s, counter_ - 1));
        const scalar sBs(globalSum(s*Bs));
        scalar theta(1);
        if (sy < 0.2*sBs)
        {
            WarningInFunction
                << "y*s is below threshold. Using damped form" << nl
                << "sy, sBs " << sy << " " << sBs << endl;

            theta = 0.8*sBs/(sBs - sy);
        }
        y = theta*y + (1 - theta)*Bs;
    }
    DebugInfo
        << "Curvature index (sy) is " << sy << endl;
}


Foam::tmp<Foam::scalarField>
Foam::LBFGS::invHessianVectorProduct(const scalarField& vector)
{
    return invHessianVectorProduct(vector, counter_);
}


Foam::tmp<Foam::scalarField>
Foam::LBFGS::invHessianVectorProduct
(
    const scalarField& vector,
    const label counter,
    const refPtr<scalarField> diag
)
{
    // Sanity checks
    tmp<scalarField> tq(tmp<scalarField>::New(activeDesignVars_.size(), Zero));
    scalarField& q = tq.ref();
    label nv = designVars_().getVars().size();
    label nav = activeDesignVars_.size();
    if (vector.size() == nv)
    {
        q.map(vector, activeDesignVars_);
    }
    else if (vector.size() == nav)
    {
        q = vector;
    }
    else
    {
        FatalErrorInFunction
            << "Size of input vector "
            << "(" << vector.size() << ") "
            << "is equal to neither the number of design variabes "
            << "(" << nv << ")"
            << " nor that of the active design variables "
            << "(" << nav << ")"
            << exit(FatalError);
    }

    if (counter)
    {
        // L-BFGS two loop recursion
        //~~~~~~~~~~~~~~~~~~~~~~~~~~
        label nSteps(min(counter, nPrevSteps_));
        label nLast(nSteps - 1);
        scalarField a(nSteps, 0.);
        scalarField r(nSteps, 0.);
        for (label i = nLast; i > -1; --i)
        {
            r[i] = 1./globalSum(y_[i]*s_[i]);
            a[i] = r[i]*globalSum(s_[i]*q);
            q -= a[i]*y_[i];
        }

        scalar gamma =
            globalSum(y_[nLast]*y_[nLast])/globalSum(y_[nLast]*s_[nLast]);
        if (diag)
        {
            q /= (gamma + diag());
        }
        else
        {
            q /= gamma;
        }

        scalarField b(activeDesignVars_.size(), Zero);
        for (label i = 0; i < nSteps; ++i)
        {
            b = r[i]*globalSum(y_[i]*q);
            q += s_[i]*(a[i] - b);
        }
    }
    else if (diag)
    {
        q /= (1 + diag());
    }

    return tq;
}


Foam::tmp<Foam::scalarField>
Foam::LBFGS::HessianVectorProduct(const scalarField& vector)
{
    return HessianVectorProduct(vector, counter_);
}


Foam::tmp<Foam::scalarField>
Foam::LBFGS::HessianVectorProduct
(
    const scalarField& vector,
    const label counter
)
{
    addProfiling(LBFGS, "LBFGS::HessianVectorProduct");
    // Sanity checks
    tmp<scalarField> tq(tmp<scalarField>::New(activeDesignVars_.size(), Zero));
    scalarField& q = tq.ref();

    scalarField source;
    if (vector.size() == designVars_().getVars().size())
    {
        source = scalarField(vector, activeDesignVars_);
    }
    else if (vector.size() == activeDesignVars_.size())
    {
        source = vector;
    }
    else
    {
        FatalErrorInFunction
            << "Size of input vector is equal to neither the number of "
            << " design variabes nor that of the active design variables"
            << exit(FatalError);
    }

    if (counter != 0)
    {
        const label nSteps(min(counter, nPrevSteps_));
        const label nLast(nSteps - 1);
        const scalar delta =
            globalSum(y_[nLast]*y_[nLast])/globalSum(y_[nLast]*s_[nLast]);

        // Product of the last matrix on the right with the input vector
        scalarField SKsource(2*nSteps, Zero);
        for (label i = 0; i < nSteps; ++i)
        {
            SKsource[i] = delta*globalSum(s_[i]*source);
            SKsource[i + nSteps] = globalSum(y_[i]*source);
        }

        // Form the middle matrix to be inverted
        SquareMatrix<scalar> M(2*nSteps, 2*nSteps, Zero);
        for (label i = 0; i < nSteps; ++i)
        {
            // Lower diagonal part
            M[nSteps + i][nSteps + i] = - globalSum(s_[i]*y_[i]);
            // Upper left part
            for (label j = 0; j < nSteps; ++j)
            {
                M[i][j] = delta*globalSum(s_[i]*s_[j]);
            }
        }

        // Upper right and lower left parts
        for (label j = 0; j < nSteps; ++j)
        {
            for (label i = j + 1; i < nSteps; ++i)
            {
                scalar value = globalSum(s_[i]*y_[j]);
                M[i][j + nSteps] = value;
                M[j + nSteps][i] = value;
            }
        }
        SquareMatrix<scalar> invM(inv(M));

        // Product of the inverted middle matrix with the right vector
        scalarField invMSource(rightMult(invM, SKsource));

        // Left vector multiplication with the rest of contributions
        // vag: parallel comms
        forAll(q, i)
        {
            for (label j = 0; j < nSteps; ++j)
            {
                q[i] -=
                    delta*s_[j][i]*invMSource[j]
                  + y_[j][i]*invMSource[j + nSteps];
            }
        }

        q += delta*source;
    }
    else
    {
        q = source;
    }

    return tq;
}


Foam::tmp<Foam::scalarField> Foam::LBFGS::HessianDiag()
{
    // Sanity checks
    const label n(activeDesignVars_.size());
    tmp<scalarField> tdiag(tmp<scalarField>::New(n, 1));
    scalarField& diag = tdiag.ref();

    if (counter_ != 0)
    {
        const label nSteps(min(counter_, nPrevSteps_));
        const label nLast(nSteps - 1);
        const scalar delta =
            globalSum(y_[nLast]*y_[nLast])/globalSum(y_[nLast]*s_[nLast]);
        diag *= delta;

        // Form the middle matrix to be inverted
        SquareMatrix<scalar> M(2*nSteps, 2*nSteps, Zero);
        for (label i = 0; i < nSteps; ++i)
        {
            // Lower diagonal part
            M[nSteps + i][nSteps + i] = - globalSum(s_[i]*y_[i]);
            // Upper left part
            for (label j = 0; j < nSteps; ++j)
            {
                M[i][j] = delta*globalSum(s_[i]*s_[j]);
            }
        }

        // Upper right and lower left parts
        for (label j = 0; j < nSteps; ++j)
        {
            for (label i = j + 1; i < nSteps; ++i)
            {
                scalar value = globalSum(s_[i]*y_[j]);
                M[i][j + nSteps] = value;
                M[j + nSteps][i] = value;
            }
        }

        // Invert the matrix
        SquareMatrix<scalar> invM(inv(M));

        // Product of the inverse of the middle matrix with the right vector
        List<scalarField> MR(2*nSteps, scalarField(n, Zero));
        for (label k = 0; k < n; ++k)
        {
            for (label i = 0; i < 2*nSteps; ++i)
            {
                for (label j = 0; j < nSteps; ++j)
                {
                    MR[i][k] +=
                        invM[i][j]*delta*s_[j][k]
                      + invM[i][j + nSteps]*y_[j][k];
                }
            }
        }

        // Part of the Hessian diagonal computed by the multiplication
        // of the above matrix with the left matrix of the recursive Hessian
        // reconstruction
        for (label k = 0; k < n; ++k)
        {
            for (label j = 0; j < nSteps; ++j)
            {
                diag[k] -=
                    delta*s_[j][k]*MR[j][k] + y_[j][k]*MR[j + nSteps][k];
            }
        }
    }

    return tdiag;
}


Foam::tmp<Foam::scalarField>
Foam::LBFGS::SR1HessianVectorProduct(const scalarField& vector)
{
    return SR1HessianVectorProduct(vector, counter_);
}


Foam::tmp<Foam::scalarField>
Foam::LBFGS::SR1HessianVectorProduct
(
    const scalarField& vector,
    const label counter
)
{
    // Sanity checks
    tmp<scalarField> tq(tmp<scalarField>::New(activeDesignVars_.size(), Zero));
    scalarField& q = tq.ref();

    scalarField source;
    if (vector.size() == designVars_().getVars().size())
    {
        source = scalarField(vector, activeDesignVars_);
    }
    else if (vector.size() == activeDesignVars_.size())
    {
        source = vector;
    }
    else
    {
        FatalErrorInFunction
            << "Size of input vector is equal to neither the number of "
            << " design variabes nor that of the active design variables"
            << exit(FatalError);
    }

    if (counter != 0)
    {
        const label nSteps(min(counter, nPrevSteps_));
        const label nLast(nSteps - 1);
        const scalar delta =
            globalSum(y_[nLast]*y_[nLast])/globalSum(y_[nLast]*s_[nLast]);

        // Product of the last matrix on the right with the input vector
        scalarField YBSsource(nSteps, Zero);
        for (label i = 0; i < nSteps; ++i)
        {
            YBSsource[i] = globalSum((y_[i] - delta*s_[i])*source);
        }

        // Form the middle matrix to be inverted
        SquareMatrix<scalar> M(nSteps, nSteps, Zero);
        for (label i = 0; i < nSteps; ++i)
        {
            // D part
            M[i][i] += globalSum(s_[i]*y_[i]);
            // (S^T)BS part
            for (label j = 0; j < nSteps; ++j)
            {
                M[i][j] -= delta*globalSum(s_[i]*s_[j]);
            }
        }

        // Upper right and lower left parts
        for (label j = 0; j < nSteps; ++j)
        {
            for (label i = j + 1; i < nSteps; ++i)
            {
                scalar value = globalSum(s_[i]*y_[j]);
                M[i][j] += value;
                M[j][i] += value;
            }
        }
        SquareMatrix<scalar> invM(inv(M));

        // Product of the inverted middle matrix with the right vector
        scalarField invMSource(rightMult(invM, YBSsource));

        // Left vector multiplication with the rest of contributions
        // vag: parallel comms
        forAll(q, i)
        {
            for (label j = 0; j < nSteps; ++j)
            {
                q[i] += (y_[j][i] - delta*s_[j][i])*invMSource[j];
            }
        }

        q += delta*source;
    }
    else
    {
        q = source;
    }

    return tq;
}


Foam::tmp<Foam::scalarField> Foam::LBFGS::SR1HessianDiag()
{
    // Sanity checks
    const label n(activeDesignVars_.size());
    tmp<scalarField> tdiag(tmp<scalarField>::New(n, 1));
    scalarField& diag = tdiag.ref();

    if (counter_ != 0)
    {
        const label nSteps(min(counter_, nPrevSteps_));
        const label nLast(nSteps - 1);
        const scalar delta =
            globalSum(y_[nLast]*y_[nLast])/globalSum(y_[nLast]*s_[nLast]);
        diag *= delta;

        // Form the middle matrix to be inverted
        SquareMatrix<scalar> M(nSteps, nSteps, Zero);
        for (label i = 0; i < nSteps; ++i)
        {
            // D part
            M[i][i] += globalSum(s_[i]*y_[i]);
            // (S^T)BS part
            for (label j = 0; j < nSteps; ++j)
            {
                M[i][j] -= delta*globalSum(s_[i]*s_[j]);
            }
        }

        // Upper right and lower left parts
        for (label j = 0; j < nSteps; ++j)
        {
            for (label i = j + 1; i < nSteps; ++i)
            {
                scalar value = globalSum(s_[i]*y_[j]);
                M[i][j] += value;
                M[j][i] += value;
            }
        }
        SquareMatrix<scalar> invM(inv(M));

        // Product of the inverse of the middle matrix with the right vector
        List<scalarField> MR(nSteps, scalarField(n, Zero));
        for (label k = 0; k < n; ++k)
        {
            for (label i = 0; i < nSteps; ++i)
            {
                for (label j = 0; j < nSteps; ++j)
                {
                    MR[i][k] += invM[i][j]*(y_[j][k] - delta*s_[j][k]);
                }
            }
        }

        // Part of the Hessian diagonal computed by the multiplication
        // of the above matrix with the left matrix of the recursive Hessian
        // reconstruction
        for (label k = 0; k < n; ++k)
        {
            for (label j = 0; j < nSteps; ++j)
            {
                diag[k] += (y_[j][k] - delta*s_[j][k])*MR[j][k];
            }
        }
    }

    return tdiag;
}


void Foam::LBFGS::updateHessian()
{
    updateVectors(objectiveDerivatives_, derivativesOld_);
}


void Foam::LBFGS::update()
{
    // In the first few iterations, use steepest descent but update the Hessian
    // matrix
    if (counter_ < nSteepestDescent_)
    {
        Info<< "Using steepest descent to update design variables" << endl;
        for (const label varI : activeDesignVars_)
        {
            correction_[varI] = -eta_*objectiveDerivatives_[varI];
        }
    }
    // else use LBFGS formula to update the design variables
    else
    {
        scalarField q(invHessianVectorProduct(objectiveDerivatives_));
        forAll(activeDesignVars_, varI)
        {
            correction_[activeDesignVars_[varI]] = -etaHessian_*q[varI];
        }
    }

    // Store fields for the next iteration
    derivativesOld_ = objectiveDerivatives_;
    correctionOld_ = correction_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LBFGS::LBFGS
(
    const fvMesh& mesh,
    const dictionary& dict,
    autoPtr<designVariables>& designVars,
    const label nConstraints,
    const word& type
)
:
    quasiNewton(mesh, dict, designVars, nConstraints, type),
    nPrevSteps_(coeffsDict(type).getOrDefault<label>("nPrevSteps", 10)),
    y_(nPrevSteps_),
    s_(nPrevSteps_),
    useSDamping_(coeffsDict(type).getOrDefault<bool>("useSDamping", false)),
    useYDamping_(coeffsDict(type).getOrDefault<bool>("useYDamping", false))
{
    // Allocate the correct sizes for y and s
    allocateVectors();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::LBFGS::writeData(Ostream& os) const
{
    // Write each component of y and s as a separate field so as to allow for
    // reading them also in binary, since PtrList does not support this
    forAll(y_, i)
    {
        y_[i].writeEntry(word("y" + Foam::name(i)), os);
        s_[i].writeEntry(word("s" + Foam::name(i)), os);
    }

    return quasiNewton::writeData(os);
}


// ************************************************************************* //
