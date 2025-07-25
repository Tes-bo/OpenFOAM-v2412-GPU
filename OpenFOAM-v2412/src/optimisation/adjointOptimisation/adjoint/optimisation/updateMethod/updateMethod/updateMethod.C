/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

#include "updateMethod.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(updateMethod, 0);
    defineRunTimeSelectionTable(updateMethod, dictionary);
}


// * * * * * * * * * *  Protected  Member Functions  * * * * * * * * * * * * //

const Foam::scalarField Foam::updateMethod::leftMult
(
    const scalarField& s,
    const SquareMatrix<scalar>& m
)
{
    if (s.size() != m.n())
    {
        FatalErrorInFunction
            << "scalar derivative and HessianInv matrix do not have the "
            << "same dimension"
            << abort(FatalError);
    }

    scalarField res(s.size(), Zero);
    forAll(s, i)
    {
        forAll(s, j)
        {
            res[i] += s[j]*m[j][i];
        }
    }

    return (res);
}


const Foam::scalarField Foam::updateMethod::rightMult
(
    const SquareMatrix<scalar>& m,
    const scalarField& s
)
{
    if (s.size() != m.n())
    {
        FatalErrorInFunction
            << "scalar derivative and HessianInv matrix do not have the "
            << "same dimension"
            << abort(FatalError);
    }

    scalarField res(s.size(), Zero);
    forAll(s, i)
    {
        forAll(s, j)
        {
            res[i] += m[i][j]*s[j];
        }
    }

    return (res);
}


Foam::SquareMatrix<Foam::scalar> Foam::updateMethod::outerProd
(
    const scalarField& a,
    const scalarField& b
)
{
    if (a.size() != b.size())
    {
        FatalErrorInFunction
            << "operands of outerProduct do not have the same dimension"
            << abort(FatalError);
    }

    SquareMatrix<scalar> res(a.size(), Zero);
    forAll(a, i)
    {
        forAll(a, j)
        {
            res[i][j]  = a[i]*b[j];
        }
    }

    return (res);
}


Foam::SquareMatrix<Foam::scalar>
Foam::updateMethod::inv(SquareMatrix<scalar> A)
{
    label n(A.n());
    SquareMatrix<scalar> invA(n, Zero);

    // LU decomposition of A
    labelList pivotIndices(n, Zero);
    LUDecompose(A, pivotIndices);
    DebugInfo
        << "LU decomposed A " << A << endl;

    // Compute inverse of A by successive back-substitutions.
    for (label j = 0; j < n; j++)
    {
        scalarField rhs(n, Zero);
        rhs[j] = scalar(1);
        LUBacksubstitute(A, pivotIndices, rhs);
        // After LUBacksubstitute, rhs contains the j-th column of the inverse
        for (label i = 0; i < n; i++)
        {
            invA[i][j] = rhs[i];
        }
    }


    /*
    // Alternative using SVD. Slower and less accurate
    tempscalarRectangularMatrix Atemp(n, n, 0);
    for (label i = 0; i < n; i++)
    {
        for (label j = 0; j < n; j++)
        {
            Atemp[i][j] = A[i][j];
        }
    }
    scalarRectangularMatrix invTemp = SVDinv(Atemp);
    scalarSquareMatrix invA(n, n, 0);
    for (label i = 0; i < n; i++)
    {
        for (label j = 0; j < n; j++)
        {
            invA[i][j] = invTemp[i][j];
        }
    }
    */

    return invA;
}


Foam::scalar Foam::updateMethod::globalSum(const scalarField& field)
{
    if (globalSum_)
    {
        return gSum(field);
    }
    return sum(field);
}


Foam::scalar Foam::updateMethod::globalSum(tmp<scalarField>& tfield)
{
    scalar value = globalSum(tfield());
    tfield.clear();
    return value;
}


Foam::label Foam::updateMethod::globalSum(const label size)
{
    label res(0);
    if (globalSum_)
    {
        res = returnReduce(size, sumOp<label>());
    }
    else
    {
        res = size;
    }
    return res;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::updateMethod::updateMethod
(
    const fvMesh& mesh,
    const dictionary& dict,
    autoPtr<designVariables>& designVars,
    const label nConstraints,
    const word& type
)
:
    localIOdictionary
    (
        IOobject
        (
            "updateMethodDict",
            mesh.time().timeName(),
            "uniform",
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        word::null // avoid type checking
    ),
    mesh_(mesh),
    dict_(dict),
    designVars_(designVars),
    nConstraints_(nConstraints),
    activeDesignVars_(designVars().activeDesignVariables()),
    objectiveDerivatives_(designVars().getVars().size(), Zero),
    constraintDerivatives_(0),
    objectiveValue_(0),
    objectiveValueOld_(nullptr),
    cValues_(0),
    correction_(readOrZeroField("correction", designVars().getVars().size())),
    cumulativeCorrection_(0),
    eta_(1),
    counter_(getOrDefault<label>("counter", Zero)),
    initialEtaSet_(false),
    correctionFolder_(mesh_.time().globalPath()/"optimisation"/"correction"),
    globalSum_(designVars_->globalSum())
{
    // Set initial eta, if present. It might be set either in the
    // optimisationDict or in the specific dictionary dedicated to the
    // updateMethod
    if (dict.readIfPresent("eta", eta_))
    {
        initialEtaSet_ = true;
    }
    else if (this->readIfPresent("eta", eta_))
    {
        initialEtaSet_ = true;
    }
}


Foam::tmp<Foam::scalarField> Foam::updateMethod::readOrZeroField
(
    const word& name,
    const label size
)
{
    return tmp<scalarField>::New(name, *this, size, IOobjectOption::LAZY_READ);
}


// * * * * * * * * * * * * * * * * Selectors  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::updateMethod> Foam::updateMethod::New
(
    const fvMesh& mesh,
    const dictionary& dict,
    autoPtr<designVariables>& designVars,
    const label nConstraints
)
{
    const word modelType(dict.get<word>("method"));

    Info<< "updateMethod type : " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "updateMethod",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<updateMethod>
        (ctorPtr(mesh, dict, designVars, nConstraints, modelType));
}


// * * * * * * * * * * * * * * *  Member Functions   * * * * * * * * * * * * //

Foam::dictionary Foam::updateMethod::coeffsDict(const word& type) const
{
    return dict_.optionalSubDict(type);
}


void Foam::updateMethod::setObjectiveDeriv(const scalarField& derivs)
{
    objectiveDerivatives_ = derivs;
}


void Foam::updateMethod::setConstraintDeriv
(
    const PtrList<scalarField>& derivs
)
{
    constraintDerivatives_ = derivs;
}


void Foam::updateMethod::setObjectiveValue(const scalar value)
{
    objectiveValue_ = value;
}


void Foam::updateMethod::setObjectiveValueOld(const scalar value)
{
    if (!objectiveValueOld_)
    {
        objectiveValueOld_.reset(new scalar(Zero));
    }
    objectiveValueOld_.ref() = value;
}


void Foam::updateMethod::setConstraintValues(const scalarField& values)
{
    cValues_ = values;
}


Foam::scalar Foam::updateMethod::getObjectiveValue() const
{
    return objectiveValue_;
}


const Foam::autoPtr<Foam::scalar>&
Foam::updateMethod::getObjectiveValueOld() const
{
    return objectiveValueOld_;
}


const Foam::scalarField& Foam::updateMethod::getConstraintValues() const
{
    return cValues_;
}


Foam::label Foam::updateMethod::getCycle() const
{
    return counter_;
}


void Foam::updateMethod::setStep(const scalar eta)
{
    eta_ = eta;
}


void Foam::updateMethod::modifyStep(const scalar multiplier)
{
    eta_ *= multiplier;
}


void Foam::updateMethod::setGlobalSum(const bool useGlobalSum)
{
    globalSum_ = useGlobalSum;
}


void Foam::updateMethod::setConstaintsNumber(const label nConstraints)
{
    nConstraints_ = nConstraints;
}


Foam::label Foam::updateMethod::nConstraints() const
{
    return nConstraints_;
}


Foam::scalarField& Foam::updateMethod::returnCorrection()
{
    return correction_;
}


void Foam::updateMethod::writeCorrection()
{
    if (Pstream::master())
    {
        // Allocate cumulativeCorrection if necessary
        if (cumulativeCorrection_.empty())
        {
            cumulativeCorrection_.setSize(correction_.size(), Zero);
        }
        // Accumulate correction
        cumulativeCorrection_ += correction_;

        fileName correctionFile
        (
            correctionFolder_/"correction" + mesh_.time().timeName()
        );
        fileName cumulativeCorrectionFile
        (
            correctionFolder_/"cumulativeCorrection" + mesh_.time().timeName()
        );

        OFstream corFile(correctionFile);
        OFstream cumulCorFile(cumulativeCorrectionFile);
        forAll(correction_, cI)
        {
            corFile
                << cI << " " << correction_[cI] << endl;
            cumulCorFile
                << cI << " " << cumulativeCorrection_[cI] << endl;
        }
    }
}


Foam::scalar Foam::updateMethod::computeMeritFunction()
{
    return objectiveValue_;
}


Foam::scalar Foam::updateMethod::meritFunctionDirectionalDerivative()
{
    return globalSum(objectiveDerivatives_*correction_);
}


bool& Foam::updateMethod::initialEtaSet()
{
    return initialEtaSet_;
}


void Foam::updateMethod::updateOldCorrection
(
    const scalarField& oldCorrection
)
{
    correction_ = oldCorrection;
}


bool Foam::updateMethod::writeData(Ostream& os) const
{
    // Insert eta if set
    if (initialEtaSet_)
    {
        os.writeEntry("eta", eta_);
    }

    os.writeEntry("counter", counter_);
    correction_.writeEntry("correction", os);

    return true;
}


bool Foam::updateMethod::writeAuxiliaryData()
{
    // Does nothing in base
    return true;
}


// ************************************************************************* //
