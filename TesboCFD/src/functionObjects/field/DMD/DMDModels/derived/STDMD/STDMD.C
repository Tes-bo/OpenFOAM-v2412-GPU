/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2022 OpenCFD Ltd.
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

#include "STDMD.H"
#include "EigenMatrix.H"
#include "QRMatrix.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace DMDModels
{
    defineTypeNameAndDebug(STDMD, 0);
    addToRunTimeSelectionTable(DMDModel, STDMD, dictionary);
}
}

const Foam::Enum
<
    Foam::DMDModels::STDMD::modeSorterType
>
Foam::DMDModels::STDMD::modeSorterTypeNames
({
    { modeSorterType::FIRST_SNAPSHOT, "firstSnapshot" },
    { modeSorterType::KIEWAT , "kiewat" },
    { modeSorterType::KOU_ZHANG , "kouZhang" }
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::DMDModels::STDMD::L2norm(const RMatrix& z) const
{
    #ifdef FULLDEBUG
    // Check if the given RectangularMatrix is effectively a column vector
    if (z.n() != 1)
    {
        FatalErrorInFunction
            << "Input matrix is not a column vector."
            << exit(FatalError);
    }
    #endif
    const bool noSqrt = true;
    scalar result = z.columnNorm(0, noSqrt);
    reduce(result, sumOp<scalar>());

    // Heuristic addition to avoid very small or zero norm
    return max(SMALL, Foam::sqrt(result));
}


Foam::RectangularMatrix<Foam::scalar> Foam::DMDModels::STDMD::orthonormalise
(
    RMatrix ez
) const
{
    List<scalar> dz(Q_.n());
    const label sz0 = ez.m();
    const label sz1 = dz.size();

    for (label i = 0; i < nGramSchmidt_; ++i)
    {
        // dz = Q_ & ez;
        dz = Zero;
        for (label i = 0; i < sz0; ++i)
        {
            for (label j = 0; j < sz1; ++j)
            {
                dz[j] += Q_(i,j)*ez(i,0);
            }
        }

        reduce(dz, sumOp<List<scalar>>());

        // ez -= Q_*dz;
        for (label i = 0; i < sz0; ++i)
        {
            scalar t = 0;
            for (label j = 0; j < sz1; ++j)
            {
                t += Q_(i,j)*dz[j];
            }
            ez(i,0) -= t;
        }
    }

    return ez;
}


void Foam::DMDModels::STDMD::expand(const RMatrix& ez, const scalar ezNorm)
{
    Info<< tab << "Expanding orthonormal basis for field: " << fieldName_
        << endl;

    // Stack a column "(ez/ezNorm)" to "Q"
    Q_.resize(Q_.m(), Q_.n() + 1);
    Q_.subColumn(Q_.n() - 1) = ez/ezNorm;

    // Stack a row (Zero) and column (Zero) to "G"
    G_.resize(G_.m() + 1);
}


void Foam::DMDModels::STDMD::updateG(const RMatrix& z)
{
    List<scalar> zTilde(Q_.n(), Zero);
    const label sz0 = z.m();
    const label sz1 = zTilde.size();

    // zTilde = Q_ & z;
    for (label i = 0; i < sz0; ++i)
    {
        for (label j = 0; j < sz1; ++j)
        {
            zTilde[j] += Q_(i,j)*z(i,0);
        }
    }

    reduce(zTilde, sumOp<List<scalar>>());

    // G_ += SMatrix(zTilde^zTilde);
    for (label i = 0; i < G_.m(); ++i)
    {
        for (label j = 0; j < G_.n(); ++j)
        {
            G_(i, j) += zTilde[i]*zTilde[j];
        }
    }
}


void Foam::DMDModels::STDMD::compress()
{
    Info<< tab << "Compressing orthonormal basis for field: " << fieldName_
        << endl;

    RMatrix q(1, 1, Zero);

    if (Pstream::master())
    {
        const bool symmetric = true;
        const EigenMatrix<scalar> EM(G_, symmetric);
        const SquareMatrix<scalar>& EVecs = EM.EVecs();
        DiagonalMatrix<scalar> EVals(EM.EValsRe());

        // Sort eigenvalues in descending order, and track indices
        const auto descend = [&](scalar a, scalar b){ return a > b; };
        const List<label> permutation(EVals.sortPermutation(descend));
        EVals.applyPermutation(permutation);
        EVals.resize(EVals.size() - 1);

        // Update "G"
        G_ = SMatrix(maxRank_, Zero);
        G_.diag(EVals);

        q.resize(Q_.n(), maxRank_);
        for (label i = 0; i < maxRank_; ++i)
        {
            q.subColumn(i) = EVecs.subColumn(permutation[i]);
        }
    }
    Pstream::broadcast(G_);
    Pstream::broadcast(q);

    // Update "Q"
    Q_ = Q_*q;
}


Foam::SquareMatrix<Foam::scalar> Foam::DMDModels::STDMD::
reducedKoopmanOperator()
{
    Info<< tab << "Computing Atilde" << endl;

    Info<< tab << "Computing local Rx" << endl;

    RMatrix Rx(1, 1, Zero);
    {
        QRMatrix<RMatrix> QRM
        (
            Qupper_,
            QRMatrix<RMatrix>::modes::ECONOMY,
            QRMatrix<RMatrix>::outputs::ONLY_R
        );
        RMatrix& R = QRM.R();
        Rx.transfer(R);
    }
    Rx.round();

    // Convenience objects A1, A2, A3 to compute "Atilde"
    RMatrix A1(1, 1, Zero);

    if (Pstream::parRun())
    {
        // Parallel direct tall-skinny QR decomposition
        // (BGD:Fig. 5) & (DGHL:Fig. 2)
        // Tests revealed that the distribution of "Q" does not affect
        // the final outcome of TSQR decomposition up to sign

        // Don't clear storage on persistent buffer
        PstreamBuffers pBufs;
        pBufs.allowClearRecv(false);

        const label myProcNo = Pstream::myProcNo();
        const label procNoInSubset = myProcNo % nAgglomerationProcs_;

        // Send Rx from sender subset neighbours to receiver subset master
        if (procNoInSubset != 0)
        {
            const label procNoSubsetMaster = myProcNo - procNoInSubset;

            UOPstream toSubsetMaster(procNoSubsetMaster, pBufs);
            toSubsetMaster << Rx;

            Rx.clear();
        }

        pBufs.finishedSends();

        // Receive Rx by receiver subset masters
        if (procNoInSubset == 0)
        {
            // Accept only subset masters possessing sender subset neighbours
            if (myProcNo + 1 < Pstream::nProcs())
            {
                for
                (
                    label nbr = myProcNo + 1;
                    (
                        nbr < myProcNo + nAgglomerationProcs_
                     && nbr < Pstream::nProcs()
                    );
                    ++nbr
                )
                {
                    RMatrix recvMtrx;

                    UIPstream fromNbr(nbr, pBufs);
                    fromNbr >> recvMtrx;

                    // Append received Rx to Rx of receiver subset masters
                    if (recvMtrx.size() > 0)
                    {
                        Rx.resize(Rx.m() + recvMtrx.m(), Rx.n());
                        Rx.subMatrix
                        (
                            Rx.m() - recvMtrx.m(),
                            Rx.n() - recvMtrx.n()
                        ) = recvMtrx;
                    }
                }
            }

            {
                // Apply interim QR decomposition
                // on Rx of receiver subset masters
                QRMatrix<RMatrix> QRM
                (
                    QRMatrix<RMatrix>::modes::ECONOMY,
                    QRMatrix<RMatrix>::outputs::ONLY_R,
                    QRMatrix<RMatrix>::pivoting::FALSE,
                    Rx
                );
                RMatrix& R = QRM.R();
                Rx.transfer(R);
            }
            Rx.round();
        }

        pBufs.clear();

        // Send interim Rx from subset masters to the master
        if (procNoInSubset == 0)
        {
            if (!Pstream::master())
            {
                UOPstream toMaster(Pstream::masterNo(), pBufs);
                toMaster << Rx;
            }
        }

        pBufs.finishedSends();


        // Receive interim Rx by the master, and apply final operations
        if (Pstream::master())
        {
            for
            (
                label nbr = Pstream::masterNo() + nAgglomerationProcs_;
                nbr < Pstream::nProcs();
                nbr += nAgglomerationProcs_
            )
            {
                RMatrix recvMtrx;

                UIPstream fromSubsetMaster(nbr, pBufs);
                fromSubsetMaster >> recvMtrx;

                // Append received Rx to Rx of the master
                if (recvMtrx.size() > 0)
                {
                    Rx.resize(Rx.m() + recvMtrx.m(), Rx.n());
                    Rx.subMatrix
                    (
                        Rx.m() - recvMtrx.m(),
                        Rx.n() - recvMtrx.n()
                    ) = recvMtrx;
                }
            }

            Info<< tab << "Computing TSQR" << endl;
            {
                QRMatrix<RMatrix> QRM
                (
                    QRMatrix<RMatrix>::modes::ECONOMY,
                    QRMatrix<RMatrix>::outputs::ONLY_R,
                    QRMatrix<RMatrix>::pivoting::FALSE,
                    Rx
                );
                RMatrix& R = QRM.R();
                Rx.transfer(R);
            }
            Rx.round();

            Info<< tab << "Computing RxInv" << endl;
            RxInv_ = MatrixTools::pinv(Rx);

            Info<< tab << "Computing A1" << endl;
            A1 = RxInv_*RMatrix(MatrixTools::pinv(Rx*(G_^Rx)));
            Rx.clear();
        }
        Pstream::broadcast(RxInv_);
        Pstream::broadcast(A1);

        Info<< tab << "Computing A2" << endl;
        SMatrix A2(Qupper_ & Qlower_);
        reduce(A2, sumOp<SMatrix>());
        Qlower_.clear();

        Info<< tab << "Computing A3" << endl;
        SMatrix A3(Qupper_ & Qupper_);
        reduce(A3, sumOp<SMatrix>());

        Info<< tab << "Computing Atilde" << endl;
        // by optimized matrix chain multiplication
        // obtained by dynamic programming memoisation
        return SMatrix(RxInv_ & ((A2*(G_*A3))*A1));
    }
    else
    {
        Info<< tab << "Computing RxInv" << endl;
        RxInv_ = MatrixTools::pinv(Rx);

        Info<< tab << "Computing A1" << endl;
        A1 = RxInv_*RMatrix(MatrixTools::pinv(Rx*(G_^Rx)));
        Rx.clear();

        Info<< tab << "Computing A2" << endl;
        SMatrix A2(Qupper_ & Qlower_);
        Qlower_.clear();

        Info<< tab << "Computing A3" << endl;
        SMatrix A3(Qupper_ & Qupper_);

        Info<< tab << "Computing Atilde" << endl;
        return SMatrix(RxInv_ & ((A2*(G_*A3))*A1));
    }
}


bool Foam::DMDModels::STDMD::eigendecomposition(SMatrix& Atilde)
{
    bool fail = false;

    // Compute eigenvalues, and clip eigenvalues with values < "minEval"
    if (Pstream::master())
    {
        Info<< tab << "Computing eigendecomposition" << endl;

        // Replace "Atilde" by eigenvalues (in-place eigendecomposition)
        const EigenMatrix<scalar> EM(Atilde);
        const DiagonalMatrix<scalar>& evalsRe = EM.EValsRe();
        const DiagonalMatrix<scalar>& evalsIm = EM.EValsIm();

        evals_.resize(evalsRe.size());
        evecs_ = RCMatrix(EM.complexEVecs());

        // Convert scalar eigenvalue pairs to complex eigenvalues
        label i = 0;
        for (auto& eval : evals_)
        {
            eval = complex(evalsRe[i], evalsIm[i]);
            ++i;
        }

        Info<< tab << "Filtering eigenvalues" << endl;

        List<complex> cp(evals_.size());
        auto it =
            std::copy_if
            (
                evals_.cbegin(),
                evals_.cend(),
                cp.begin(),
                [&](const complex& x){ return mag(x) > minEval_; }
            );
        cp.resize(std::distance(cp.begin(), it));

        if (cp.size() == 0)
        {
            WarningInFunction
                << "No eigenvalue with mag(eigenvalue) larger than "
                << "minEval = " << minEval_ << " was found." << nl
                << "    Input field may contain only zero-value elements," << nl
                << "    e.g. no-slip velocity condition on a given patch." << nl
                << "    Otherwise, please decrease the value of minEval." << nl
                << "    Skipping further dynamics/mode computations." << nl
                << endl;

            fail = true;
        }
        else
        {
            evals_ = cp;
        }
    }
    Pstream::broadcast(fail);

    if (fail)
    {
        return false;
    }

    Pstream::broadcast(evals_);
    Pstream::broadcast(evecs_);

    return true;
}


void Foam::DMDModels::STDMD::frequencies()
{
    if (Pstream::master())
    {
        Info<< tab << "Computing frequencies" << endl;

        freqs_.resize(evals_.size());

        // Frequency equation (K:Eq. 81)
        auto frequencyEquation =
            [&](const complex& eval)
            {
                return Foam::log(max(eval, complex(SMALL))).imag()/(twoPi*dt_);
            };

        // Compute frequencies
        std::transform
        (
            evals_.cbegin(),
            evals_.cend(),
            freqs_.begin(),
            frequencyEquation
        );

        Info<< tab << "Computing frequency indices" << endl;

        // Initialise iterator by the first search
        auto margin = [&](const scalar& x){ return (fMin_ <= x && x < fMax_); };

        auto it = std::find_if(freqs_.cbegin(), freqs_.cend(), margin);

        while (it != freqs_.end())
        {
            freqsi_.append(std::distance(freqs_.cbegin(), it));
            it = std::find_if(std::next(it), freqs_.cend(), margin);
        }
    }
    Pstream::broadcast(freqs_);
    Pstream::broadcast(freqsi_);
}


void Foam::DMDModels::STDMD::amplitudes()
{
    IOField<scalar> snapshot0
    (
        IOobject
        (
            "snapshot0_" + name_ + "_" + fieldName_,
            timeName0_,
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    // RxInv^T*(Qupper^T*snapshot0)
    // tmp = (Qupper^T*snapshot0)
    List<scalar> tmp(Qupper_.n(), Zero);
    const label sz0 = snapshot0.size();
    const label sz1 = tmp.size();

    for (label i = 0; i < sz0; ++i)
    {
        for (label j = 0; j < sz1; ++j)
        {
            tmp[j] += Qupper_(i,j)*snapshot0[i];
        }
    }
    snapshot0.clear();

    // R = RxInv^T*tmp
    List<scalar> R(Qupper_.n(), Zero);
    for (label i = 0; i < sz1; ++i)
    {
        for (label j = 0; j < R.size(); ++j)
        {
            R[j] += RxInv_(i,j)*tmp[i];
        }
    }
    tmp.clear();

    reduce(R, sumOp<List<scalar>>());

    if (Pstream::master())
    {
        Info<< tab << "Computing amplitudes" << endl;

        amps_.resize(R.size());
        const RCMatrix pEvecs(MatrixTools::pinv(evecs_));

        // amps_ = pEvecs*R;
        for (label i = 0; i < amps_.size(); ++i)
        {
            for (label j = 0; j < R.size(); ++j)
            {
                amps_[i] += pEvecs(i,j)*R[j];
            }
        }
    }
    Pstream::broadcast(amps_);
}


void Foam::DMDModels::STDMD::magnitudes()
{
    if (Pstream::master())
    {
        Info<< tab << "Computing magnitudes" << endl;

        mags_.resize(amps_.size());

        Info<< tab << "Sorting modes with ";

        switch (modeSorter_)
        {
            case modeSorterType::FIRST_SNAPSHOT:
            {
                Info<< "method of first snapshot" << endl;

                std::transform
                (
                    amps_.cbegin(),
                    amps_.cend(),
                    mags_.begin(),
                    [&](const complex& val){ return mag(val); }
                );
                break;
            }

            case modeSorterType::KIEWAT:
            {
                Info<< "modified weighted amplitude scaling method" << endl;

                //  Eigendecomposition returns evecs with
                //  the unity norm, hence "modeNorm = 1"
                const scalar modeNorm = 1;
                const scalar pr = 1;
                List<scalar> w(step_);
                std::iota(w.begin(), w.end(), 1);
                w = sin(twoPi/step_*(w - 1 - 0.25*step_))*pr + pr;

                forAll(mags_, i)
                {
                    mags_[i] = sorter(w, amps_[i], evals_[i], modeNorm);
                }

                break;
            }

            case modeSorterType::KOU_ZHANG:
            {
                Info<< "weighted amplitude scaling method" << endl;

                const scalar modeNorm = 1;
                const List<scalar> w(step_, 1);

                forAll(mags_, i)
                {
                    mags_[i] = sorter(w, amps_[i], evals_[i], modeNorm);
                }

                break;
            }

            default:
                break;
        }

        Info<< tab << "Computing magnitude indices" << endl;

        magsi_ = freqsi_;

        auto descend =
            [&](const label i1, const label i2)
            {
                return !(mags_[i1] < mags_[i2]);
            };

        std::sort(magsi_.begin(), magsi_.end(), descend);
    }
    Pstream::broadcast(mags_);
    Pstream::broadcast(magsi_);
}


Foam::scalar Foam::DMDModels::STDMD::sorter
(
    const List<scalar>& weight,
    const complex& amplitude,
    const complex& eigenvalue,
    const scalar modeNorm
) const
{
    // Omit eigenvalues with very large or very small mags
    if (!(mag(eigenvalue) < GREAT && mag(eigenvalue) > VSMALL))
    {
        Info<< "    Returning zero magnitude for mag(eigenvalue) = "
            << mag(eigenvalue) << endl;

        return 0;
    }

    // Omit eigenvalue-STDMD step combinations that pose a risk of overflow
    if (mag(eigenvalue)*step_ > sortLimiter_)
    {
        Info<< "    Returning zero magnitude for"
            << " mag(eigenvalue) = " << mag(eigenvalue)
            << " current index = " << step_
            << " sortLimiter = " << sortLimiter_
            << endl;

        return 0;
    }

    scalar magnitude = 0;

    for (label j = 0; j < step_; ++j)
    {
        magnitude += weight[j]*modeNorm*mag(amplitude*pow(eigenvalue, j + 1));
    }

    return magnitude;
}


bool Foam::DMDModels::STDMD::dynamics()
{
    SMatrix Atilde(reducedKoopmanOperator());
    G_.clear();

    if(!eigendecomposition(Atilde))
    {
        return false;
    }

    Atilde.clear();

    frequencies();

    amplitudes();

    magnitudes();

    return true;
}


bool Foam::DMDModels::STDMD::modes()
{
    Info<< tab << "Computing modes" << endl;

    bool processed = false;
    processed = processed || modes<scalar>();
    processed = processed || modes<vector>();
    processed = processed || modes<sphericalTensor>();
    processed = processed || modes<symmTensor>();
    processed = processed || modes<tensor>();

    if (!processed)
    {
        return false;
    }

    return true;
}


void Foam::DMDModels::STDMD::writeToFile(const word& fileName) const
{
    Info<< tab << "Writing objects of dynamics" << endl;

    // Write objects of dynamics
    {
        autoPtr<OFstream> osPtr =
            newFileAtTime
            (
                fileName + "_" + fieldName_,
                mesh_.time().timeOutputValue()
            );
        OFstream& os = osPtr.ref();

        writeFileHeader(os);

        for (const auto& i : labelRange(0, freqs_.size()))
        {
            os  << freqs_[i] << tab
                << mags_[i] << tab
                << amps_[i].real() << tab
                << amps_[i].imag() << tab
                << evals_[i].real() << tab
                << evals_[i].imag() << endl;
        }
    }

    Info<< tab << "Ends output processing for field: " << fieldName_
        << endl;
}


void Foam::DMDModels::STDMD::writeFileHeader(Ostream& os) const
{
    writeHeader(os, "DMD output");
    writeCommented(os, "Frequency");
    writeTabbed(os, "Magnitude");
    writeTabbed(os, "Amplitude (real)");
    writeTabbed(os, "Amplitude (imag)");
    writeTabbed(os, "Eigenvalue (real)");
    writeTabbed(os, "Eigenvalue (imag)");
    os  << endl;
}


void Foam::DMDModels::STDMD::filter()
{
    Info<< tab << "Filtering objects of dynamics" << endl;

    // Filter objects according to magsi
    filterIndexed(evals_, magsi_);
    filterIndexed(evecs_, magsi_);
    filterIndexed(freqs_, magsi_);
    filterIndexed(amps_, magsi_);
    filterIndexed(mags_, magsi_);

    // Clip objects if need be (assuming objects have the same len)
    if (freqs_.size() > nModes_)
    {
        evals_.resize(nModes_);
        evecs_.resize(evecs_.m(), nModes_);
        freqs_.resize(nModes_);
        amps_.resize(nModes_);
        mags_.resize(nModes_);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::DMDModels::STDMD::STDMD
(
    const fvMesh& mesh,
    const word& name,
    const dictionary& dict
)
:
    DMDModel(mesh, name, dict),
    modeSorter_
    (
        modeSorterTypeNames.getOrDefault
        (
            "modeSorter",
            dict,
            modeSorterType::FIRST_SNAPSHOT
        )
    ),
    Q_(),
    G_(),
    Qupper_(1, 1, Zero),
    Qlower_(1, 1, Zero),
    RxInv_(1, 1, Zero),
    evals_(Zero),
    evecs_(1, 1, Zero),
    freqs_(Zero),
    freqsi_(),
    amps_(Zero),
    mags_(Zero),
    magsi_(Zero),
    patches_
    (
        dict.getOrDefault<wordRes>
        (
            "patches",
            dict.found("patch") ? wordRes(1,dict.get<word>("patch")) : wordRes()
        )
    ),
    fieldName_(dict.get<word>("field")),
    timeName0_(),
    minBasis_(0),
    minEval_(0),
    dt_(0),
    sortLimiter_(500),
    fMin_(0),
    fMax_(pTraits<label>::max),
    nGramSchmidt_(5),
    maxRank_(pTraits<label>::max),
    step_(0),
    nModes_(pTraits<label>::max),
    nAgglomerationProcs_(20),
    empty_(false)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::DMDModels::STDMD::read(const dictionary& dict)
{
    writeToFile_ = dict.getOrDefault<bool>("writeToFile", true);

    modeSorter_ =
        modeSorterTypeNames.getOrDefault
        (
            "modeSorter",
            dict,
            modeSorterType::FIRST_SNAPSHOT
        );

    minBasis_ =
        dict.getCheckOrDefault<scalar>("minBasis", 1e-8, scalarMinMax::ge(0));

    minEval_ =
        dict.getCheckOrDefault<scalar>("minEVal", 1e-8, scalarMinMax::ge(0));

    dt_ =
        dict.getCheckOrDefault<scalar>
        (
            "interval",
            (
                dict.getCheck<label>("executeInterval", labelMinMax::ge(1))
               *mesh_.time().deltaTValue()
            ),
            scalarMinMax::ge(0)
        );

    sortLimiter_ =
        dict.getCheckOrDefault<scalar>("sortLimiter", 500, scalarMinMax::ge(1));

    nGramSchmidt_ =
        dict.getCheckOrDefault<label>("nGramSchmidt", 5, labelMinMax::ge(1));

    maxRank_ =
        dict.getCheckOrDefault<label>
        (
            "maxRank",
            pTraits<label>::max,
            labelMinMax::ge(1)
        );

    nModes_ =
        dict.getCheckOrDefault<label>
        (
            "nModes",
            pTraits<label>::max,
            labelMinMax::ge(1)
        );

    fMin_ = dict.getCheckOrDefault<label>("fMin", 0, labelMinMax::ge(0));

    fMax_ =
        dict.getCheckOrDefault<label>
        (
            "fMax",
            pTraits<label>::max,
            labelMinMax::ge(fMin_ + 1)
        );

    nAgglomerationProcs_ =
        dict.getCheckOrDefault<label>
        (
            "nAgglomerationProcs",
            20,
            labelMinMax::ge(1)
        );

    Info<< tab << "Settings are read for:" << nl
        << tab << "    field: " << fieldName_ << nl
        << tab << "    modeSorter: " << modeSorterTypeNames[modeSorter_] << nl
        << tab << "    nModes: " << nModes_ << nl
        << tab << "    maxRank: " << maxRank_ << nl
        << tab << "    nGramSchmidt: " << nGramSchmidt_ << nl
        << tab << "    fMin: " << fMin_ << nl
        << tab << "    fMax: " << fMax_ << nl
        << tab << "    minBasis: " << minBasis_ << nl
        << tab << "    minEVal: " << minEval_ << nl
        << tab << "    sortLimiter: " << sortLimiter_ << nl
        << tab << "    nAgglomerationProcs: " << nAgglomerationProcs_ << nl
        << endl;

    return true;
}


bool Foam::DMDModels::STDMD::initialise(const RMatrix& z)
{
    const scalar norm = L2norm(z);

    if (mag(norm) > 0)
    {
        // First-processed snapshot required by the mode-sorting
        // algorithms at the final output computations (K:p. 43)
        {
            const label nSnap = z.m()/2;
            timeName0_ =
                mesh_.time().timeName(mesh_.time().startTime().value());

            if (nSnap == 0)
            {
                empty_ = true;
            }

            IOField<scalar> snapshot0
            (
                IOobject
                (
                    "snapshot0_" + name_ + "_" + fieldName_,
                    timeName0_,
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                nSnap
            );

            std::copy
            (
                z.cbegin(),
                z.cbegin() + nSnap,
                snapshot0.begin()
            );

            const IOstreamOption streamOpt
            (
                mesh_.time().writeFormat(),
                mesh_.time().writeCompression()
            );

            fileHandler().writeObject(snapshot0, streamOpt, true);
        }

        Q_ = z/norm;
        G_ = SMatrix(1);
        G_(0,0) = sqr(norm);

        ++step_;

        return true;
    }

    return false;
}


bool Foam::DMDModels::STDMD::update(const RMatrix& z)
{
    {
        //- Working copy of the augmented snapshot matrix "z"
        //- being used in the classical Gram-Schmidt process
        const RMatrix ez(orthonormalise(z));

        // Check basis for "z" and, if necessary, expand "Q" and "G"
        const scalar ezNorm = L2norm(ez);
        const scalar zNorm = L2norm(z);

        if (ezNorm/zNorm > minBasis_)
        {
            expand(ez, ezNorm);
        }
    }

    // Update "G" before the potential orthonormal basis compression
    updateG(z);

    // Compress the orthonormal basis if required
    if (Q_.n() >= maxRank_)
    {
        compress();
    }

    ++step_;

    return true;
}


bool Foam::DMDModels::STDMD::fit()
{
    // DMD statistics and mode evaluation (K:Fig. 16)
    const label nSnap = Q_.m()/2;

    // Move upper and lower halves of "Q" to new containers
    Qupper_ = RMatrix(Q_.subMatrix(0, 0, max(nSnap, 1)));
    Qlower_ = RMatrix(Q_.subMatrix(nSnap, 0, max(nSnap, 1)));
    Q_.clear();

    if (!dynamics())
    {
        return true;
    }

    modes();

    if (Pstream::master() && writeToFile_)
    {
        writeToFile(word("dynamics"));

        filter();

        writeToFile(word("filtered_dynamics"));
    }

    step_ = 0;

    return true;
}


// ************************************************************************* //
