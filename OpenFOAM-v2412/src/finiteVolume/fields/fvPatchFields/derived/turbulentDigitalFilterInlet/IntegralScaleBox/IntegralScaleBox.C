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

\*---------------------------------------------------------------------------*/

#include "IntegralScaleBox.H"
#include "cartesianCS.H"
#include "OBJstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Type>
const Foam::Enum
<
    typename Foam::turbulence::IntegralScaleBox<Type>::kernelType
>
Foam::turbulence::IntegralScaleBox<Type>::kernelTypeNames
({
    { kernelType::GAUSSIAN, "Gaussian" },
    { kernelType::EXPONENTIAL , "exponential" }
});


template<class Type>
int Foam::turbulence::IntegralScaleBox<Type>::debug = 0;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::coordinateSystem>
Foam::turbulence::IntegralScaleBox<Type>::calcCoordinateSystem
(
    const dictionary& dict
) const
{
    return coordinateSystem::NewIfPresent(dict);
}


template<class Type>
void Foam::turbulence::IntegralScaleBox<Type>::calcCoordinateSystem()
{
    // Get patch normal direction into the domain
    const vector nf((-gAverage(p_.nf())).normalise());

    // Find the second local coordinate direction
    direction minCmpt = 0;
    scalar minMag = mag(nf[minCmpt]);
    for (direction cmpt = 1; cmpt < pTraits<vector>::nComponents; ++cmpt)
    {
        const scalar s = mag(nf[cmpt]);
        if (s < minMag)
        {
            minMag = s;
            minCmpt = cmpt;
        }
    }

    // Create the second local coordinate direction
    vector e2(Zero);
    e2[minCmpt] = 1;

    // Remove normal component
    e2 -= (nf&e2)*nf;

    // Create the local coordinate system - default e3-e1 order
    csysPtr_.reset
    (
        new coordSystem::cartesian
        (
            Zero,   // origin
            nf^e2,  // e3
            nf      // e1
        )
    );
}


template<class Type>
Foam::Vector2D<Foam::vector>
Foam::turbulence::IntegralScaleBox<Type>::calcBoundBox() const
{
    // Convert patch points into local coordinate system
    const pointField localPos
    (
        csysPtr_->localPosition
        (
            pointField
            (
                p_.patch().points(),
                p_.patch().meshPoints()
            )
        )
    );

    // Calculate bounding-box span and min
    const bool globalReduction = true;
    const boundBox bb(localPos, globalReduction);

    return Vector2D<vector>(bb.span(), bb.min());
}


template<class Type>
Foam::Vector2D<Foam::scalar>
Foam::turbulence::IntegralScaleBox<Type>::calcDelta() const
{
    return Vector2D<scalar>
    (
        boundingBoxSpan_[1]/n_.x(),
        boundingBoxSpan_[2]/n_.y()
    );
}


template<class Type>
Foam::labelList Foam::turbulence::IntegralScaleBox<Type>::calcSpans() const
{
    if (!Pstream::master())
    {
        return labelList();
    }

    labelList spans(pTraits<TypeL>::nComponents, label(1));
    const Vector<label> slice(label(1), n_.x(), n_.y());
    const TypeL L(convert(L_));

    label j = 0;
    if (fsm_)
    {
        j = pTraits<Type>::nComponents;
    }

    for (label i = j; i < pTraits<TypeL>::nComponents; ++i)
    {
        const label slicei = label(i/pTraits<Type>::nComponents);

        const label n = ceil(L[i]);
        const label twiceN = 4*n;

        spans[i] = slice[slicei] + twiceN;
    }

    return spans;
}


template<class Type>
Foam::scalarListList
Foam::turbulence::IntegralScaleBox<Type>::calcKernel() const
{
    if (!Pstream::master())
    {
        return scalarListList();
    }

    scalarListList kernel
    (
        pTraits<TypeL>::nComponents,
        scalarList(1, scalar(1))
    );

    const TypeL L(convert(L_));

    label i = 0;
    if (fsm_)
    {
        i = pTraits<Type>::nComponents;
    }

    for (direction dir = i; dir < pTraits<TypeL>::nComponents; ++dir)
    {
        // The smallest kernel width larger than integral scale
        // (KSJ:'n' in Eq. 15)
        const label n = ceil(L[dir]);

        // Full kernel-width (except mid-zero) according to the condition
        // (KSJ:Eq. 15 whereat N is minimum = 2n)
        const label twiceN = 4*n;

        // Initialise kernel-coeffs containers with full kernel-width size
        // Extra elem is mid-zero within [-N, N]
        kernel[dir] = scalarList(twiceN + 1, Zero);

        // First element: -N within [-N, N]
        const scalar initElem = -2*scalar(n);

        // Container initialised with [-N, N] (KSJ:p. 658, item-b)
        std::iota
        (
            kernel[dir].begin(),
            kernel[dir].end(),
            initElem
        );

        // Compute kernel coefficients (KSJ:Eq. 14 (Gaussian))
        scalarList kTemp(kernel[dir]);
        scalar kSum = 0;

        // Model constant shaping the autocorrelation function (KSJ:Eq. 14)
        const scalar C = -0.5*constant::mathematical::pi;

        if (kernelType_)
        {
            const scalar nSqr = n*n;

            kTemp = sqr(Foam::exp(C*sqr(kTemp)/nSqr));
            kSum = Foam::sqrt(sum(kTemp));

            kernel[dir] = Foam::exp(C*sqr(kernel[dir])/nSqr)/kSum;
        }
        else
        {
            kTemp = sqr(Foam::exp(C*mag(kTemp)/n));
            kSum = Foam::sqrt(sum(kTemp));

            kernel[dir] = Foam::exp(C*mag(kernel[dir])/n)/kSum;
        }
    }

    return kernel;
}


template<class Type>
Foam::scalarListList Foam::turbulence::IntegralScaleBox<Type>::calcBox()
{
    if (!Pstream::master())
    {
        return scalarListList();
    }

    scalarListList box(pTraits<Type>::nComponents, scalarList());

    // Initialise: Remaining convenience factors for (e1 e2 e3)
    for (direction dir = 0; dir < pTraits<Type>::nComponents; ++dir)
    {
        scalarList& randomSet = box[dir];

        randomSet = scalarList
        (
            spans_[dir]
           *spans_[dir+pTraits<TypeL>::nComponents/3]
           *spans_[dir+2*(pTraits<TypeL>::nComponents/3)]
        );

        if (randomSet.size() > 1e8)
        {
            WarningInFunction
                << "Size of random-number set is relatively high:" << nl
                << "    size = " << randomSet.size() << nl
                << "    Please consider to use the forward-stepwise method."
                << endl;
        }

        // Initialise: Integral-scale box content with random-number
        // sets obeying the standard normal distribution
        std::generate
        (
            randomSet.begin(),
            randomSet.end(),
            [&]{ return rndGen_.GaussNormal<scalar>(); }
        );
    }

    return box;
}


template<class Type>
Foam::pointField
Foam::turbulence::IntegralScaleBox<Type>::calcPatchPoints() const
{
    if (!Pstream::master())
    {
        return pointField();
    }

    // List of vertex points of the virtual patch in local coordinate system
    const label nx = n_.x();
    const label ny = n_.y();
    const label nPoints = (nx + 1)*(ny + 1);
    pointField points(nPoints, Zero);

    label pointi = 0;
    for (label j = 0; j <= ny; ++j)
    {
        for (label i = 0; i <= nx; ++i)
        {
            const point p
            (
                boundingBoxMin_[0],
                boundingBoxMin_[1] + i*delta_.x(),
                boundingBoxMin_[2] + j*delta_.y()
            );
            points[pointi] = p;
            ++pointi;
        }
    }

    points = csysPtr_->globalPosition(points);

    return points;
}


template<class Type>
Foam::faceList Foam::turbulence::IntegralScaleBox<Type>::calcPatchFaces() const
{
    if (!Pstream::master())
    {
        return faceList();
    }

    // List of faces of the virtual patch
    const label nx = n_.x();
    const label ny = n_.y();
    const label nFaces = nx*ny;
    faceList faces(nFaces);

    label m = 0;
    for (label j = 0; j < ny; ++j)
    {
        for (label i = 0; i < nx; ++i)
        {
            const label k = j*(nx+1) + i;
            faces[m] = face({k, k+(nx+1), k+(nx+2), k+1});
            ++m;
        }
    }

    return faces;
}


template<class Type>
void Foam::turbulence::IntegralScaleBox<Type>::calcPatch()
{
    if (debug && Pstream::master())
    {
        const auto& tm = p_.patch().boundaryMesh().mesh().time();
        OBJstream os(tm.path()/"patch.obj");
        os.write(patchFaces_, patchPoints_, false);
    }

    if (!patchPtr_)
    {
        patchPtr_.reset
        (
            new primitivePatch
            (
                SubList<face>(patchFaces_, patchFaces_.size()),
                patchPoints_
            )
        );
    }
}


template<class Type>
typename Foam::turbulence::IntegralScaleBox<Type>::TypeL
Foam::turbulence::IntegralScaleBox<Type>::convert
(
    const typename Foam::turbulence::IntegralScaleBox<Type>::TypeL& L
) const
{
    TypeL Ls(L);

    const scalar deltaT =
        p_.patch().boundaryMesh().mesh().time().deltaTValue();

    for (direction dir = 0; dir < pTraits<Type>::nComponents; ++dir)
    {
        //  (KSJ:Eq. 13)
        // Integral time scales
        Ls[dir] /= deltaT;
        // Integral length scales
        Ls[dir + Ls.size()/3] /= delta_.x();
        Ls[dir + 2*(Ls.size()/3)] /= delta_.y();
    }

    return Ls;
}


template<class Type>
Foam::scalar Foam::turbulence::IntegralScaleBox<Type>::calcC1
(
    const vector& L
) const
{
    constexpr scalar c1 = -0.25*constant::mathematical::pi;
    return Foam::exp(c1/L.x());
}


template<class Type>
Foam::vector Foam::turbulence::IntegralScaleBox<Type>::calcC1
(
    const tensor& L
) const
{
    constexpr scalar c1 = -0.25*constant::mathematical::pi;

    vector C1(Zero);
    forAll(C1, i)
    {
        C1[i] = Foam::exp(c1/L.x()[i]);
    }

    return C1;
}


template<class Type>
Foam::scalar Foam::turbulence::IntegralScaleBox<Type>::calcC2
(
    const vector& L
) const
{
    constexpr scalar c2 = -0.5*constant::mathematical::pi;
    return Foam::sqrt(scalar(1) - Foam::exp(c2/L.x()));
}


template<class Type>
Foam::vector Foam::turbulence::IntegralScaleBox<Type>::calcC2
(
    const tensor& L
) const
{
    constexpr scalar c2 = -0.5*constant::mathematical::pi;

    vector C2(Zero);
    forAll(C2, i)
    {
        C2[i] = Foam::sqrt(scalar(1) - Foam::exp(c2/L.x()[i]));
    }

    return C2;
}


template<class Type>
void Foam::turbulence::IntegralScaleBox<Type>::updateC1C2()
{
    if (p_.patch().boundaryMesh().mesh().time().isAdjustTimeStep())
    {
        C1_ = calcC1(convert(L_));
        C2_ = calcC2(convert(L_));
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::turbulence::IntegralScaleBox<Type>::IntegralScaleBox
(
    const fvPatch& p
)
:
    p_(p),
    patchPtr_(nullptr),
    csysPtr_(nullptr),
    kernelType_(kernelType::GAUSSIAN),
    rndGen_(0),
    n_(Zero),
    delta_(Zero),
    boundingBoxSpan_(Zero),
    boundingBoxMin_(Zero),
    L_(Zero),
    spans_(Zero),
    box_(Zero),
    kernel_(Zero),
    patchPoints_(Zero),
    patchFaces_(Zero),
    fsm_(false),
    C1_(Zero),
    C2_(Zero),
    slice_(Zero)
{}


template<class Type>
Foam::turbulence::IntegralScaleBox<Type>::IntegralScaleBox
(
    const fvPatch& p,
    const IntegralScaleBox& b
)
:
    p_(p),
    patchPtr_(nullptr),
    csysPtr_(b.csysPtr_.clone()),
    kernelType_(b.kernelType_),
    rndGen_(b.rndGen_),
    n_(b.n_),
    delta_(b.delta_),
    boundingBoxSpan_(b.boundingBoxSpan_),
    boundingBoxMin_(b.boundingBoxMin_),
    L_(b.L_),
    spans_(b.spans_),
    box_(b.box_),
    kernel_(b.kernel_),
    patchPoints_(b.patchPoints_),
    patchFaces_(b.patchFaces_),
    fsm_(b.fsm_),
    C1_(b.C1_),
    C2_(b.C2_),
    slice_(b.slice_)
{}


template<class Type>
Foam::turbulence::IntegralScaleBox<Type>::IntegralScaleBox
(
    const fvPatch& p,
    const dictionary& dict
)
:
    p_(p),
    patchPtr_(nullptr),
    csysPtr_(calcCoordinateSystem(dict)),
    kernelType_
    (
        kernelTypeNames.getOrDefault
        (
            "kernelType",
            dict,
            kernelType::GAUSSIAN
        )
    ),
    rndGen_(time(0)),
    n_(dict.get<Vector2D<label>>("n")),
    delta_(Zero),
    boundingBoxSpan_(Zero),
    boundingBoxMin_(Zero),
    L_(dict.get<TypeL>("L")),
    spans_(Zero),
    box_(Zero),
    kernel_(Zero),
    patchPoints_(Zero),
    patchFaces_(Zero),
    fsm_(dict.getOrDefault("fsm", false)),
    C1_(Zero),
    C2_(Zero),
    slice_(Zero)
{
    if (cmptMin(L_) < ROOTVSMALL)
    {
        FatalIOErrorInFunction(dict)
            << "Integral scale set contains a very small input" << nl
            << "    L = " << L_
            << exit(FatalIOError);
    }

    if (min(n_.x(), n_.y()) <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "Number of faces on box inlet plane has non-positive input"
            << "    n = " << n_
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::turbulence::IntegralScaleBox<Type>::IntegralScaleBox
(
    const IntegralScaleBox& b
)
:
    p_(b.p_),
    patchPtr_(nullptr),
    csysPtr_(b.csysPtr_.clone()),
    kernelType_(b.kernelType_),
    rndGen_(b.rndGen_),
    n_(b.n_),
    delta_(b.delta_),
    boundingBoxSpan_(b.boundingBoxSpan_),
    boundingBoxMin_(b.boundingBoxMin_),
    L_(b.L_),
    spans_(b.spans_),
    box_(b.box_),
    kernel_(b.kernel_),
    patchPoints_(b.patchPoints_),
    patchFaces_(b.patchFaces_),
    fsm_(b.fsm_),
    C1_(b.C1_),
    C2_(b.C2_),
    slice_(b.slice_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::turbulence::IntegralScaleBox<Type>::initialise()
{
    if (!csysPtr_)
    {
        calcCoordinateSystem();
    }

    if (debug && csysPtr_)
    {
        Info<< "Local coordinate system:" << nl
            << "    - origin        = " << csysPtr_->origin() << nl
            << "    - e1-axis       = " << csysPtr_->e1() << nl
            << "    - e2-axis       = " << csysPtr_->e2() << nl
            << "    - e3-axis       = " << csysPtr_->e3() << nl << endl;
    }

    {
        const Vector2D<vector> bb(calcBoundBox());

        boundingBoxSpan_ = bb.x();

        boundingBoxMin_ = bb.y();
    }

    delta_ = calcDelta();

    spans_ = calcSpans();

    kernel_ = calcKernel();

    box_ = calcBox();

    patchPoints_ = calcPatchPoints();

    patchFaces_ = calcPatchFaces();

    calcPatch();

    if (fsm_)
    {
        C1_ = calcC1(convert(L_));

        C2_ = calcC2(convert(L_));

        slice_ = Field<Type>(p_.size(), Zero);
    }
}


template<class Type>
void Foam::turbulence::IntegralScaleBox<Type>::shift()
{
    for (direction dir = 0; dir < pTraits<Type>::nComponents; ++dir)
    {
        scalarList& slice = box_[dir];

        // Slice span: span of each inlet-normal slice of integral-scale box
        // e.g. for U: (Lyu*Lzu, Lyv*Lzv, Lyw*Lzw)
        const label sliceSpan =
            spans_[dir + pTraits<TypeL>::nComponents/3]
           *spans_[dir + 2*(pTraits<TypeL>::nComponents/3)];

        // Shift forward from the back to the front
        inplaceRotateList(slice, sliceSpan);
    }
}


template<class Type>
void Foam::turbulence::IntegralScaleBox<Type>::refill()
{
    for (direction dir = 0; dir < pTraits<Type>::nComponents; ++dir)
    {
        scalarList& slice = box_[dir];

        const label sliceSpan =
            spans_[dir + pTraits<TypeL>::nComponents/3]
           *spans_[dir + 2*(pTraits<TypeL>::nComponents/3)];

        // Refill the back with a new random-number set
        for (label i = 0; i < sliceSpan; ++i)
        {
            slice[i] = rndGen_.GaussNormal<scalar>();
        }
    }
}


template<class Type>
Foam::Field<Type>
Foam::turbulence::IntegralScaleBox<Type>::convolve() const
{
    Field<Type> outFld(n_.x()*n_.y(), Zero);

    for (direction dir = 0; dir < pTraits<Type>::nComponents; ++dir)
    {
        const scalarList& in = box_[dir];

        Field<scalar> out(n_.x()*n_.y(), Zero);

        const scalarList& kernel1 = kernel_[dir];
        const scalarList& kernel2 =
            kernel_[dir + pTraits<TypeL>::nComponents/3];
        const scalarList& kernel3 =
            kernel_[dir + 2*(pTraits<TypeL>::nComponents/3)];

        const label szkernel1 = kernel1.size();
        const label szkernel2 = kernel2.size();
        const label szkernel3 = kernel3.size();

        const label sz1 = spans_[dir];
        const label sz2 = spans_[dir + pTraits<TypeL>::nComponents/3];
        const label sz3 = spans_[dir + 2*(pTraits<TypeL>::nComponents/3)];
        const label sz23 = sz2*sz3;
        const label sz123 = sz1*sz23;

        const label validSlice2 = sz2 - (szkernel2 - 1);
        const label validSlice3 = sz3 - (szkernel3 - 1);

        // Convolution summation - Along 1st direction
        scalarField tmp(sz123, Zero);
        {
            const label filterCentre = label(szkernel2/label(2));
            const label endIndex = sz2 - filterCentre;
            label i0 = 0;
            label i1 = 0;

            for (label i = 0; i < sz1; ++i)
            {
                for (label j = 0; j < sz3; ++j)
                {
                    i1 += filterCentre;

                    for (label k = filterCentre; k < endIndex; ++k)
                    {
                        label q = 0;

                        for (label p = szkernel2 - 1; p >= 0; --p, ++q)
                        {
                            tmp[i1] += in[i0 + q]*kernel2[p];
                        }
                        ++i0;
                        ++i1;
                    }
                    i0 += 2*filterCentre;
                    i1 += filterCentre;
                }
            }
        }

        // Convolution summation - Along 2nd direction
        {
            const scalarField tmp2(tmp);
            const label filterCentre = label(szkernel3/label(2));
            const label endIndex = sz3 - filterCentre;
            label i1 = 0;
            label i2 = 0;

            for (label i = 0; i < sz1; ++i)
            {
                const label sl = i*sz23;

                for (label j = 0; j < sz2; ++j)
                {
                    i1 = j + sl;
                    i2 = i1;

                    for (label k = 0; k < endIndex - filterCentre; ++k)
                    {
                        tmp[i1] = 0;

                        for (label p = szkernel3 - 1, q = 0; p >= 0; --p, ++q)
                        {
                            tmp[i1] += tmp2[i2 + q*sz2]*kernel3[p];
                        }
                        i1 += sz2;
                        i2 += sz2;
                    }
                    i1 += (sz2 + filterCentre);
                    i2 += (sz2 + filterCentre);
                }
            }
        }

        // Convolution summation - Along 3rd direction
        {
            label i1 = (szkernel2 - label(1))/label(2);
            label i2 = (szkernel2 - label(1))/label(2);
            label i3 = 0;

            for (label i = 0; i < validSlice3; ++i)
            {
                for (label j = 0; j < validSlice2; ++j)
                {
                    scalar sum = 0;
                    i1 = i2 + j;

                    for (label k = szkernel1 - 1; k >= 0; --k)
                    {
                        sum += tmp[i1]*kernel1[k];
                        i1 += sz23;
                    }
                    out[i3] = sum;
                    ++i3;
                }
                i2 += sz2;
            }
        }

        outFld.replace(dir, out);
    }

    return outFld;
}


template<class Type>
void Foam::turbulence::IntegralScaleBox<Type>::correlate
(
    scalarField& fld
)
{
    updateC1C2();

    fld *= C2_;
    fld += slice_*C1_;

    // Store current field for the next time-step
    slice_ = fld;
}


template<class Type>
void Foam::turbulence::IntegralScaleBox<Type>::correlate
(
    vectorField& fld
)
{
    updateC1C2();

    for (direction dir = 0; dir < pTraits<vector>::nComponents; ++dir)
    {
        fld.replace
        (
            dir,
            slice_.component(dir)*C1_[dir] + fld.component(dir)*C2_[dir]
        );
    }

    // Store current field for the next time-step
    slice_ = fld;
}


template<class Type>
void Foam::turbulence::IntegralScaleBox<Type>::write
(
    Ostream& os
) const
{
    os.writeEntryIfDifferent<bool>("fsm", false, fsm_);
    os.writeEntry("n", n_);
    os.writeEntry("L", L_);
    os.writeEntry("kernelType", kernelTypeNames[kernelType_]);
    if (csysPtr_)
    {
        csysPtr_->writeEntry(os);
    }
}


// ************************************************************************* //
