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

Class
    Foam::turbulence::IntegralScaleBox

Description
    Class to describe the integral-scale container being used in the
    \c turbulentDigitalFilterInletFvPatchField boundary condition.

SourceFiles
    IntegralScaleBox.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_IntegralScaleBox_H
#define Foam_IntegralScaleBox_H

#include "Random.H"
#include "primitivePatch.H"
#include "boundBox.H"
#include "coordinateSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace turbulence
{

/*---------------------------------------------------------------------------*\
                      Class IntegralScaleBox Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class IntegralScaleBox
{
    // Private Typedefs

        //- Type of input set of integral scales
        typedef typename std::conditional
        <
            std::is_same<Type, scalar>::value,
            vector,
            tensor
        >::type TypeL;


    // Private Enumerations

        //- Options for the correlation function form
        enum kernelType : bool
        {
            GAUSSIAN = true,     //!< "Standard Gaussian correlation function"
            EXPONENTIAL = false  //!< "Exponential correlation function"
        };

        //- Names for kernelType
        static const Enum<kernelType> kernelTypeNames;


    // Private Data

        //- Const reference to patch
        const fvPatch& p_;

        //- Virtual patch of the integral-scale box
        autoPtr<primitivePatch> patchPtr_;

        //- Local coordinate system
        autoPtr<coordinateSystem> csysPtr_;

        //- Correlation function form of the kernel
        const enum kernelType kernelType_;

        //- Random number generator
        Random rndGen_;

        //- Number of faces on integral-scale box inlet plane (<e2> <e3>)
        const Vector2D<label> n_;

        //- Uniform cell size on integral-scale box inlet plane [m]
        Vector2D<scalar> delta_;

        //- Span of patch bounding box in local coordinate system
        vector boundingBoxSpan_;

        //- Min point of patch bounding box in local coordinate system
        vector boundingBoxMin_;

        //- Integral scale set in local coordinates
        const TypeL L_;

        //- Random-number box ranges to be filtered by corresponding kernels
        //  e.g. for U: (Lxu, Lxv, Lxw, Lyu, Lyv, Lyw, Lzu, Lzv, Lzw)
        labelList spans_;

        //- Random-number box
        //  First List: Components of physical variable, e.g. (u,v,w)
        //  Second List: Containing all random numbers in a continuous list
        //  for each component of Cartesian coordinate system, e.g. (x,y,z)
        //  e.g. for U, box_[0] contains sets for u-component in x-y-z order
        scalarListList box_;

        //- Kernel coefficients corresponding to L [-]
        scalarListList kernel_;

        //- Vertices of virtual patch of the integral-scale box
        pointField patchPoints_;

        //- Faces of virtual patch of the integral-scale box
        faceList patchFaces_;

        // Forward-stepwise method information

            //- Flag to enable the forward-stepwise method
            const bool fsm_;

            //- One of the two exponential functions in FSM
            Type C1_;

            //- One of the two exponential functions in FSM
            Type C2_;

            //- Kernel-applied previous-time-step slice field in FSM
            Field<Type> slice_;


    // Private Member Functions

        //- Initialise the local coordinate system with dictionary
        autoPtr<coordinateSystem> calcCoordinateSystem
        (
            const dictionary& dict
        ) const;

        //- Initialise the local coordinate system with patch properties
        void calcCoordinateSystem();

        //- Initialise the bounding box in local coordinates
        Vector2D<vector> calcBoundBox() const;

        //- Initialise the object: delta
        Vector2D<scalar> calcDelta() const;

        //- Initialise the object: spans
        labelList calcSpans() const;

        //- Return the kernel once per simulation
        scalarListList calcKernel() const;

        //- Initialise the object: box
        scalarListList calcBox();

        //- Initialise the object: patchPoints
        pointField calcPatchPoints() const;

        //- Initialise the object: patchFaces
        faceList calcPatchFaces() const;

        //- Initialise the object: virtual patch
        void calcPatch();

        //- Convert integral scales in seconds/meters
        //- to time-step units/integral-scale box inlet cell-size units
        TypeL convert(const TypeL& L) const;

        // Forward-stepwise method

            //- Initialise the object: C1 for scalar fields
            scalar calcC1(const vector& L) const;

            //- Initialise the object: C1 for vector fields
            vector calcC1(const tensor& L) const;

            //- Initialise the object: C2 for scalar fields
            scalar calcC2(const vector& L) const;

            //- Initialise the object: C2 for vector fields
            vector calcC2(const tensor& L) const;

            //- Update C1 and C2 coefficients
            void updateC1C2();


        //- No copy assignment
        void operator=(const IntegralScaleBox&) = delete;


public:

    // Constructors

        //- Construct from patch
        explicit IntegralScaleBox(const fvPatch& p);

        //- Copy construct with patch
        IntegralScaleBox(const fvPatch& p, const IntegralScaleBox& b);

        //- Construct from patch and dictionary
        IntegralScaleBox(const fvPatch& p, const dictionary& dict);

        //- Copy construct
        IntegralScaleBox(const IntegralScaleBox& b);


    // Public Data

        //- Flag to activate debug statements
        static int debug;


    // Member Functions

    // Access

        //- Return const reference to integral-scale box inlet patch
        const primitivePatch& patch()
        {
            if (!patchPtr_)
            {
                calcPatch();
            }

            return *patchPtr_;
        }

        //- Return the object: fsm
        bool fsm() const noexcept
        {
            return fsm_;
        }


    // Evaluation

        //- Initialise integral-scale box properties
        //  Avoid constructor level evaluations
        void initialise();

        //- Discard current time-step integral-scale box slice (the closest
        //- to the patch) by shifting from the back to the front
        void shift();

        //- Add a new integral-scale box slice to the rear of the box
        void refill();

        //- Embed two-point correlations, i.e. L
        //  Apply three-dimensional "valid"-type separable
        //  convolution summation algorithm
        Field<Type> convolve() const;

        // Forward-stepwise method

            //- Apply forward-stepwise correlation for scalar fields
            void correlate(scalarField& fld);

            //- Apply forward-stepwise correlation for vector fields
            void correlate(vectorField& fld);


    // I-O

        //- Write integral-scale box settings
        void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace turbulence
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "IntegralScaleBox.C"
#endif

#endif

// ************************************************************************* //
