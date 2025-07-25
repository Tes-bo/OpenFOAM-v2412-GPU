/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::NURBS3DCurve

Description
    A NURBS 3D curve

SourceFiles
    NURBS3DCurve.C

\*---------------------------------------------------------------------------*/

#ifndef NURBS3DCurve_H
#define NURBS3DCurve_H

#include "NURBSbasis.H"
#include "vectorField.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class NURBS3DCurve Declaration
\*---------------------------------------------------------------------------*/

class NURBS3DCurve
:
    public vectorField
{
    // Private Data

        enum nrmOrientation
        {
            ALIGNED = 1,
            OPPOSED = -1
        };

        //- NURBS basis function
        List<vector> CPs_;
        List<scalar> weights_;
        scalarList u_;
        word name_;
        const NURBSbasis& basis_;

        vector givenInitNrm_;

        label nrmOrientation_;

private:

    // Private Member Functions

        label sgn(const scalar val) const;

        scalar abs(const scalar val) const;

        label mod(const label x, const label interval) const;

        void setUniformU();

        // Bound parametric coor to certain limits.
        bool bound
        (
            scalar& u,
            const scalar minVal,
            const scalar maxVal
        ) const;

        // Function structure assumes that U is from 0->1 and pre-sized to nPts.
        void setEquidistantU
        (
            scalarList& U,
            const label lenAcc,
            const label maxIter,
            const label spacingCorrInterval,
            const scalar tolerance
        ) const;


public:

    // Constructors

        //- Construct from control points, weights and basis order and
        //- parametric coordinates
        NURBS3DCurve
        (
            const NURBSbasis& basis,
            const List<vector>& CPs,
            const List<scalar>& weights,
            const scalarField& u,
            const label nPts,
            const word name="NURBS3DCurve"
        );

        //- Construct from control points, weights and basis order.
        //  Uniform coordinate distribution is implied
        NURBS3DCurve
        (
            const NURBSbasis& basis,
            const List<vector>& CPs,
            const List<scalar>& weights,
            const label nPts,
            const word name="NURBS3DCurve"
        );

        //- Construct from control points and basis order
        //  Uniform coordinate distribution and unit weights are implied
        NURBS3DCurve
        (
            const NURBSbasis& basis,
            const List<vector>& CPs,
            const label nPts,
            const word name="NURBS3DCurve"
        );

        //- Construct as copy
        NURBS3DCurve(const NURBS3DCurve&);


    //- Destructor
    ~NURBS3DCurve() = default;


    // Member Functions

    // Set Functions

        //- Take a given normal and use to determine if
        //- NURBS normals should be reversed.
        //- Computation taken from u = 0.
        void setNrm3DOrientation
        (
            const vector& givenNrm,
            const vector& givenTan
        );

        void setNrm2DOrientation
        (
            const vector& givenNrm,
            const scalar zVal
        );

        //- Flip the orientation of the nrm.
        void flipNrmOrientation();

        //- Set CPs.
        void setCPs(const List<vector>& CPs);

        //- Set weights.
        void setWeights(const List<scalar>& weights);

        //- Set name.
        void setName(const word& name);

        //- Build curve.
        void buildCurve();

        //- Invert CPs order and re-build curve.
        void invert();

        //- Insert a knot by re-computing the control points.
        //  The basis' insertKnot function must have beem called first to
        //  correctly provide the required info.
        void insertKnot
        (
            const scalarField& oldKnots,
            const scalar uBar,
            const label kInsert
        );

        //- Make curve points equidistant in cartesian space.
        void makeEquidistant
        (
            const label lenAcc = 25,
            const label maxIter = 10,
            const label spacingCorrInterval = -1,
            const scalar tolerance = 1.e-5
        );


    // Point Calc Functions

        //- Curve point cartesian coordinates at ui.
        vector curvePoint(const scalar u) const;

        //- Find curve point which is closest to given point
        //- using Newton-Raphson. Returns the param coordinate.
        scalar findClosestCurvePoint
        (
            const vector& targetPoint,
            const label maxIter = 1000,
            const scalar tolerance = 1.e-13
        );

        //- Find curve point which is closest to given point
        //- using Newton-Raphson. Returns the param coordinate.
        scalar findClosestCurvePoint
        (
            const vector& targetPoint,
            const scalar initGuess,
            const label maxIter = 1000,
            const scalar tolerance = 1.e-13
        );

        //- Find the normal to the curve, with the option of forcing a z-plane.
        const vector nrm3D(const vector& refTan, const scalar u) const;

        const vector nrm2D(const scalar zVal, const scalar u) const;

        //- Generate points on the curve which are evenly spaced in cartesian
        //- coordinate distances.
        scalarList genEquidistant
        (
            const label nPts = 100,
            const label lenAcc = 25,
            const label maxIter = 10,
            const label spacingCorrInterval=-1,
            const scalar tolerance = 1.e-5
        );


    // Location Functions

        //- Check if given parametric coordinate u and CP are linked.
        bool checkRange
        (
            const scalar u,
            const label CPI,
            const label degree
        ) const;

        bool checkRange
        (
            const scalar u,
            const label CPI
        ) const;

        //- Calculate Length from starting to ending indices via
        //- computational evaluation using trapezoid rule.
        scalar length
        (
            const label uIStart,
            const label uIEnd
        ) const;

        //- Calculate length from starting to ending parametric coordinates via
        //- computational evaluation using trapezoid rule.
        scalar length
        (
            const scalar uStart,
            const scalar uEnd,
            const label nPts
        ) const;

        //- Calculate length for the entire curve length.
        scalar length() const;


    // Derivative Functions

        //- Curve derivative wrt u at point ui.
        vector curveDerivativeU(const scalar u) const;

        //- Curve second derivative wrt u at point ui.
        vector curveDerivativeUU(const scalar u) const;

        //- Curve derivative wrt b at point ui:
        //- scalar since dx/dX = dy/dY = dz/dZ.
        scalar curveDerivativeCP
        (
            const scalar u,
            const label CPI
        );

        //- Curve derivative wrt CPII at point u.
        vector curveDerivativeWeight
        (
            const scalar u,
            const label CPI
        );

        //- Calculate length from starting to ending indices via
        //- computational evaluation using trapezoid rule
        scalar lengthDerivativeU
        (
            const scalar uStart,
            const scalar uEnd,
            const label nPts
        );


    // Access Functions

        //- Get basis function.
        inline const NURBSbasis& getBasisFunction() const
        {
            return basis_;
        }

        //- Get CPs.
        inline const List<vector>& getCPs() const
        {
            return CPs_;
        }

        //- Get weights.
        inline const List<scalar>& getWeights() const
        {
            return weights_;
        }

        //- Get parametric coordinates.
        inline const List<scalar>& getParametricCoordinates() const
        {
            return u_;
        }

        //- Get name.
        inline const word& getName() const
        {
            return name_;
        }

        //- Return the nrm sgn relation to the u=0 nrm.
        inline label nrmOrientation() const
        {
            return nrmOrientation_;
        }

        //- Return the initial normal given to compare to the Curve's normals.
        inline const vector& givenInitNrm() const
        {
            return givenInitNrm_;
        }


    // Write Functions

        //- Write curve to file.
        void write();

        void write
        (
            const word fileName
        );

        void write
        (
            const fileName dirName,
            const fileName fileName
        );

        void writeWParses();

        void writeWParses
        (
            const word fileName
        );

        void writeWParses
        (
            const fileName dirName,
            const fileName fileName
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
