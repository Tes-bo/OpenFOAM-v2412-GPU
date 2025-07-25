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
    Foam::NURBS3DSurface

Description
    A NURBS 3D surface

SourceFiles
    NURBS3DSurface.C

\*---------------------------------------------------------------------------*/

#ifndef NURBS3DSurface_H
#define NURBS3DSurface_H

#include "NURBSbasis.H"
#include "vectorField.H"
#include "vector2DField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class NURBS3DSurface Declaration
\*---------------------------------------------------------------------------*/

class NURBS3DSurface
:
    public vectorField
{
    // Private Data

        enum paramType
        {
            PARAMU=0,
            PARAMV,
            NPARAMS
        };

        enum nrmOrientation
        {
            ALIGNED=1,
            OPPOSED=-1
        };

        //- Held in sets of const v.
        List<vector> CPs_;

        //- Held in sets of const u.
        scalarList u_;

        scalarList v_;

        scalarList weights_;

        label nUPts_;

        label nVPts_;

        word name_;

        NURBSbasis uBasis_;

        NURBSbasis vBasis_;

        vector givenInitNrm_;

        labelList CPsUCPIs_;

        labelList CPsVCPIs_;

        label nrmOrientation_;

        autoPtr<labelList> boundaryCPIDs_;

        autoPtr<labelList> whichBoundaryCPID_;


    // Private Member Functions

        label sgn(const scalar val) const;

        scalar abs(const scalar val) const;

        label mod(const label x, const label interval) const;

        void setCPUVLinking();

        void setUniformUV
        (
            scalarList& u,
            scalarList& v,
            const label nUPts,
            const label nVPts
        ) const;

        void setUniformUV();

        bool bound
        (
            scalar& u,
            scalar& v,
            const scalar minVal=1e-7,
            const scalar maxVal=0.999999
        ) const;

        bool boundDirection
        (
            scalar& u,
            const scalar minVal=1e-7,
            const scalar maxVal=0.999999
        ) const;

        // Function structure assumes that range specified for param R
        // should range from 0->1 and is sized correctly for nPts.
        // Pts are set in equidistant manner along fixed value of other param S.
        void setEquidistantR
        (
            scalarList& R,
            const scalar SHeld,
            const label paramR,
            const label lenAcc,
            const label maxIter,
            const label spacingCorrInterval,
            const scalar tolerance
        ) const;

public:

    // Constructors

        //- Construct from number of control points and basis functions
        NURBS3DSurface
        (
            const List<vector>& CPs,
            const label nPointsU,
            const label nPointsV,
            const NURBSbasis& uBasis,
            const NURBSbasis& vBasis,
            const word name="NURBS3DSurface"
        );

        //- Construct from number of control points, weights and basis functions
        NURBS3DSurface
        (
            const List<vector>& CPs,
            const List<scalar>& weights,
            const label nPointsU,
            const label nPointsV,
            const NURBSbasis& uBasis,
            const NURBSbasis& vBasis,
            const word name="NURBS3DSurface"
        );

        //- Construct from control points,  basis degree and number of points
        NURBS3DSurface
        (
            const List<vector>& CPs,
            const label nPointsU,
            const label nPointsV,
            const label uDegree,
            const label vDegree,
            const label nCPsU,
            const label nCPsV,
            const word name="NURBS3DSurface"
        );

        //- Construct from control points,  basis degree, knots and number of
        //  points
        NURBS3DSurface
        (
            const List<vector>& CPs,
            const label nPointsU,
            const label nPointsV,
            const label uDegree,
            const label vDegree,
            const label nCPsU,
            const label nCPsV,
            const scalarField& knotsU,
            const scalarField& knotsV,
            const word name="NURBS3DSurface"
        );

        //- Construct from control points, weights, basis degree, and number of
        //  points
        NURBS3DSurface
        (
            const List<vector>& CPs,
            const List<scalar>& weights,
            const label nPointsU,
            const label nPointsV,
            const label uDegree,
            const label vDegree,
            const label nCPsU,
            const label nCPsV,
            const word name="NURBS3DSurface"
        );

        //- Construct from control points, weights, basis degree, knots and
        //  number of points
        NURBS3DSurface
        (
            const List<vector>& CPs,
            const List<scalar>& weights,
            const label nPointsU,
            const label nPointsV,
            const label uDegree,
            const label vDegree,
            const label nCPsU,
            const label nCPsV,
            const scalarField& knotsU,
            const scalarField& knotsV,
            const word name="NURBS3DSurface"
        );

        //- Construct as copy
        NURBS3DSurface(const NURBS3DSurface&);


    //- Destructor

        ~NURBS3DSurface() = default;

    // Member Functions

        // Set Functions

            void setNrmOrientation
            (
                const vector& givenNrm,
                const scalar u,
                const scalar v
            );

            //- Flip the orientation of the nrm.
            void flipNrmOrientation();

            void setCPs(const List<vector>& CPs);

            void setWeights(const scalarList& weights);

            void setName(const word& name);

            void buildSurface();

            void invertU();

            void invertV();

            void invertUV();

            void makeEquidistant
            (
                const label lenAcc=25,
                const label maxIter=10,
                const label spacingCorrInterval=-1,
                const scalar tolerance=1.e-5
            );

    // Calculation Functions

        vector surfacePoint( const scalar& u, const scalar& v);

        scalarList findClosestSurfacePoint
        (
            const vector& targetPoint,
            const label maxIter=100,
            const scalar tolerance=1.e-6
        );

        tmp<vector2DField> findClosestSurfacePoint
        (
            const vectorField& targetPoints,
            const label maxIter=100,
            const scalar tolerance=1.e-6
        );

        scalarList findClosestSurfacePoint
        (
            const vector& targetPoint,
            const scalar& uInitGuess,
            const scalar& vInitGuess,
            const label maxIter=100,
            const scalar tolerance=1.e-6
        );

        const vector nrm(scalar u, scalar v);

        //- Generate points on the surface which are evenly spaced in cartesian
        //  coordinate distances.
        List<scalarList> genEquidistant
        (
            const label nUPts=100,
            const label nVPts=100,
            const label lenAcc=25,
            const label maxIter=10,
            const label spacingCorrInterval=-1,
            const scalar tolerance=1.e-5
        );

    // Location Functions

        bool checkRangeU
        (
            const scalar u,
            const label CPI,
            const label uDegree
        ) const;

        bool checkRangeU(const scalar u, const label CPI) const;

        bool checkRangeV
        (
            const scalar v,
            const label CPI,
            const label vDegree
        ) const;

        bool checkRangeV(const scalar v, const label CPI) const;

        bool checkRangeUV
        (
            const scalar v,
            const scalar u,
            const label CPI,
            const label uDegree,
            const label vDegree
        ) const;

        bool checkRangeUV
        (
            const scalar v,
            const scalar u,
            const label CPI
        ) const;

        scalar lengthU
        (
            const label vIConst,
            const label uIStart,
            const label uIEnd
        ) const;

        scalar lengthU
        (
            const scalar vConst,
            const scalar uStart,
            const scalar uEnd,
            const label nPts
        ) const;

        scalar lengthU(const label vIConst) const;

        scalar lengthU(const scalar vConst) const;

        scalar lengthV
        (
            const label uIConst,
            const label vIStart,
            const label vIEnd
        ) const;

        scalar lengthV
        (
            const scalar uConst,
            const scalar vStart,
            const scalar vEnd,
            const label nPts
        ) const;

        scalar lengthV(const label uIConst) const;

        scalar lengthV(const scalar uConst) const;

    // Derivative Functions

        //- Surface derivative wrt u at point u,v.
        vector surfaceDerivativeU(const scalar u, const scalar v) const;

        //- Surface derivative wrt v at point u,v.
        vector surfaceDerivativeV(const scalar u, const scalar v) const;

        //- Surface second derivative wrt u and v at point u,v.
        vector surfaceDerivativeUV(const scalar u, const scalar v) const;

        //- Surface second derivative wrt u at point u,v.
        vector surfaceDerivativeUU(const scalar u, const scalar v) const;

        //- Surface second derivative wrt v at point u,v.
        vector surfaceDerivativeVV(const scalar u, const scalar v) const;

        //- Surface derivative wrt the weight of CPI at point u,v.
        scalar surfaceDerivativeCP
        (
            const scalar u,
            const scalar v,
            const label CPI
        ) const;

        //- Surface derivative wrt WI at point u,v.
        vector surfaceDerivativeW
        (
            const scalar u,
            const scalar v,
            const label CPI
        ) const;

        //- Surface derivative wrt u length along v=const contour range.
        scalar lengthDerivativeU
        (
            const scalar vConst,
            const scalar uStart,
            const scalar uEnd,
            const label nPts
        ) const;

        //- Surface derivative wrt v length along u=const contour range.
        scalar lengthDerivativeV
        (
            const scalar uConst,
            const scalar vStart,
            const scalar vEnd,
            const label nPts
        ) const;


    // Access Functions

        //- Get basis function.
        inline const NURBSbasis& getBasisFunctionU() const
        {
            return uBasis_;
        }

        inline const NURBSbasis& getBasisFunctionV() const
        {
            return vBasis_;
        }

        //- Get CPs.
        inline const List<vector>& getCPs() const
        {
            return CPs_;
        }

        //- Get weights.
        inline const scalarList& getWeights() const
        {
            return weights_;
        }

        //- Get parametric coordinates.
        inline const scalarList& getParametricCoordinatesU() const
        {
            return u_;
        }

        inline const scalarList& getParametricCoordinatesV() const
        {
            return v_;
        }

        //- Get name.
        inline const word& getName() const
        {
            return name_;
        }

        //- Get number of point in u direction
        inline label getNPtsU() const
        {
            return nUPts_;
        }

        //- Get number of point in u direction
        inline label getNPtsV() const
        {
            return nVPts_;
        }

        //- Get IDs of boundary control points
        const labelList& getBoundaryCPIDs();

        const labelList& getBoundaryCPIs();

        //- Get the boundary CP ID based on the global CP ID
        const label& whichBoundaryCPI(const label& globalCPI);

        //- Return the nrm sgn relation to the u=0 nrm.
        inline label nrmOrientation() const
        {
            return nrmOrientation_;
        }

        //- Return the initial nrmal given to compare to the Curve's nrmals.
        inline const vector& givenInitNrm() const
        {
            return givenInitNrm_;
        }

        //- Return ID in u direction for a given cp ID
        const labelList& getCPsUCPIs() const
        {
            return CPsUCPIs_;
        }

        //- Return ID in v direction for a given cp ID
        const labelList& getCPsVCPIs() const
        {
            return CPsVCPIs_;
        }

    // Write Functions

        //- Write curve to file
        void write();

        void write(const word fileName);

        void write(const fileName dirName, const fileName fileName);

        void writeWParses();

        void writeWParses(const word fileName);

        void writeWParses(const fileName dirName, const fileName fileName);

        void writeVTK(const fileName vtkDirName, const fileName vtkFileName);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
