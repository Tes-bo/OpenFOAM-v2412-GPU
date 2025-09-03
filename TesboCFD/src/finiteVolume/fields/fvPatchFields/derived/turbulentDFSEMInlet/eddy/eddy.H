/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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
    Foam::eddy

Description
    Class to describe eddies for the turbulentDFSEMInletFvPatchVectorField
    boundary condition.

SourceFiles
    eddy.C
    eddyI.H
    eddyIO.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentDFSEMInletFvPatchVectorField_eddy_H
#define turbulentDFSEMInletFvPatchVectorField_eddy_H

#include "vector.H"
#include "point.H"
#include "tensor.H"
#include "Random.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class eddy;

bool operator==(const eddy& a, const eddy& b);
bool operator!=(const eddy& a, const eddy& b);
Istream& operator>>(Istream& is, eddy& e);
Ostream& operator<<(Ostream& os, const eddy& e);


/*---------------------------------------------------------------------------*\
                            Class eddy Declaration
\*---------------------------------------------------------------------------*/

class eddy
{
    // Private Data

        static label Gamma2Values[8];
        static UList<label> Gamma2;

        //- Patch face index that spawned the eddy
        label patchFaceI_;

        //- Reference position
        point position0_;

        //- Distance from reference position in normal direction
        scalar x_;

        //- Integral-length scales in 3-D space
        vector sigma_;

        //- Time-averaged intensity
        vector alpha_;

        //- Coordinate system transformation from local to global axes
        //  X-direction aligned with max stress eigenvalue
        tensor Rpg_;

        //- Model coefficient c1
        scalar c1_;

        //- Index of streamwise direction (0,1,2)
        label dir1_;


    // Private Member Functions

        //- Set the eddy scales: length, intensity
        bool setScales
        (
            const scalar sigmaX,
            const label gamma2,
            const vector& e,
            const vector& lambda,
            vector& sigma,
            vector& alpha
        ) const;

        //- Return a number with zero mean and unit variance
        inline scalar epsi(Random& rndGen) const;


public:

    // Constructors

        //- Construct null
        eddy();

        //- Construct from Istream
        eddy(Istream& is);

        //- Construct from components
        eddy
        (
            const label patchFaceI,     // patch face index
            const point& position0,     // reference position
            const scalar x,             // distance from reference position
            const scalar sigmaX,        // integral-length scale
            const symmTensor& R,        // Reynolds stress tensor
            Random& rndGen
        );

        //- Copy construct
        eddy(const eddy& e);


    // Public Data

        //- Flag to activate debug statements
        static int debug;


    // Public Member Functions

        // Access

            //- Return the patch face index that spawned the eddy
            inline label patchFaceI() const noexcept;

            //- Return the reference position
            inline const point& position0() const noexcept;

            //- Return the distance from the reference position
            inline scalar x() const noexcept;

            //- Return the length scales in 3-D space
            inline const vector& sigma() const noexcept;

            //- Return the time-averaged intensity
            inline const vector& alpha() const noexcept;

            //- Return the coordinate system transformation from local
            //- principal to global axes
            inline const tensor& Rpg() const noexcept;

            //- Return the model coefficient c1
            inline scalar c1() const noexcept;

            //- Return the eddy position
            inline point position(const vector& n) const;

            //- Return the index of the streamwise direction (0,1,2)
            label dir1() const noexcept { return dir1_; }

            //- Return random vector of -1 and 1's
            inline vector epsilon(Random& rndGen) const;


        // Helper functions

            //- Volume
            inline scalar volume() const;

            //- Move the eddy
            inline void move(const scalar dx);

            //- Eddy bounds
            inline boundBox bounds(const bool global = true) const;


        // Evaluation

            //- Return the fluctuating velocity contribution at local point xp
            vector uPrime(const point& xp, const vector& n) const;


        // Writing

            //- Write the eddy centre in OBJ format
            void writeCentreOBJ(const vector& n, Ostream& os) const;

            //- Write the eddy surface in OBJ format
            //  Returns the number of points used to describe the eddy surface
            label writeSurfaceOBJ
            (
                const label pointOffset,
                const vector& n,
                Ostream& os
            ) const;


    // Member Operators

        void operator=(const eddy& e);


    // Friend Operators

        friend bool operator==(const eddy& a, const eddy& b)
        {
            return
                a.patchFaceI_ == b.patchFaceI_
             && a.position0_ == b.position0_
             && a.x_ == b.x_
             && a.sigma_ == b.sigma_
             && a.alpha_ == b.alpha_
             && a.Rpg_ == b.Rpg_
             && a.c1_ == b.c1_
             && a.dir1_ == b.dir1_;
        }

        friend bool operator!=(const eddy& a, const eddy& b)
        {
            return !(a == b);
        }


    // IOstream Operators

        friend Istream& operator>>(Istream& is, eddy& e);
        friend Ostream& operator<<(Ostream& os, const eddy& e);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "eddyI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
