/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2018-2022 OpenCFD Ltd.
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
    Foam::triad

Description
    Representation of a 3D Cartesian coordinate system as a Vector of
    row vectors.

See also
    Foam::quaternion

SourceFiles
    triadI.H
    triad.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_triad_H
#define Foam_triad_H

#include "vector.H"
#include "tensor.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class triad;
class quaternion;
inline Istream& operator>>(Istream&, triad&);
inline Ostream& operator<<(Ostream&, const triad&);

/*---------------------------------------------------------------------------*\
                           Class triad Declaration
\*---------------------------------------------------------------------------*/

class triad
:
    public Vector<vector>
{
public:

    // Constructors

        //- Default construct as 'unset'
        inline triad();

        //- Construct from components
        inline triad(const Vector<vector>& vv);

        //- Construct from coordinate axes (row vectors)
        inline triad(const vector& x, const vector& y, const vector& z);

        //- Construct from a tensor
        inline triad(const tensor& t);

        //- Construct from a primary axis with the other two unset
        inline triad(const vector& pa);

        //- Construct from a quaternion
        triad(const quaternion& q);

        //- Construct from Istream
        inline triad(Istream&);


    // Static Data Members

        static const triad I;
        static const triad unset;


    // Member Functions

        //- Is the vector in the direction d set
        inline bool set(const direction d) const;

        //- Are all the vector set
        inline bool set() const;

        //- Return the primary direction of the vector v
        static inline direction primaryDirection(const vector& v);

        //- Return the vector orthogonal to the two provided
        static inline vector orthogonal(const vector& v1, const vector& v2);

        //- Inplace orthogonalise so that it is ortho-normal
        void orthogonalise();

        //- Normalise each set axis vector to have a unit magnitude
        //  Uses vector::normalise with handling of small magnitudes.
        void normalise();

        //- Align this triad with the given vector v
        //  by rotating the most aligned axis to be coincident with v
        void align(const vector& v);

        //- Sort the axes such that they are closest to the x, y and z axes
        triad sortxyz() const;

        //- Convert to a quaternion
        operator quaternion() const;

        //- Return transpose
        inline triad T() const;


    // Column-vector access.

        //- Extract vector for column 0
        inline vector cx() const;

        //- Extract vector for column 1
        inline vector cy() const;

        //- Extract vector for column 2
        inline vector cz() const;


    // Member Operators

        inline void operator=(const Vector<vector>& vv);

        inline void operator=(const tensor& t);

        //- Add the triad t2 to this triad
        //  without normalising or orthogonalising
        void operator+=(const triad& t2);


    // IOstream Operators

        friend Istream& operator>>(Istream&, triad&);
        friend Ostream& operator<<(Ostream&, const triad&);


    // Housekeeping

        //- Same as orthogonalise
        void orthogonalize() { this->orthogonalise(); }

        //- Same as normalise
        void normalize() { this->normalise(); }
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Contiguous data for triad
template<> struct is_contiguous<triad> : std::true_type {};

//- Contiguous 'scalar' data for triad
template<> struct is_contiguous_scalar<triad> : std::true_type {};


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

//- Return a quantity of the difference between two triads
scalar diff(const triad& A, const triad& B);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "triadI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
