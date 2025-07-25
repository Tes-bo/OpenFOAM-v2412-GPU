/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::labelledTri

Description
    A triFace with additional (region) index.

SourceFiles
    labelledTriI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_labelledTri_H
#define Foam_labelledTri_H

#include "triFace.H"
#include "ListListOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class labelledTri;
inline Istream& operator>>(Istream&, labelledTri&);
inline Ostream& operator<<(Ostream&, const labelledTri&);

/*---------------------------------------------------------------------------*\
                         Class labelledTri Declaration
\*---------------------------------------------------------------------------*/

class labelledTri
:
    public triFace
{
    // Private Data

        //- The object index (region)
        label index_;


    // Private Member Functions

        //- Assign from 3 or 4 labels, default region is 0.
        template<class ListType>
        inline void assignList(const ListType& list);


public:

    // Constructors

        //- Default construct, with invalid point labels and region (-1).
        inline labelledTri();

        //- Construct from triFace
        //- and optional region index (0 if unspecified)
        inline labelledTri
        (
            const triFace& tri,
            const label region = 0
        );

        //- Construct from three point labels
        //- and optional region index (0 if unspecified)
        inline labelledTri
        (
            const label p0,
            const label p1,
            const label p2,
            const label region = 0
        );

        //- Construct from a list of 3 or 4 labels. Default region is 0.
        inline explicit labelledTri(const labelUList& list);

        //- Construct from a list of 3 or 4 labels. Default region is 0.
        inline explicit labelledTri(std::initializer_list<label>);

        //- Copy construct from a subset of point labels
        //- and optional region index (0 if unspecified)
        inline labelledTri
        (
            const labelUList& list,
            const FixedList<label, 3>& triIndices,
            const label region = 0
        );

        //- Copy construct from a subset of point labels
        //- with region index from input labelledTri
        inline labelledTri
        (
            const labelUList& list,
            const labelledTri& triIndices
        );

        //- Construct from Istream
        inline labelledTri(Istream& is);


    // Member Functions

        //- Return the index (eg, the region)
        label index() const noexcept
        {
            return index_;
        }

        //- Non-const access to the index (eg, the region)
        label& index() noexcept
        {
            return index_;
        }

        //- Set the index (eg, the region)
        void setIndex(const label idx) noexcept
        {
            index_ = idx;
        }

        //- Return the region index
        label region() const noexcept
        {
            return index_;
        }

        //- Non-const access to the region index
        label& region() noexcept
        {
            return index_;
        }


    // IOstream Operators

        friend Istream& operator>>(Istream&, labelledTri&);
        friend Ostream& operator<<(Ostream&, const labelledTri&);
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Contiguous data for labelledTri
template<> struct is_contiguous<labelledTri> : std::true_type {};

//- Contiguous label data for labelledTri
template<> struct is_contiguous_label<labelledTri> : std::true_type {};


//- Specialization to offset faces, used in ListListOps::combineOffset
template<>
struct offsetOp<labelledTri>
{
    labelledTri operator()(const labelledTri& x, const label offset) const
    {
        labelledTri f(x);
        f += offset;
        return f;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "labelledTriI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
