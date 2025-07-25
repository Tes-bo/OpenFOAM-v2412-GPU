/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

Class
    Foam::interpolation2DTable

Description
    2D table interpolation.
    The data must be in ascending order in both dimensions x and y.

SourceFiles
    interpolation2DTable.C

\*---------------------------------------------------------------------------*/

#ifndef interpolation2DTable_H
#define interpolation2DTable_H

#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class interpolation2DTable Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class interpolation2DTable
:
    public List<Tuple2<scalar, List<Tuple2<scalar, Type>>>>
{
    // Private Data

        //- Handling for out-of-bound values
        bounds::normalBounding bounding_;

        //- File name
        fileName fileName_;

        //- Table reader
        autoPtr<tableReader<Type>> reader_;


    // Private Member Functions

        //- Read the table of data from file
        void readTable();

        //- Interpolated value within 1-D list
        Type interpolateValue
        (
            const List<Tuple2<scalar, Type>>& list,
            scalar lookupValue
        ) const;

        //- Return an X index from the matrix
        template<class BinaryOp>
        label Xi
        (
            const BinaryOp& bop,
            const scalar valueX,
            const bool reverse
        ) const;


public:

    // Public Data Types

        //- The element data type
        typedef Tuple2<scalar, List<Tuple2<scalar, Type>>> value_type;

        //- Convenience typedef
        typedef List<Tuple2<scalar, List<Tuple2<scalar, Type>>>> table;


    // Constructors

        //- Default construct
        interpolation2DTable();

        //- Construct from components
        interpolation2DTable
        (
            const List<Tuple2<scalar, List<Tuple2<scalar, Type>>>>& values,
            const bounds::normalBounding bounding,
            const fileName& fName
        );

        //- Construct given the name of the file containing the table of data
        explicit interpolation2DTable(const fileName& fName);

        //- Construct by reading file name and outOfBounds from dictionary
        explicit interpolation2DTable(const dictionary& dict);

        //- Copy construct
        interpolation2DTable(const interpolation2DTable& tbl);


    // Member Functions

        //- Check that list is monotonically increasing
        //  Exit with a FatalError if there is a problem
        void check() const;

        //- Write
        void write(Ostream& os) const;


    // Member Operators

        //- Copy assignment
        void operator=(const interpolation2DTable<Type>& rhs);

        //- Return an interpolated value
        Type operator()(const scalar valueX, const scalar valueY) const;


    // Housekeeping

        //- Deprecated(2019-08) check list is monotonically increasing
        //  \deprecated(2019-08) - older name for check() method
        void checkOrder() const { check(); }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "interpolation2DTable.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
