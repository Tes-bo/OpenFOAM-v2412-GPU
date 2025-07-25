/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020-2021 OpenCFD Ltd.
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
    Foam::Function1Types::TableBase

Description
    Base class for table with bounds handling, interpolation and integration

SourceFiles
    TableBase.C

\*---------------------------------------------------------------------------*/

#ifndef Function1Types_TableBase_H
#define Function1Types_TableBase_H

#include "tableBounds.H"
#include "Function1.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class interpolationWeights;

namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                        Class TableBase Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class TableBase
:
    public Function1<Type>
{
protected:

    // Protected Data

        //- Handling for out-of-bound values
        const bounds::repeatableBounding bounding_;

        //- Interpolation type
        const word interpolationScheme_;

        //- Table data
        List<Tuple2<scalar, Type>> table_;

        //- Extracted values
        mutable autoPtr<scalarField> tableSamplesPtr_;

        //- Interpolator method
        mutable autoPtr<interpolationWeights> interpolatorPtr_;

        //- Cached indices and weights
        mutable labelList currentIndices_;

        mutable scalarField currentWeights_;


    // Protected Member Functions

        //- Return (demand driven) interpolator
        const interpolationWeights& interpolator() const;

        //- No copy assignment
        void operator=(const TableBase<Type>&) = delete;


public:

    // Constructors

        //- Construct from dictionary - note table is not populated
        TableBase
        (
            const word& name,
            const dictionary& dict,
            const objectRegistry* obrPtr = nullptr
        );

        //- Copy constructor. Note: steals interpolator, tableSamples
        explicit TableBase(const TableBase<Type>& tbl);


    //- Destructor
    virtual ~TableBase() = default;


    // Member Functions

        //- Check the table for size and consistency
        void initialise();

        //- Check minimum table bounds
        bool checkMinBounds(const scalar x, scalar& xDash) const;

        //- Check maximum table bounds
        bool checkMaxBounds(const scalar x, scalar& xDash) const;

        //- Convert time
        virtual void userTimeToTime(const Time& t);

        //- Return Table value
        virtual Type value(const scalar x) const;

        //- Integrate between two (scalar) values
        virtual Type integrate(const scalar x1, const scalar x2) const;

        //- Return the reference values
        virtual tmp<scalarField> x() const;

        //- Return the dependent values
        virtual tmp<Field<Type>> y() const;

        //- Write all table data in dictionary format
        virtual void writeData(Ostream& os) const;

        //- Write keywords only in dictionary format.
        //  Used for non-inline table types
        virtual void writeEntries(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "TableBase.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
