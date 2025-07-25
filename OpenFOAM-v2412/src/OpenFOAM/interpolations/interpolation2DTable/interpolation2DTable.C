/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "openFoamTableReader.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class Type>
void Foam::interpolation2DTable<Type>::readTable()
{
    fileName fName(fileName_);
    fName.expand();

    // Read data from file
    reader_()(fName, *this);

    if (this->empty())
    {
        FatalErrorInFunction
            << "table read from " << fName << " is empty" << nl
            << exit(FatalError);
    }

    // Check that the data are in ascending order
    check();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::interpolation2DTable<Type>::interpolation2DTable()
:
    List<value_type>(),
    bounding_(bounds::normalBounding::WARN),
    fileName_("fileNameIsUndefined"),
    reader_(nullptr)
{}


template<class Type>
Foam::interpolation2DTable<Type>::interpolation2DTable
(
    const List<Tuple2<scalar, List<Tuple2<scalar, Type>>>>& values,
    const bounds::normalBounding bounding,
    const fileName& fName
)
:
    List<value_type>(values),
    bounding_(bounding),
    fileName_(fName),
    reader_(nullptr)
{}


template<class Type>
Foam::interpolation2DTable<Type>::interpolation2DTable(const fileName& fName)
:
    List<value_type>(),
    bounding_(bounds::normalBounding::WARN),
    fileName_(fName),
    reader_(new openFoamTableReader<Type>())
{
    readTable();
}


template<class Type>
Foam::interpolation2DTable<Type>::interpolation2DTable(const dictionary& dict)
:
    List<value_type>(),
    bounding_
    (
        bounds::normalBoundingNames.getOrDefault
        (
            "outOfBounds",
            dict,
            bounds::normalBounding::WARN,
            true  // Failsafe behaviour
        )
    ),
    fileName_(dict.get<fileName>("file")),
    reader_(tableReader<Type>::New(dict))
{
    readTable();
}


template<class Type>
Foam::interpolation2DTable<Type>::interpolation2DTable
(
     const interpolation2DTable& tbl
)
:
    List<value_type>(tbl),
    bounding_(tbl.bounding_),
    fileName_(tbl.fileName_),
    reader_(tbl.reader_.clone())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::interpolation2DTable<Type>::interpolateValue
(
    const List<Tuple2<scalar, Type>>& list,
    scalar lookupValue
) const
{
    return interpolationTable<Type>::interpolateValue
    (
        list,
        lookupValue,
        bounds::repeatableBounding(bounding_)
    );
}


template<class Type>
template<class BinaryOp>
Foam::label Foam::interpolation2DTable<Type>::Xi
(
    const BinaryOp& bop,
    const scalar valueX,
    const bool reverse
) const
{
    const List<value_type>& t = *this;

    label limitI = 0;
    if (reverse)
    {
        limitI = t.size() - 1;
    }

    if (bop(valueX, t[limitI].first()))
    {
        switch (bounding_)
        {
            case bounds::normalBounding::ERROR:
            {
                FatalErrorInFunction
                    << "value (" << valueX << ") out of bounds" << nl
                    << exit(FatalError);
                break;
            }
            case bounds::normalBounding::WARN:
            {
                WarningInFunction
                    << "value (" << valueX << ") out of bounds" << nl;

                // Behaviour as per CLAMP
                return limitI;
                break;
            }
            case bounds::normalBounding::CLAMP:
            {
                return limitI;
                break;
            }
            default:
            {
                FatalErrorInFunction
                    << "Unhandled bounding type " << int(bounding_)
                    << abort(FatalError);
            }
        }
    }

    label i = 0;
    if (reverse)
    {
        const label nX = t.size();
        i = 0;
        while ((i < nX) && (valueX > t[i].first()))
        {
            ++i;
        }
    }
    else
    {
        i = t.size() - 1;
        while ((i > 0) && (valueX < t[i].first()))
        {
            --i;
        }
    }

    return i;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::interpolation2DTable<Type>::operator=
(
    const interpolation2DTable<Type>& rhs
)
{
    if (this == &rhs)
    {
        return;
    }

    static_cast<List<value_type>&>(*this) = rhs;
    bounding_ = rhs.bounding_;
    fileName_ = rhs.fileName_;
    reader_.reset(rhs.reader_.clone());
}


template<class Type>
Type Foam::interpolation2DTable<Type>::operator()
(
    const scalar valueX,
    const scalar valueY
) const
{
    const List<value_type>& t = *this;

    // Assumes all of the list in Y being equal length
    const label nX = t.size();

    if (nX == 0)
    {
        WarningInFunction
            << "Cannot interpolate zero-sized table - returning zero" << nl;

        return Zero;
    }
    else if (nX == 1)
    {
        // Only 1 column (in X) - simply interpolate to find Y value
        return interpolateValue(t.first().second(), valueY);
    }


    // Find low and high indices in the X range that bound valueX
    const label lo = Xi(lessOp<scalar>(), valueX, false);
    const label hi = Xi(greaterOp<scalar>(), valueX, true);

    if (lo == hi)
    {
        return interpolateValue(t[lo].second(), valueY);
    }


    // Normal interpolation

    const Type y0(interpolateValue(t[lo].second(), valueY));
    const Type y1(interpolateValue(t[hi].second(), valueY));

    const scalar& x0 = t[lo].first();
    const scalar& x1 = t[hi].first();

    return (y0 + (y1 - y0)*(valueX - x0)/(x1 - x0));
}


template<class Type>
void Foam::interpolation2DTable<Type>::check() const
{
    const List<value_type>& list = *this;

    scalar prevValue(0);

    label i = 0;
    for (const auto& item : list)
    {
        const scalar& currValue = item.first();

        // Avoid duplicate values (divide-by-zero error)
        if (i && currValue <= prevValue)
        {
            FatalErrorInFunction
                << "out-of-order value: "
                << currValue << " at index " << i << nl
                << exit(FatalError);
        }
        prevValue = currValue;
        ++i;
    }
}


template<class Type>
void Foam::interpolation2DTable<Type>::write(Ostream& os) const
{
    os.writeEntry("file", fileName_);
    os.writeEntry("outOfBounds", bounds::normalBoundingNames[bounding_]);

    os  << *this;
}


// ************************************************************************* //
