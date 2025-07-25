/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "TableBase.H"
#include "Time.H"
#include "interpolationWeights.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
const Foam::interpolationWeights&
Foam::Function1Types::TableBase<Type>::interpolator() const
{
    if (!interpolatorPtr_)
    {
        // Re-work table into linear list
        tableSamplesPtr_.reset(new scalarField(table_.size()));
        auto& samples = *tableSamplesPtr_;
        forAll(table_, i)
        {
            samples[i] = table_[i].first();
        }
        interpolatorPtr_ = interpolationWeights::New
        (
            interpolationScheme_,
            samples
        );
    }

    return *interpolatorPtr_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::TableBase<Type>::TableBase
(
    const word& name,
    const dictionary& dict,
    const objectRegistry* obrPtr
)
:
    Function1<Type>(name, dict, obrPtr),
    bounding_
    (
        bounds::repeatableBoundingNames.getOrDefault
        (
            "outOfBounds",
            dict,
            bounds::repeatableBounding::CLAMP,
            true  // Failsafe behaviour
        )
    ),
    interpolationScheme_
    (
        dict.getOrDefault<word>
        (
            "interpolationScheme",
            "linear",
            keyType::LITERAL
        )
    ),
    table_(),
    tableSamplesPtr_(nullptr),
    interpolatorPtr_(nullptr)
{}


template<class Type>
Foam::Function1Types::TableBase<Type>::TableBase(const TableBase<Type>& tbl)
:
    Function1<Type>(tbl),
    bounding_(tbl.bounding_),
    interpolationScheme_(tbl.interpolationScheme_),
    table_(tbl.table_),
    tableSamplesPtr_(nullptr),
    interpolatorPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::Function1Types::TableBase<Type>::initialise()
{
    if (!table_.size())
    {
        FatalErrorInFunction
            << "Table for entry " << this->name() << " is invalid (empty)"
            << nl << exit(FatalError);
    }

    scalar prevValue(0);

    label i = 0;
    for (const auto& item : table_)
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
bool Foam::Function1Types::TableBase<Type>::checkMinBounds
(
    const scalar x,
    scalar& xDash
) const
{
    const scalar minLimit = table_.first().first();
    const scalar maxLimit = table_.last().first();

    if (x < minLimit)
    {
        switch (bounding_)
        {
            case bounds::repeatableBounding::ERROR:
            {
                FatalErrorInFunction
                    << "value (" << x << ") less than lower "
                    << "bound (" << minLimit << ")" << nl
                    << exit(FatalError);
                break;
            }
            case bounds::repeatableBounding::WARN:
            {
                WarningInFunction
                    << "value (" << x << ") less than lower "
                    << "bound (" << minLimit << ")" << nl
                    << "    Continuing with the first entry" << endl;

                // Behaviour as per CLAMP
                xDash = minLimit;
                return true;
                break;
            }
            case bounds::repeatableBounding::CLAMP:
            {
                xDash = minLimit;
                return true;
                break;
            }
            case bounds::repeatableBounding::REPEAT:
            {
                // Adjust x to >= minX
                const scalar span = maxLimit - minLimit;
                xDash = fmod(x - minLimit, span) + minLimit;
                break;
            }
        }
    }
    else
    {
        xDash = x;
    }

    return false;
}


template<class Type>
bool Foam::Function1Types::TableBase<Type>::checkMaxBounds
(
    const scalar x,
    scalar& xDash
) const
{
    const scalar minLimit = table_.first().first();
    const scalar maxLimit = table_.last().first();

    if (x > maxLimit)
    {
        switch (bounding_)
        {
            case bounds::repeatableBounding::ERROR:
            {
                FatalErrorInFunction
                    << "value (" << x << ") greater than upper "
                    << "bound (" << maxLimit << ")" << nl
                    << exit(FatalError);
                break;
            }
            case bounds::repeatableBounding::WARN:
            {
                WarningInFunction
                    << "value (" << x << ") greater than upper "
                    << "bound (" << maxLimit << ")" << nl
                    << "    Continuing with the last entry" << endl;

                // Behaviour as per CLAMP
                xDash = maxLimit;
                return true;
                break;
            }
            case bounds::repeatableBounding::CLAMP:
            {
                xDash = maxLimit;
                return true;
                break;
            }
            case bounds::repeatableBounding::REPEAT:
            {
                // Adjust x to >= minX
                const scalar span = maxLimit - minLimit;
                xDash = fmod(x - minLimit, span) + minLimit;
                break;
            }
        }
    }
    else
    {
        xDash = x;
    }

    return false;
}


template<class Type>
void Foam::Function1Types::TableBase<Type>::userTimeToTime(const Time& t)
{
    for (auto& item : table_)
    {
        item.first() = t.userTimeToTime(item.first());
    }

    tableSamplesPtr_.clear();
    interpolatorPtr_.clear();
}


template<class Type>
Type Foam::Function1Types::TableBase<Type>::value(const scalar x) const
{
    scalar xDash = x;

    if (checkMinBounds(x, xDash))
    {
        return table_.first().second();
    }

    if (checkMaxBounds(xDash, xDash))
    {
        return table_.last().second();
    }

    // Use interpolator
    interpolator().valueWeights(xDash, currentIndices_, currentWeights_);

    Type t(currentWeights_[0]*table_[currentIndices_[0]].second());
    for (label i = 1; i < currentIndices_.size(); i++)
    {
        t += currentWeights_[i]*table_[currentIndices_[i]].second();
    }

    return t;
}


template<class Type>
Type Foam::Function1Types::TableBase<Type>::integrate
(
    const scalar x1,
    const scalar x2
) const
{
    // Use interpolator
    interpolator().integrationWeights(x1, x2, currentIndices_, currentWeights_);

    Type sum(currentWeights_[0]*table_[currentIndices_[0]].second());
    for (label i = 1; i < currentIndices_.size(); i++)
    {
       sum += currentWeights_[i]*table_[currentIndices_[i]].second();
    }

    return sum;
}


template<class Type>
Foam::tmp<Foam::scalarField> Foam::Function1Types::TableBase<Type>::x() const
{
    auto tfld = tmp<scalarField>::New(table_.size(), Zero);
    auto& fld = tfld.ref();

    forAll(table_, i)
    {
        fld[i] = table_[i].first();
    }

    return tfld;
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::Function1Types::TableBase<Type>::y() const
{
    auto tfld = tmp<Field<Type>>::New(table_.size(), Zero);
    auto& fld = tfld.ref();

    forAll(table_, i)
    {
        fld[i] = table_[i].second();
    }

    return tfld;
}


template<class Type>
void Foam::Function1Types::TableBase<Type>::writeEntries(Ostream& os) const
{
    if (bounds::repeatableBounding::CLAMP != bounding_)
    {
        os.writeEntry
        (
            "outOfBounds",
            bounds::repeatableBoundingNames[bounding_]
        );
    }

    os.writeEntryIfDifferent<word>
    (
        "interpolationScheme",
        "linear",
        interpolationScheme_
    );
}


template<class Type>
void Foam::Function1Types::TableBase<Type>::writeData(Ostream& os) const
{
    Function1<Type>::writeData(os);

    os  << nl << indent << table_;
    os.endEntry();

    writeEntries(os);
}


// ************************************************************************* //
