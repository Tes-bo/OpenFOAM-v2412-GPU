/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020-2022 OpenCFD Ltd.
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

#include "uniformInterpolationTable.H"
#include "Time.H"
#include "IOstream.H"
#include "IOdictionary.H"

// * * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * //

template<class Type>
void Foam::uniformInterpolationTable<Type>::checkTable() const
{
    if (size() < 2)
    {
        FatalErrorInFunction
            << "Table " << name() << ": must have at least 2 values." << nl
            << "Table size = " << size() << nl
            << "    min, interval width = " << x0_ << ", " << dx_ << nl
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::uniformInterpolationTable<Type>::uniformInterpolationTable
(
    const IOobject& io,
    bool readFields
)
:
    IOobject(io),
    List<scalar>(2, Zero),
    x0_(0.0),
    dx_(1.0),
    log10_(false),
    bound_(false)
{
    if (readFields)
    {
        IOdictionary dict(io);

        dict.readEntry("data", *this);
        dict.readEntry("x0", x0_);
        dict.readEntry("dx", dx_);
        dict.readIfPresent("log10", log10_);
        dict.readIfPresent("bound", bound_);
    }

    checkTable();
}


template<class Type>
Foam::uniformInterpolationTable<Type>::uniformInterpolationTable
(
    const word& tableName,
    const objectRegistry& db,
    const dictionary& dict,
    const bool initialiseOnly
)
:
    IOobject
    (
        tableName,
        db.time().constant(),
        db,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
        // No register: could be used by multiple patches if used in BCs
    ),
    List<scalar>(2, Zero),
    x0_(dict.get<scalar>("x0")),
    dx_(dict.get<scalar>("dx")),
    log10_(dict.getOrDefault<Switch>("log10", false)),
    bound_(dict.getOrDefault<Switch>("bound", false))
{
    if (initialiseOnly)
    {
        const scalar xMax = dict.get<scalar>("xMax");
        const label nIntervals = static_cast<label>(xMax - x0_)/dx_ + 1;
        this->setSize(nIntervals);
    }
    else
    {
        dict.readEntry("data", *this);
    }

    checkTable();
}


template<class Type>
Foam::uniformInterpolationTable<Type>::uniformInterpolationTable
(
    const uniformInterpolationTable& uit
)
:
    IOobject(uit),
    List<scalar>(uit),
    x0_(uit.x0_),
    dx_(uit.dx_),
    log10_(uit.log10_),
    bound_(uit.bound_)
{
    checkTable();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::uniformInterpolationTable<Type>::~uniformInterpolationTable()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
Type Foam::uniformInterpolationTable<Type>::interpolate(scalar x) const
{
    if (bound_)
    {
        x = clamp(x, x0_, (xMax() - SMALL*dx_));
    }
    else
    {
        if (x < x0_)
        {
            FatalErrorInFunction
                << "Supplied value is less than minimum table value:" << nl
                << "xMin=" << x0_ << ", xMax=" << xMax() << ", x=" << x << nl
                << exit(FatalError);
        }

        if (x > xMax())
        {
            FatalErrorInFunction
                << "Supplied value is greater than maximum table value:" << nl
                << "xMin=" << x0_ << ", xMax=" << xMax() << ", x=" << x << nl
                << exit(FatalError);
        }
    }

    const label i = static_cast<label>((x - x0_)/dx_);

    const scalar xLo = x0_ + i*dx_;

    Type fx = (x - xLo)/dx_*(operator[](i+1) - operator[](i)) + operator[](i);

    if (debug)
    {
        Info<< "Table: " << name() << ", x=" << x
            << ", x_lo=" << xLo << ", x_hi=" << xLo + dx_
            << ", f(x_lo)=" << operator[](i) << ", f(x_hi)=" << operator[](i+1)
            << ", f(x)=" << fx << endl;
    }

    return fx;
}


template<class Type>
Type Foam::uniformInterpolationTable<Type>::interpolateLog10
(
    scalar x
) const
{
    if (log10_)
    {
        if (x > 0)
        {
            x = ::log10(x);
        }
        else if (bound_ && (x <= 0))
        {
            x = x0_;
        }
        else
        {
            FatalErrorInFunction
                << "Table " << name() << nl
                << "Supplied value must be greater than 0 when in log10 mode"
                << nl << "x=" << x << nl << exit(FatalError);
        }
    }

    return interpolate(x);
}


template<class Type>
void Foam::uniformInterpolationTable<Type>::write() const
{
    IOdictionary dict(*this);

    dict.add("data", static_cast<const List<scalar>&>(*this));
    dict.add("x0", x0_);
    dict.add("dx", dx_);
    if (log10_)
    {
        dict.add("log10", log10_);
    }
    if (bound_)
    {
        dict.add("bound", bound_);
    }

    dict.regIOobject::writeObject
    (
        IOstreamOption(IOstreamOption::ASCII, dict.time().writeCompression()),
        true
    );
}


// ************************************************************************* //
