/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

#include "lumpedPointState.H"
#include "unitConversion.H"
#include "EulerCoordinateRotation.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::lumpedPointState::inputFormatType
>
Foam::lumpedPointState::formatNames
({
    { inputFormatType::PLAIN, "plain" },
    { inputFormatType::DICTIONARY, "dictionary" },
});


Foam::scalar Foam::lumpedPointState::visLength = 0.1;

bool Foam::lumpedPointState::visUnused = true;


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//! \cond fileScope
static Foam::string getLineNoComment
(
    Foam::ISstream& is,
    const char comment = '#'
)
{
    Foam::string line;
    do
    {
        is.getLine(line);
    }
    while ((line.empty() || line[0] == comment) && is.good());

    return line;
}
//! \endcond


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::lumpedPointState::calcRotations() const
{
    rotationPtr_.reset(new tensorField(angles_.size()));

    auto rotIter = rotationPtr_->begin();

    for (const vector& angles : angles_)
    {
        *rotIter =
            coordinateRotations::euler::rotation(order_, angles, degrees_);

        ++rotIter;
    }
}


void Foam::lumpedPointState::readDict
(
    const dictionary& dict,
    const quaternion::eulerOrder rotOrder,
    const bool degrees
)
{
    dict.readEntry("points", points_);
    dict.readEntry("angles", angles_);
    order_ =
        quaternion::eulerOrderNames.getOrDefault
        (
            "rotationOrder",
            dict,
            rotOrder
        );

    degrees_ = dict.getOrDefault("degrees", degrees);

    rotationPtr_.reset(nullptr);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lumpedPointState::lumpedPointState()
:
    points_(),
    angles_(),
    order_(quaternion::eulerOrder::ZXZ),
    degrees_(false),
    rotationPtr_(nullptr)
{}


Foam::lumpedPointState::lumpedPointState(const lumpedPointState& rhs)
:
    points_(rhs.points_),
    angles_(rhs.angles_),
    order_(rhs.order_),
    degrees_(rhs.degrees_),
    rotationPtr_(nullptr)
{}


Foam::lumpedPointState::lumpedPointState
(
    const pointField& pts,
    const vectorField& ang,
    const quaternion::eulerOrder rotOrder,
    const bool degrees
)
:
    points_(pts),
    angles_(ang),
    order_(rotOrder),
    degrees_(degrees),
    rotationPtr_(nullptr)
{
    if (points_.size() != angles_.size())
    {
        #ifdef FULLDEBUG
        FatalErrorInFunction
            << "Have " << points_.size() << " points but "
            << angles_.size() << " angles" << nl
            << exit(FatalError);
        #else
        WarningInFunction
            << "Have " << points_.size() << " points but "
            << angles_.size() << " angles, resizing angles to match" << nl;
        #endif

        angles_.resize(points_.size(), Zero);
    }
}


Foam::lumpedPointState::lumpedPointState
(
    const pointField& pts,
    const quaternion::eulerOrder rotOrder,
    const bool degrees
)
:
    points_(pts),
    angles_(points_.size(), Zero),
    order_(rotOrder),
    degrees_(degrees),
    rotationPtr_(nullptr)
{}


Foam::lumpedPointState::lumpedPointState
(
    tmp<pointField>& pts,
    const quaternion::eulerOrder rotOrder,
    const bool degrees
)
:
    points_(pts),
    angles_(points_.size(), Zero),
    order_(rotOrder),
    degrees_(false),
    rotationPtr_(nullptr)
{}


Foam::lumpedPointState::lumpedPointState
(
    const dictionary& dict,
    const quaternion::eulerOrder rotOrder,
    const bool degrees
)
:
    points_(),
    angles_(),
    order_(rotOrder),
    degrees_(degrees),
    rotationPtr_(nullptr)
{
    readDict(dict, rotOrder, degrees);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::lumpedPointState::operator=(const lumpedPointState& rhs)
{
    points_  = rhs.points_;
    angles_  = rhs.angles_;
    order_   = rhs.order_;
    degrees_ = rhs.degrees_;

    rotationPtr_.reset(nullptr);
}


void Foam::lumpedPointState::operator+=(const point& origin)
{
    for (point& p : points_)
    {
        p += origin;
    }
}


void Foam::lumpedPointState::scalePoints(const scalar scaleFactor)
{
    if (scaleFactor > 0)
    {
        points_ *= scaleFactor;
    }
}


void Foam::lumpedPointState::relax
(
    const scalar alpha,
    const lumpedPointState& prev
)
{
    points_ = prev.points_ + alpha*(points_ - prev.points_);

    scalar convert = 1.0;
    if (degrees_ != prev.degrees_)
    {
        if (prev.degrees_)
        {
            // Was degrees, now radians
            convert = degToRad();
        }
        else
        {
            // Was radians, now degrees
            convert = radToDeg();
        }
    }

    angles_ = convert*prev.angles_ + alpha*(angles_ - convert*prev.angles_);

    rotationPtr_.reset(nullptr);
}


bool Foam::lumpedPointState::readPlain
(
    Istream& is,
    const quaternion::eulerOrder rotOrder,
    const bool degrees
)
{
    // Assume generic input stream so we can do line-based input
    ISstream& iss = dynamic_cast<ISstream&>(is);

    string line = getLineNoComment(iss);

    label count;
    {
        IStringStream isstr(line);
        isstr >> count;
    }

    points_.resize(count);
    angles_.resize(count);

    count = 0;
    forAll(points_, i)
    {
        iss.getLine(line);
        IStringStream isstr(line);

        isstr
            >> points_[count].x() >> points_[count].y() >> points_[count].z()
            >> angles_[count].x() >> angles_[count].y() >> angles_[count].z();

        ++count;
    }

    points_.resize(count);
    angles_.resize(count);
    order_ = quaternion::eulerOrder::ZXZ;
    degrees_ = false;

    rotationPtr_.reset(nullptr);

    return count;
}


bool Foam::lumpedPointState::readData
(
    Istream& is,
    const quaternion::eulerOrder rotOrder,
    const bool degrees
)
{
    dictionary dict(is);
    readDict(dict, rotOrder, degrees);

    return points_.size();
}


bool Foam::lumpedPointState::writeData(Ostream& os) const
{
    writeDict(os);
    return true;
}


void Foam::lumpedPointState::writeDict(Ostream& os) const
{
    os.writeEntry("points", points_);
    os.writeEntry("angles", angles_);
    if (order_ != quaternion::eulerOrder::ZXZ)
    {
        os.writeEntry("order", quaternion::eulerOrderNames[order_]);
    }
    if (degrees_)
    {
        os.writeEntry("degrees", "true");
    }
}


void Foam::lumpedPointState::writePlain(Ostream& os) const
{
    os  <<"# input for OpenFOAM\n"
        <<"# N, points, angles\n"
        << points_.size() << "\n";

    forAll(points_, i)
    {
        const vector& p = points_[i];

        os  << p.x() << ' ' << p.y() << ' ' << p.z();

        if (i < angles_.size())
        {
            os  << ' ' << angles_[i].x()
                << ' ' << angles_[i].y()
                << ' ' << angles_[i].z() << '\n';
        }
        else
        {
            os  << " 0 0 0\n";
        }
    }
}


bool Foam::lumpedPointState::readData
(
    const inputFormatType fmt,
    const fileName& file,
    const quaternion::eulerOrder rotOrder,
    const bool degrees
)
{
    bool ok = false;
    if (Pstream::master())
    {
        IFstream is(file);

        if (fmt == inputFormatType::PLAIN)
        {
            ok = this->readPlain(is, rotOrder, degrees);
        }
        else
        {
            ok = this->readData(is, rotOrder, degrees);
        }
    }

    if (Pstream::parRun())
    {
        // Broadcast master data to everyone

        Pstream::broadcasts
        (
            UPstream::worldComm,
            ok,
            degrees_,
            points_,
            angles_
        );
    }
    rotationPtr_.reset(nullptr);

    return ok;
}


// ************************************************************************* //
