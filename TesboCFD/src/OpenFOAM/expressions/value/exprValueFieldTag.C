/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023-2024 OpenCFD Ltd.
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

#include "exprValueFieldTag.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::expressions::exprValueFieldTag::empty() const noexcept
{
    return
    (
        uniformity_ == Foam::Detail::ListPolicy::uniformity::EMPTY
    );
}


bool Foam::expressions::exprValueFieldTag::is_uniform() const noexcept
{
    return
    (
        uniformity_ == Foam::Detail::ListPolicy::uniformity::UNIFORM
    );
}


bool Foam::expressions::exprValueFieldTag::is_nonuniform() const noexcept
{
    return
    (
        uniformity_ == Foam::Detail::ListPolicy::uniformity::NONUNIFORM
        // Extra safety for direct reductions?
        // || uniformity_ == Foam::Detail::ListPolicy::uniformity::MIXED
    );
}


const Foam::expressions::exprValue&
Foam::expressions::exprValueFieldTag::value() const noexcept
{
    return value_;
}


void Foam::expressions::exprValueFieldTag::set_empty()
{
    uniformity_ = Foam::Detail::ListPolicy::uniformity::EMPTY;
    value_ = Foam::zero{};
}


void Foam::expressions::exprValueFieldTag::set_nouniform()
{
    uniformity_ = Foam::Detail::ListPolicy::uniformity::NONUNIFORM;
    value_ = Foam::zero{};
}


int Foam::expressions::exprValueFieldTag::compare
(
    const exprValueFieldTag& rhs
) const
{
    if (uniformity_ != rhs.uniformity_)
    {
        // First compare by uniformity
        return (int(uniformity_) - int(rhs.uniformity_));
    }
    if (this == &rhs)
    {
        // Identical objects
        return 0;
    }

    return value_.compare(rhs.value_);
}


bool Foam::expressions::exprValueFieldTag::equal
(
    const exprValueFieldTag& rhs
) const
{
    return (value_ == rhs.value_);
}


void Foam::expressions::exprValueFieldTag::reduce()
{
    if (!UPstream::is_parallel())
    {
        // Nothing to do
        return;
    }

    // Two-stage reduction
    // ~~~~~~~~~~~~~~~~~~~
    //
    // Fields will usually be non-uniform somewhere, so first check with
    // the cheapest option (bit-wise Allreduce).
    // Only if they are uniform (with/without empty) do we actually
    // need to compare values.

    typedef unsigned char bitmask_type;

    bitmask_type shape = static_cast<bitmask_type>(uniformity_);

    // Step 1
    // ~~~~~~
    Foam::reduce
    (
        shape,
        bitOrOp<bitmask_type>{},
        UPstream::msgType(),  // ignored
        UPstream::worldComm
    );

    // Step 2
    // ~~~~~~
    if
    (
        shape == static_cast<bitmask_type>
        (
            Foam::Detail::ListPolicy::uniformity::EMPTY
        )
    )
    {
        // no-op (empty everywhere)
        value_ = Foam::zero{};
    }
    else if
    (
        shape == static_cast<bitmask_type>
        (
            Foam::Detail::ListPolicy::uniformity::UNIFORM
        )
    )
    {
        // Ranks are locally uniform (or empty), need to check values too
        Foam::reduce
        (
            *this,
            exprValueFieldTag::combineOp{},
            UPstream::msgType(),
            UPstream::worldComm
        );
    }
    else
    {
        // Field is global non-empty and not uniform
        set_nouniform();
    }
}


Foam::expressions::exprValueFieldTag
Foam::expressions::exprValueFieldTag::returnReduce
(
    const exprValueFieldTag& tag
)
{
    exprValueFieldTag work(tag);
    work.reduce();
    return work;
}


void Foam::expressions::exprValueFieldTag::combine
(
    const exprValueFieldTag& b
)
{
    exprValueFieldTag& a = *this;

    if (b.empty())
    {
        // no-op
        return;
    }
    else if (a.empty())
    {
        a = b;
    }
    else if (a.is_uniform() && (!b.is_uniform() || !a.equal(b)))
    {
        // Handle two cases:
        // 1. uniform / non-uniform
        // 2. uniform / uniform, but with different values
        a.set_nouniform();
    }

    // No meaningful value if it is not uniform.
    // So use zero, but keep the type
    if (!a.is_uniform())
    {
        a.value_ = Foam::zero{};
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

bool Foam::expressions::exprValueFieldTag::
operator==(const exprValueFieldTag& rhs) const
{
    if (uniformity_ != rhs.uniformity_)
    {
        // Uniformity must match
        return false;
    }
    else if (this == &rhs)
    {
        return true;
    }

    return (value_ == rhs.value_);
}


bool Foam::expressions::exprValueFieldTag::
operator<(const exprValueFieldTag& rhs) const
{
    return (this->compare(rhs) < 0);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

void Foam::expressions::exprValueFieldTag::read(Istream& is)
{
    label uniformTag;

    is.readBegin("fieldTag");

    is >> uniformTag;
    uniformity_ = int(uniformTag);
    value_.read(is);

    is.readEnd("fieldTag");
}


void Foam::expressions::exprValueFieldTag::write(Ostream& os) const
{
    os  << token::BEGIN_LIST
        << label(uniformity_) << token::SPACE;
    value_.write(os, false);  // No pruning
    os  << token::END_LIST;
}


void Foam::expressions::exprValueFieldTag::print(Ostream& os) const
{
    os  << "{ uniform:"
        << label(uniformity_)
        << " type:" << label(value_.typeCode())
        << " value: " << value_ << " }";
}


Foam::Istream& Foam::operator>>
(
    Istream& is,
    expressions::exprValueFieldTag& tag
)
{
    tag.read(is);
    return is;
}


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const expressions::exprValueFieldTag& tag
)
{
    tag.write(os);
    return os;
}


// ************************************************************************* //
