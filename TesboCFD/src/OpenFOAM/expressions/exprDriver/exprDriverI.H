/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2010-2018 Bernhard Gschaider
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::expressions::exprDriver::setArgument(scalar val) noexcept
{
    arg1Value_ = val;
}


inline Foam::scalar Foam::expressions::exprDriver::argValue() const noexcept
{
    return arg1Value_;
}


inline bool Foam::expressions::exprDriver::hasVariable
(
    const word& name
) const
{
    return variables_.found(name);
}


inline const Foam::expressions::exprResult&
Foam::expressions::exprDriver::variable
(
    const word& name
) const
{
    return variables_[name];
}


inline Foam::expressions::exprResult&
Foam::expressions::exprDriver::variable
(
    const word& name
)
{
    return variables_[name];
}


template<class Type>
inline void Foam::expressions::exprDriver::addUniformVariable
(
    const word& name,
    const Type& val
)
{
    exprResult result;
    result.setSingleValue(val);
    variables_.set(name, std::move(result));
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::expressions::exprDriver::evaluate
(
    const expressions::exprString& expr,
    bool wantPointData
)
{
    parse(expr);

    return getResult<Type>(wantPointData);
}



template<class Type>
inline Type Foam::expressions::exprDriver::evaluateUniform
(
    const expressions::exprString& expr,
    bool wantPointData
)
{
    parse(expr);

    // noWarn = true
    return result_.getUniform(1, true).cref<Type>().first();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline int Foam::expressions::exprDriver::stashedTokenId() const noexcept
{
    return stashedTokenId_;
}


inline int Foam::expressions::exprDriver::resetStashedTokenId
(
    int tokenId
) const noexcept
{
    int old = stashedTokenId_;
    stashedTokenId_ = tokenId;
    return old;
}


inline bool Foam::expressions::exprDriver::debugScanner() const noexcept
{
    return debugScanner_;
}


inline bool Foam::expressions::exprDriver::debugParser() const noexcept
{
    return debugParser_;
}


inline bool Foam::expressions::exprDriver::searchRegistry() const noexcept
{
    return (searchCtrl_ & searchControls::SEARCH_REGISTRY);
}


inline bool Foam::expressions::exprDriver::searchFiles() const noexcept
{
    return (searchCtrl_ & searchControls::SEARCH_FILES);
}


inline bool Foam::expressions::exprDriver::cacheReadFields() const noexcept
{
    return (searchCtrl_ & searchControls::CACHE_READ_FIELDS);
}


// ************************************************************************* //
