/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2018 Bernhard Gschaider
    Copyright (C) 2019 OpenCFD Ltd.
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

template<class T>
bool Foam::expressions::exprResultStack::pushChecked
(
    const exprResult& result
)
{
    if (!isType<T>())
    {
        return false;
    }

    // The value to push
    T val(Zero);

    const Field<T>& resultField = result.cref<T>();

    if (!resultField.empty())
    {
        val = resultField.front();
    }

    this->ref<T>().push_back(val);

    return true;
}


template<class T>
bool Foam::expressions::exprResultStack::popChecked
(
    exprResult& result
)
{
    if (!isType<T>())
    {
        return false;
    }

    // The popped value
    T val(Zero);

    Field<T>& oldField = this->ref<T>();

    if (!oldField.empty())
    {
        val = oldField.back();
        oldField.pop_back();
    }

    result.setSingleValue(val);

    return true;
}


// ************************************************************************* //
