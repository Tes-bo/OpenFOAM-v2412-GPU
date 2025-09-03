/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::expressions::fieldExpr::parseDriver::getField
(
    const word& name
) const
{
    bool hasPointData = false;

    refPtr<expressions::exprResult> tvar;

    if (hasVariable(name) && variable(name).isType<Type>())
    {
        tvar.cref(variable(name));
        hasPointData = tvar().isPointData();
    }


    if (tvar)
    {
        const auto& var = tvar.cref();
        const Field<Type>& vals = var.cref<Type>();

        const label len = (hasPointData ? this->pointSize() : this->size());

        if (returnReduceAnd(vals.size() == len))
        {
            // Return a copy of the field
            return tmp<Field<Type>>::New(vals);
        }

        if (!var.isUniform())
        {
            WarningInFunction
                << "Variable " << name
                << " is nonuniform and does not fit the size"
                << ". Using average" << endl;
        }

        return tmp<Field<Type>>::New(this->size(), gAverage(vals));
    }

    return nullptr;
}


// ************************************************************************* //
