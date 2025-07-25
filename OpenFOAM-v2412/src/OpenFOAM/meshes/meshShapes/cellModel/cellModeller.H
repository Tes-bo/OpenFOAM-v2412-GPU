/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Class
    Foam::cellModeller

Description
    Compatibility definitions of static cellModel lookups.

    Superseded (NOV-2017) by cellModel methods.

See Also
    Foam::cellModel

\*---------------------------------------------------------------------------*/

#ifndef cellModeller_H
#define cellModeller_H

#include "cellModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

struct cellModeller
{
    //- Deprecated(2017-11) equivalent to cellModel::ptr static method.
    //  \deprecated(2017-11) use cellModel::ptr instead
    FOAM_DEPRECATED_FOR(2017-11, "cellModel::ptr() static method")
    static const cellModel* lookup(const word& modelName)
    {
        return cellModel::ptr(modelName);
    }

    //- Deprecated(2017-11) equivalent to cellModel::ptr static method.
    //  \deprecated(2017-11) use cellModel::ptr instead
    FOAM_DEPRECATED_FOR(2017-11, "cellModel::ptr() static method")
    static const cellModel* lookup(const label modelIndex)
    {
        return cellModel::ptr(modelIndex);
    }
};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
