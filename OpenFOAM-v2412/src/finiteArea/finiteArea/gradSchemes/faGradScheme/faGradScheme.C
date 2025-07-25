/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

#include "fa.H"
#include "objectRegistry.H"
#include "solution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class Type>
tmp<gradScheme<Type>> gradScheme<Type>::New
(
    const faMesh& mesh,
    Istream& schemeData
)
{
    if (fa::debug)
    {
        InfoInFunction
            << "constructing gradScheme<Type>"
            << endl;
    }

    if (schemeData.eof())
    {
        FatalIOErrorInFunction(schemeData)
            << "Grad scheme not specified" << nl << nl
            << "Valid grad schemes are :" << nl
            << IstreamConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    const word schemeName(schemeData);

    auto* ctorPtr = IstreamConstructorTable(schemeName);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            schemeData,
            "grad",
            schemeName,
            *IstreamConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return ctorPtr(mesh, schemeData);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector, Type>::type,
        faPatchField,
        areaMesh
    >
>
gradScheme<Type>::grad
(
    const GeometricField<Type, faPatchField, areaMesh>& vsf,
    const word& name
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;
    typedef GeometricField<GradType, faPatchField, areaMesh> GradFieldType;

    GradFieldType* pgGrad =
        mesh().thisDb().template getObjectPtr<GradFieldType>(name);

    if (!this->mesh().cache(name)) // || this->mesh().changing()
    {
        // Delete any old occurrences to avoid double registration
        if (pgGrad && pgGrad->ownedByRegistry())
        {
            solution::cachePrintMessage("Deleting", name, vsf);
            delete pgGrad;
        }

        solution::cachePrintMessage("Calculating", name, vsf);
        return calcGrad(vsf, name);
    }


    if (!pgGrad)
    {
        solution::cachePrintMessage("Calculating and caching", name, vsf);

        pgGrad = calcGrad(vsf, name).ptr();
        regIOobject::store(pgGrad);
    }
    else
    {
        if (pgGrad->upToDate(vsf))
        {
            solution::cachePrintMessage("Reusing", name, vsf);
        }
        else
        {
            solution::cachePrintMessage("Updating", name, vsf);
            delete pgGrad;

            pgGrad = calcGrad(vsf, name).ptr();
            regIOobject::store(pgGrad);
        }
    }

    return *pgGrad;
}


template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector, Type>::type,
        faPatchField,
        areaMesh
    >
>
gradScheme<Type>::grad
(
    const GeometricField<Type, faPatchField, areaMesh>& vsf
) const
{
    return grad(vsf, "grad(" + vsf.name() + ')');
}


template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector, Type>::type,
        faPatchField,
        areaMesh
    >
>
gradScheme<Type>::grad
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tvsf
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;
    typedef GeometricField<GradType, faPatchField, areaMesh> GradFieldType;

    tmp<GradFieldType> tgrad = grad(tvsf());
    tvsf.clear();
    return tgrad;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
