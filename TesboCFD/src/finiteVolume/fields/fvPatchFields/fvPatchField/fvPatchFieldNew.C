/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::fvPatchField<Type>> Foam::fvPatchField<Type>::New
(
    const word& patchFieldType,
    const word& actualPatchType,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
{
    DebugInFunction
        << "patchFieldType = " << patchFieldType
        << " [" << actualPatchType
        << "] : " << p.type() << " name = " << p.name() << endl;

    auto* ctorPtr = patchConstructorTable(patchFieldType);

    if (!ctorPtr)
    {
        FatalErrorInLookup
        (
            "patchField",
            patchFieldType,
            *patchConstructorTablePtr_
        ) << exit(FatalError);
    }

    auto* patchTypeCtor = patchConstructorTable(p.type());

    if (actualPatchType.empty() || actualPatchType != p.type())
    {
        if (patchTypeCtor)
        {
            return patchTypeCtor(p, iF);
        }
        else
        {
            return ctorPtr(p, iF);
        }
    }


    tmp<fvPatchField<Type>> tpfld(ctorPtr(p, iF));

    // If constraint type: override and store patchType
    if (patchTypeCtor)
    {
        tpfld.ref().patchType() = actualPatchType;
    }
    return tpfld;
}


template<class Type>
Foam::tmp<Foam::fvPatchField<Type>> Foam::fvPatchField<Type>::New
(
    const word& patchFieldType,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
{
    return New(patchFieldType, word::null, p, iF);
}


template<class Type>
Foam::tmp<Foam::fvPatchField<Type>> Foam::fvPatchField<Type>::New
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
{
    const word patchFieldType(dict.get<word>("type"));

    word actualPatchType;
    dict.readIfPresent("patchType", actualPatchType, keyType::LITERAL);

    DebugInFunction
        << "patchFieldType = " << patchFieldType
        << " [" << actualPatchType
        << "] : " << p.type() << " name = " << p.name() << endl;

    auto* ctorPtr = dictionaryConstructorTable(patchFieldType);

    if (!ctorPtr)
    {
        if (!fvPatchFieldBase::disallowGenericPatchField)
        {
            ctorPtr = dictionaryConstructorTable("generic");
        }

        if (!ctorPtr)
        {
            FatalIOErrorInFunction(dict)
                << "Unknown patchField type " << patchFieldType
                << " for patch type " << p.type() << nl << nl
                << "Valid patchField types :" << endl
                << dictionaryConstructorTablePtr_->sortedToc()
                << exit(FatalIOError);
        }
    }

    if (actualPatchType.empty() || actualPatchType != p.type())
    {
        auto* patchTypeCtor = dictionaryConstructorTable(p.type());

        if (patchTypeCtor && patchTypeCtor != ctorPtr)
        {
            FatalIOErrorInFunction(dict)
                << "inconsistent patch and patchField types for\n"
                   "    patch type " << p.type()
                << " and patchField type " << patchFieldType
                << exit(FatalIOError);
        }
    }

    return ctorPtr(p, iF, dict);
}


template<class Type>
Foam::tmp<Foam::fvPatchField<Type>> Foam::fvPatchField<Type>::New
(
    const fvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& pfMapper
)
{
    DebugInFunction
        << "patchFieldType = " << ptf.type()
        << " : " << p.type() << " name = " << p.name() << endl;

    auto* ctorPtr = patchMapperConstructorTable(ptf.type());

    if (!ctorPtr)
    {
        FatalErrorInLookup
        (
            "patchField",
            ptf.type(),
            *patchMapperConstructorTablePtr_
        ) << exit(FatalError);
    }

    return ctorPtr(ptf, p, iF, pfMapper);
}


// ************************************************************************* //
