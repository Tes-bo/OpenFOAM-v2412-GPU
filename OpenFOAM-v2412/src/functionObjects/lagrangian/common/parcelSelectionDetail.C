/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2019 OpenCFD Ltd.
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

#include "cloud.H"
#include "parcelSelectionDetail.H"
#include "scalarPredicates.H"
#include "labelField.H"
#include "scalarField.H"
#include "pointField.H"
#include "ListListOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::Detail::parcelSelection::actionType
>
Foam::Detail::parcelSelection::actionNames
({
    { actionType::ALL, "all" },
    { actionType::CLEAR, "clear" },
    { actionType::INVERT, "invert" },
    { actionType::USE, "use" },
    { actionType::ADD, "add" },
    { actionType::SUBTRACT, "subtract" },
    { actionType::SUBSET, "subset" },
    { actionType::IGNORE, "ignore" },
});


const Foam::Enum
<
    Foam::Detail::parcelSelection::sourceType
>
Foam::Detail::parcelSelection::sourceNames
({
    { sourceType::FIELD, "field" },
    { sourceType::STRIDE, "stride" },
});


const Foam::Enum
<
    Foam::Detail::parcelSelection::logicType
> Foam::Detail::parcelSelection::logicNames
({
    { logicType::AND, "and" },
    { logicType::OR, "or" },
});


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
    template<class Type, class Predicate, class AccessOp>
    static void apply
    (
        bitSet& selection,
        const Detail::parcelSelection::actionType action,
        const Predicate& accept,
        const UList<Type>& list,
        const AccessOp& aop
    )
    {
        using actionType = Detail::parcelSelection::actionType;

        const label len = selection.size();

        switch (action)
        {
            case actionType::ADD:
            case actionType::USE:
            {
                if (actionType::USE == action)
                {
                    // USE = CLEAR + ADD  (ie, only use this selection)
                    selection = false;
                }

                for (label parceli = 0; parceli < len; ++parceli)
                {
                    if (accept(aop(list[parceli])))
                    {
                        selection.set(parceli);
                    }
                }
            }
            break;

            case actionType::SUBTRACT:
            {
                for (label parceli = 0; parceli < len; ++parceli)
                {
                    if (accept(aop(list[parceli])))
                    {
                        selection.unset(parceli);
                    }
                }
            }
            break;

            case actionType::SUBSET:
            {
                for (const label parceli : selection)
                {
                    if (!accept(aop(list[parceli])))
                    {
                        selection.unset(parceli);
                    }
                }
            }
            break;

            default:
                break;
        }
    }

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::Detail::parcelSelection::calculateFilter
(
    const objectRegistry& obrTmp,
    const bool log
)
{
    if (parcelSelect_.empty())
    {
        parcelAddr_.clear();
        return false;
    }

    // Start with all parcels unselected

    // Number of parcels (locally)
    const auto* pointsPtr = cloud::findIOPosition(obrTmp);
    label nParcels = pointsPtr->size();

    parcelAddr_.reset();
    parcelAddr_.resize(nParcels);

    reduce(nParcels, sumOp<label>());

    Log << "Applying parcel filtering to " << nParcels << " parcels" << nl;

    if (!nParcels)
    {
        parcelAddr_.clear();
        return false;
    }

    // The unary function type(s) for testing a scalar.
    // Allocate 3 slots.
    // 0 is the test
    // 1,2 are for storage of composite tests (eg, and/or logic)
    predicates::scalars tests(3);

    for (const entry& dEntry : parcelSelect_)
    {
        if (!dEntry.isDict())
        {
            WarningInFunction
                << "Ignoring non-dictionary entry "
                << dEntry << endl;
            continue;
        }

        const dictionary& dict = dEntry.dict();

        // A very limited number of sources (stride, field)
        // and actions (all add subtract subset) so handle manually

        auto action = actionNames.get("action", dict);

        // These ones we do directly
        switch (action)
        {
            case actionType::ALL:
                Log << "- select all" << nl;
                parcelAddr_ = true;
                continue;
                break;

            case actionType::CLEAR:
                Log << "- clear" << nl;
                parcelAddr_ = false;
                continue;
                break;

            case actionType::INVERT:
                Log << "- invert" << nl;
                parcelAddr_.flip();
                continue;
                break;

            case actionType::IGNORE:
                continue;
                break;

            default:
                break;
        }

        // The others need a source
        // Need a source
        const auto source = sourceNames.get("source", dict);

        switch (source)
        {
            case sourceType::STRIDE:
            {
                const label stride = dict.get<label>("stride");

                const labelField& ids =
                    obrTmp.lookupObject<labelField>("origId");

                Log << "- " << actionNames[action]
                    << " stride " << stride << nl;

                if (stride <= 0)
                {
                    WarningInFunction
                        << nl
                        << "Ignoring bad value for stride=" << stride << nl
                        << endl;
                }
                else if (stride == 1)
                {
                    // More efficient handling of stride 1, but should
                    // not really be using stride == 1.
                    switch (action)
                    {
                        case actionType::ADD:
                            parcelAddr_ = true;
                            break;

                        case actionType::SUBTRACT:
                            parcelAddr_ = false;
                            break;

                        default:
                            break;
                    }
                }
                else
                {
                    // Using stride > 1
                    apply
                    (
                        parcelAddr_,
                        action,
                        [=](const label id) -> bool { return !(id % stride); },
                        ids,
                        accessOp<label>()  // pass-through
                    );
                }
            }
            break;

            case sourceType::FIELD:
            {
                const word fieldName(dict.get<word>("field"));

                const auto* labelFld =
                    obrTmp.findObject<labelField>(fieldName);

                const auto* scalarFld =
                    obrTmp.findObject<scalarField>(fieldName);

                const auto* vectorFld =
                    obrTmp.findObject<vectorField>(fieldName);


                Log << "- " << actionNames[action] << " field " << fieldName;

                if (!labelFld && !scalarFld && !vectorFld)
                {
                    WarningInFunction
                        << nl
                        << "No scalar/vector parcel field: " << fieldName
                        << " ignoring selection" << nl
                        << endl;
                    continue;
                }

                const entry& e = dict.lookupEntry("accept", keyType::LITERAL);

                ITstream& is = e.stream();

                if (4 == is.size())
                {
                    // 4 tokens:
                    //  -> (op val)

                    Tuple2<word,scalar> expr1(is);
                    e.checkITstream(is);

                    tests.first() = predicates::scalars::operation(expr1);

                    Log << " : " << expr1;
                }
                else if (9 == is.size())
                {
                    // 9 tokens:
                    //  ->  (op val) and (op val)
                    //  ->  (op val) or (op val)

                    Tuple2<word,scalar> expr1(is);
                    word logicName(is);
                    Tuple2<word,scalar> expr2(is);
                    e.checkITstream(is);

                    logicType logic = logicNames[logicName];
                    tests[1] = predicates::scalars::operation(expr1);
                    tests[2] = predicates::scalars::operation(expr2);

                    switch (logic)
                    {
                        case logicType::AND:
                            tests.first() =
                                predicates::scalars::andOp(tests[1], tests[2]);
                            break;

                        case logicType::OR:
                            tests.first() =
                                predicates::scalars::orOp(tests[1], tests[2]);
                            break;
                    }

                    Log << " : " << expr1 << ' ' << logicName << ' ' << expr2;
                }
                else
                {
                    action = actionType::IGNORE;

                    // Use the following to always provoke an error.
                    e.checkITstream(is);
                }

                if (labelFld)
                {
                    apply
                    (
                        parcelAddr_,
                        action,
                        tests.first(),
                        *labelFld,
                        accessOp<label>()  // pass-through
                    );
                }
                else if (scalarFld)
                {
                    apply
                    (
                        parcelAddr_,
                        action,
                        tests.first(),
                        *scalarFld,
                        accessOp<scalar>() // pass-through
                    );
                }
                else if (vectorFld)
                {
                    apply
                    (
                        parcelAddr_,
                        action,
                        tests.first(),
                        *vectorFld,
                        [](const vector& val) -> scalar
                        {
                            return Foam::mag(val);
                        }
                    );
                }

                Log << endl;
            }
            break;

            default:
                break;
        }
    }

    return true;
}


// ************************************************************************* //
