/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd.
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

#include "Reaction.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * Static Data * * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::label Foam::Reaction<ReactionThermo>::nUnNamedReactions = 0;


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class ReactionThermo>
void Foam::Reaction<ReactionThermo>::reactionStr
(
    OStringStream& reaction,
    const speciesTable& species,
    const List<specieCoeffs>& reactCoeffs
)
{
    for (label i = 0; i < reactCoeffs.size(); ++i)
    {
        const specieCoeffs& coeff = reactCoeffs[i];

        if (i)
        {
            reaction << " + ";
        }
        if (mag(coeff.stoichCoeff - 1) > SMALL)
        {
            reaction << coeff.stoichCoeff;
        }
        reaction << species[coeff.index];
        if (mag(coeff.exponent - coeff.stoichCoeff) > SMALL)
        {
            reaction << "^" << coeff.exponent;
        }
    }
}


template<class ReactionThermo>
void Foam::Reaction<ReactionThermo>::reactionStrLeft
(
    OStringStream& reaction
) const
{
    reactionStr(reaction, species_, lhs_);
}


template<class ReactionThermo>
void Foam::Reaction<ReactionThermo>::reactionStrRight
(
    OStringStream& reaction
) const
{
    reactionStr(reaction, species_, rhs_);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ReactionThermo>
Foam::label Foam::Reaction<ReactionThermo>::getNewReactionID()
{
    return nUnNamedReactions++;
}


template<class ReactionThermo>
Foam::string Foam::Reaction<ReactionThermo>::reactionStr
(
    OStringStream& reaction
) const
{
    reactionStrLeft(reaction);
    reaction << " = ";
    reactionStrRight(reaction);
    return reaction.str();
}


template<class ReactionThermo>
void Foam::Reaction<ReactionThermo>::setThermo
(
    const ReactionTable<ReactionThermo>& thermoDatabase
)
{

    typename ReactionThermo::thermoType rhsThermo
    (
        rhs_[0].stoichCoeff
        *(*thermoDatabase[species_[rhs_[0].index]]).W()
        *(*thermoDatabase[species_[rhs_[0].index]])
    );

    for (label i=1; i<rhs_.size(); ++i)
    {
        rhsThermo +=
            rhs_[i].stoichCoeff
        *(*thermoDatabase[species_[rhs_[i].index]]).W()
        *(*thermoDatabase[species_[rhs_[i].index]]);
    }

    typename ReactionThermo::thermoType lhsThermo
    (
        lhs_[0].stoichCoeff
       *(*thermoDatabase[species_[lhs_[0].index]]).W()
       *(*thermoDatabase[species_[lhs_[0].index]])
    );

    for (label i=1; i<lhs_.size(); ++i)
    {
        lhsThermo +=
            lhs_[i].stoichCoeff
           *(*thermoDatabase[species_[lhs_[i].index]]).W()
           *(*thermoDatabase[species_[lhs_[i].index]]);
    }

    ReactionThermo::thermoType::operator=(lhsThermo == rhsThermo);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


template<class ReactionThermo>
Foam::Reaction<ReactionThermo>::Reaction
(
    const speciesTable& species,
    const List<specieCoeffs>& lhs,
    const List<specieCoeffs>& rhs,
    const ReactionTable<ReactionThermo>& thermoDatabase,
    bool initReactionThermo
)
:
    ReactionThermo::thermoType(*thermoDatabase[species[0]]),
    name_("un-named-reaction-" + Foam::name(getNewReactionID())),
    species_(species),
    lhs_(lhs),
    rhs_(rhs)
{
    if (initReactionThermo)
    {
        setThermo(thermoDatabase);
    }
}


template<class ReactionThermo>
Foam::Reaction<ReactionThermo>::Reaction
(
    const Reaction<ReactionThermo>& r,
    const speciesTable& species
)
:
    ReactionThermo::thermoType(r),
    name_(r.name() + "Copy"),
    species_(species),
    lhs_(r.lhs_),
    rhs_(r.rhs_)
{}


template<class ReactionThermo>
Foam::Reaction<ReactionThermo>::specieCoeffs::specieCoeffs
(
    const speciesTable& species,
    Istream& is,
    bool failUnknownSpecie
)
{
    token t(is);
    if (t.isNumber())
    {
        stoichCoeff = t.number();
        is >> t;
    }
    else
    {
        stoichCoeff = 1.0;
    }

    exponent = stoichCoeff;

    if (t.isWord())
    {
        word specieName = t.wordToken();

        const size_t i = specieName.find('^');

        if (i != word::npos)
        {
            exponent = atof(specieName.substr(i + 1).c_str());
            specieName.resize(i);
        }

        // Lookup specie name: -1 if not found
        index = species.find(specieName);

        if (failUnknownSpecie && index < 0)
        {
            FatalErrorInFunction
                << "Unknown specie " << specieName << nl
                << "Not in " << flatOutput(species) << endl
                << exit(FatalError);
        }
    }
    else
    {
        FatalIOErrorInFunction(is)
            << "Expected a word but found " << t.info()
            << exit(FatalIOError);
    }
}


template<class ReactionThermo>
void Foam::Reaction<ReactionThermo>::setLRhs
(
    Istream& is,
    const speciesTable& species,
    List<specieCoeffs>& lhs,
    List<specieCoeffs>& rhs,
    bool failUnknownSpecie
)
{
    DynamicList<specieCoeffs> dlrhs;

    bool parsingRight = false;

    while (is.good())
    {
        dlrhs.push_back(specieCoeffs(species, is, failUnknownSpecie));

        if (dlrhs.back().index < 0)
        {
            dlrhs.pop_back();
        }

        if (is.good())
        {
            token t(is);
            if (t.isPunctuation())
            {
                if (t == token::ADD)
                {
                }
                else if (t == token::ASSIGN)
                {
                    if (parsingRight)
                    {
                        FatalErrorInFunction
                            << "Multiple '=' in reaction equation" << endl
                            << exit(FatalError);
                    }

                    lhs = dlrhs;
                    dlrhs.clear();
                    parsingRight = true;
                }
                else
                {
                    FatalErrorInFunction
                        << "Unknown punctuation token '" << t
                        << "' in reaction equation" << endl
                        << exit(FatalError);
                }
            }
            else
            {
                rhs = dlrhs;
                is.putBack(t);
                return;
            }
        }
        else if (parsingRight)
        {
            if (!dlrhs.empty())
            {
                rhs = dlrhs;
            }
            return;
        }
    }

    FatalIOErrorInFunction(is)
        << "Cannot continue reading reaction data from stream"
        << exit(FatalIOError);
}


template<class ReactionThermo>
Foam::Reaction<ReactionThermo>::Reaction
(
    const speciesTable& species,
    const ReactionTable<ReactionThermo>& thermoDatabase,
    const dictionary& dict,
    bool initReactionThermo,
    bool failUnknownSpecie
)
:
    ReactionThermo::thermoType(*thermoDatabase[species[0]]),
    name_(dict.dictName()),
    species_(species)
{
    ICharStream reactionIs(dict.getString("reaction"));

    setLRhs
    (
        reactionIs,
        species_,
        lhs_,
        rhs_,
        failUnknownSpecie
    );

    if (initReactionThermo)
    {
        setThermo(thermoDatabase);
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::autoPtr<Foam::Reaction<ReactionThermo>>
Foam::Reaction<ReactionThermo>::New
(
    const speciesTable& species,
    const ReactionTable<ReactionThermo>& thermoDatabase,
    const dictionary& dict
)
{
    const word reactionTypeName(dict.get<word>("type"));

    auto* ctorPtr = dictionaryConstructorTable(reactionTypeName);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "reaction",
            reactionTypeName,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<Reaction<ReactionThermo>>
    (
        ctorPtr(species, thermoDatabase, dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ReactionThermo>
void Foam::Reaction<ReactionThermo>::write(Ostream& os) const
{
    OStringStream reaction;
    os.writeEntry("reaction", reactionStr(reaction));
}


template<class ReactionThermo>
Foam::scalar Foam::Reaction<ReactionThermo>::kf
(
    const scalar p,
    const scalar T,
    const scalarField& c
) const
{
    return 0.0;
}


template<class ReactionThermo>
Foam::scalar Foam::Reaction<ReactionThermo>::kr
(
    const scalar kfwd,
    const scalar p,
    const scalar T,
    const scalarField& c
) const
{
    return 0.0;
}


template<class ReactionThermo>
Foam::scalar Foam::Reaction<ReactionThermo>::kr
(
    const scalar p,
    const scalar T,
    const scalarField& c
) const
{
    return 0.0;
}


template<class ReactionThermo>
const Foam::speciesTable& Foam::Reaction<ReactionThermo>::gasSpecies() const
{
    NotImplemented;
    return NullObjectRef<speciesTable>();
}


template<class ReactionThermo>
const Foam::List<typename Foam::Reaction<ReactionThermo>::specieCoeffs>&
Foam::Reaction<ReactionThermo>::glhs() const
{
    NotImplemented;
    return List<specieCoeffs>::null();
}


template<class ReactionThermo>
const Foam::List<typename Foam::Reaction<ReactionThermo>::specieCoeffs>&
Foam::Reaction<ReactionThermo>::grhs() const
{
    NotImplemented;
    return List<specieCoeffs>::null();
}


// ************************************************************************* //
