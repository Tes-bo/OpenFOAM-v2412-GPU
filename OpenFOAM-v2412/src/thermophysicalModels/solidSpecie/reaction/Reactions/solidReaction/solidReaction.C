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

#include "solidReaction.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::solidReaction<ReactionThermo>::solidReaction
(
    const Reaction<ReactionThermo>& reaction,
    const speciesTable& pyrolisisGases,
    const List<specieCoeffs>& glhs,
    const List<specieCoeffs>& grhs
)
:
    Reaction<ReactionThermo>(reaction),
    pyrolisisGases_(pyrolisisGases),
    glhs_(glhs),
    grhs_(grhs)
{}


template<class ReactionThermo>
Foam::solidReaction<ReactionThermo>::solidReaction
(
    const solidReaction<ReactionThermo>& r,
    const speciesTable& pyrolisisGases
)
:
    Reaction<ReactionThermo>(r),
    pyrolisisGases_(pyrolisisGases),
    glhs_(r.glhs_),
    grhs_(r.grhs_)
{}


template<class ReactionThermo>
Foam::solidReaction<ReactionThermo>::solidReaction
(
    const speciesTable& species,
    const ReactionTable<ReactionThermo>& thermoDatabase,
    const dictionary& dict
)
:
    Reaction<ReactionThermo>
    (
        species,
        thermoDatabase,
        dict,
        false,  // initReactionThermo = false
        false   // failUnknownSpecie = false
    ),
    pyrolisisGases_(dict.parent().parent().lookup("gaseousSpecies")),
    glhs_(),
    grhs_()
{
    ICharStream reactionIs(dict.getString("reaction"));

    this->setLRhs
    (
        reactionIs,
        pyrolisisGases_,
        glhs_,
        grhs_,
        false   // failUnknownSpecie = false
    );

    speciesTable allSpecies(species);
    for (const word& gasName : pyrolisisGases_)
    {
        allSpecies.push_uniq(gasName);
    }
    List<specieCoeffs> dummyLhs;
    List<specieCoeffs> dummyRhs;

    // Rescan (and fail) if a species is neither gas nor solid
    reactionIs.rewind();
    this->setLRhs
    (
        reactionIs,
        allSpecies,
        dummyLhs,
        dummyRhs
        // failUnknownSpecie = true
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ReactionThermo>
const Foam::List<typename Foam::solidReaction<ReactionThermo>::specieCoeffs>&
Foam::solidReaction<ReactionThermo>::glhs() const
{
    return glhs_;
}


template<class ReactionThermo>
const Foam::List<typename Foam::Reaction<ReactionThermo>::specieCoeffs>&
Foam::solidReaction<ReactionThermo>::grhs() const
{
    return grhs_;
}


template<class ReactionThermo>
const Foam::speciesTable& Foam::solidReaction<ReactionThermo>::
gasSpecies() const
{
    return pyrolisisGases_;
}


template<class ReactionThermo>
void Foam::solidReaction<ReactionThermo>::write(Ostream& os) const
{
    OStringStream reaction;
    os.writeEntry("reaction", solidReactionStr(reaction));
}


template<class ReactionThermo>
Foam::string Foam::solidReaction<ReactionThermo>::solidReactionStr
(
    OStringStream& reaction
) const
{
    this->reactionStrLeft(reaction);
    if (!glhs().empty())
    {
        reaction << " + ";
        solidReactionStrLeft(reaction);
    }

    reaction << " = ";

    this->reactionStrRight(reaction);
    if (!grhs().empty())
    {
        reaction << " + ";
        solidReactionStrRight(reaction);
    }
    return reaction.str();
}


template<class ReactionThermo>
void Foam::solidReaction<ReactionThermo>::solidReactionStrLeft
(
    OStringStream& reaction
) const
{
    Reaction<ReactionThermo>::reactionStr(reaction, gasSpecies(), glhs());
}


template<class ReactionThermo>
void Foam::solidReaction<ReactionThermo>::solidReactionStrRight
(
    OStringStream& reaction
) const
{
    Reaction<ReactionThermo>::reactionStr(reaction, gasSpecies(), grhs());
}


// ************************************************************************* //
