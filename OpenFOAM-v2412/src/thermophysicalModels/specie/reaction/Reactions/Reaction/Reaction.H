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

Class
    Foam::Reaction

Description
    Simple extension of ReactionThermo to handle reaction kinetics in addition
    to the equilibrium thermodynamics already handled.

SourceFiles
    ReactionI.H
    Reaction.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_Reaction_H
#define Foam_Reaction_H

#include "speciesTable.H"
#include "HashPtrTable.H"
#include "scalarField.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations

template<class ThermoType>
using ReactionTable = HashPtrTable<ThermoType>;

template<class ReactionThermo>
class Reaction;

template<class ReactionThermo>
inline Ostream& operator<<(Ostream&, const Reaction<ReactionThermo>&);


/*---------------------------------------------------------------------------*\
                           Class Reaction Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo>
class Reaction
:
    public ReactionThermo::thermoType
{
public:

    // Static Data

        //- Number of un-named reactions
        static label nUnNamedReactions;


    // Public Classes

        //- Hold specie index and its coefficients in the
        //- reaction rate expression
        struct specieCoeffs
        {
            label index;
            scalar stoichCoeff;
            scalar exponent;

            specieCoeffs()
            :
                index(-1),
                stoichCoeff(0),
                exponent(1)
            {}

            specieCoeffs
            (
                const speciesTable& species,
                Istream& is,
                bool failUnknownSpecie = true
            );

            bool operator==(const specieCoeffs& sc) const noexcept
            {
                return index == sc.index;
            }

            bool operator!=(const specieCoeffs& sc) const noexcept
            {
                return index != sc.index;
            }

            friend Ostream& operator<<(Ostream& os, const specieCoeffs& sc)
            {
                os  << sc.index << token::SPACE
                    << sc.stoichCoeff << token::SPACE
                    << sc.exponent;
                return os;
            }
        };


private:

    // Private Data

        //- Name of reaction
        const word name_;

        //- List of specie names present in reaction system
        const speciesTable& species_;

        //- Specie info for the left-hand-side of the reaction
        List<specieCoeffs> lhs_;

        //- Specie info for the right-hand-side of the reaction
        List<specieCoeffs> rhs_;


    // Private Member Functions

        //- Return string representation of reaction
        string reactionStr(OStringStream& reaction) const;

        //- Construct reaction thermo
        void setThermo(const ReactionTable<ReactionThermo>& thermoDatabase);

        //- Return new reaction ID for un-named reactions
        label getNewReactionID();

        //- No copy assignment
        void operator=(const Reaction<ReactionThermo>&) = delete;


protected:

    // Protected Member Functions

        //- Add string representation for given species reaction coeffs
        static void reactionStr
        (
            OStringStream& reaction,
            const speciesTable& species,
            const List<specieCoeffs>& reactCoeffs
        );

        //- Add string representation of the left of the reaction
        void reactionStrLeft(OStringStream& reaction) const;

        //- Add string representation of the right of the reaction
        void reactionStrRight(OStringStream& reaction) const;


public:

    //- Runtime type information
    TypeName("Reaction");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            Reaction,
            dictionary,
            (
                const speciesTable& species,
                const ReactionTable<ReactionThermo>& thermoDatabase,
                const dictionary& dict
            ),
            (species, thermoDatabase, dict)
        );


    // Constructors

        //- Construct from components
        Reaction
        (
            const speciesTable& species,
            const List<specieCoeffs>& lhs,
            const List<specieCoeffs>& rhs,
            const ReactionTable<ReactionThermo>& thermoDatabase,
            bool initReactionThermo = true
        );

        //- Construct as copy given new speciesTable
        Reaction
        (
            const Reaction<ReactionThermo>&,
            const speciesTable& species
        );

        //- Construct from dictionary
        // NOTE: initReactionThermo is used by solidReaction where there is no
        // need of setting a lhs - rhs thermo type for each reaction. This is
        // needed for mechanism with reversible reactions
        Reaction
        (
            const speciesTable& species,
            const ReactionTable<ReactionThermo>& thermoDatabase,
            const dictionary& dict,
            bool initReactionThermo = true,
            bool failUnknownSpecie = true
        );

        //- Construct and return a clone
        virtual autoPtr<Reaction<ReactionThermo>> clone() const
        {
            return autoPtr<Reaction<ReactionThermo>>::New(*this);
        }

        //- Construct and return a clone with new speciesTable
        virtual autoPtr<Reaction<ReactionThermo>> clone
        (
            const speciesTable& species
        ) const
        {
            return autoPtr<Reaction<ReactionThermo>>::New(*this, species);
        }


    // Selectors

        //- Return a pointer to new patchField created on freestore from dict
        static autoPtr<Reaction<ReactionThermo>> New
        (
            const speciesTable& species,
            const ReactionTable<ReactionThermo>& thermoDatabase,
            const dictionary& dict
        );


    //- Destructor
    virtual ~Reaction() = default;


    // Member Functions

    // Access

        //- Name of reaction
        const word& name() const noexcept { return name_; }

        //- Access to specie list
        const speciesTable& species() const noexcept { return species_; }

        //- Specie info for the left-hand-side of the reaction
        const List<specieCoeffs>& lhs() const noexcept { return lhs_; }

        //- Specie info for the right-hand-side of the reaction
        const List<specieCoeffs>& rhs() const noexcept { return rhs_; }

        //- Access to gas specie list
        virtual const speciesTable& gasSpecies() const;

        //- Access to gas components of the reaction
        virtual const List<specieCoeffs>& grhs() const;
        virtual const List<specieCoeffs>& glhs() const;


        // Manipulation

            //- Construct the left- and right-hand-side reaction coefficients
            void setLRhs
            (
                Istream&,
                const speciesTable&,
                List<specieCoeffs>& lhs,
                List<specieCoeffs>& rhs,
                bool failUnknownSpecie = true
            );


        // Reaction rate coefficients

            //- Forward rate constant
            virtual scalar kf
            (
                const scalar p,
                const scalar T,
                const scalarField& c
            ) const;

            //- Reverse rate constant from the given forward rate constant
            virtual scalar kr
            (
                const scalar kfwd,
                const scalar p,
                const scalar T,
                const scalarField& c
            ) const;

            //- Reverse rate constant.
            //  Note this evaluates the forward rate constant and divides by the
            //  equilibrium constant
            virtual scalar kr
            (
                const scalar p,
                const scalar T,
                const scalarField& c
            ) const;


        //- Write
        virtual void write(Ostream& os) const;


    // Ostream Operator

        friend Ostream& operator<< <ReactionThermo>
        (
            Ostream&,
            const Reaction<ReactionThermo>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ReactionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Reaction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
