/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

Class
    Foam::LocalInteraction

Group
    grpLagrangianIntermediatePatchInteractionSubModels

Description
    Patch interaction specified on a patch-by-patch basis

\*---------------------------------------------------------------------------*/

#ifndef LocalInteraction_H
#define LocalInteraction_H

#include "PatchInteractionModel.H"
#include "patchInteractionDataList.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class LocalInteraction Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class LocalInteraction
:
    public PatchInteractionModel<CloudType>
{
    // Private data

        //- List of participating patches
        const patchInteractionDataList patchData_;


        // Bookkeeping for particle fates

            //- Number of parcels escaped
            List<List<label>> nEscape_;

            //- Mass of parcels escaped
            List<List<scalar>> massEscape_;

            //- Number of parcels stuck to patches
            List<List<label>> nStick_;

            //- Mass of parcels stuck to patches
            List<List<scalar>> massStick_;


        //- Flag to output data as fields
        bool writeFields_;

        //- InjectorId to index map, when outputting escaped/stick/...
        //  particles sorted by injectorID
        Map<label> injIdToIndex_;

        //- Mass escape field
        autoPtr<volScalarField> massEscapePtr_;

        //- Mass stick field
        autoPtr<volScalarField> massStickPtr_;


protected:

    // Protected Member Functions

        //- Output file header information
        virtual void writeFileHeader(Ostream& os);


public:

    //- Runtime type information
    TypeName("localInteraction");


    // Constructors

        //- Construct from dictionary
        LocalInteraction(const dictionary& dict, CloudType& owner);

        //- Construct copy from owner cloud and patch interaction model
        LocalInteraction(const LocalInteraction<CloudType>& pim);

        //- Construct and return a clone using supplied owner cloud
        virtual autoPtr<PatchInteractionModel<CloudType>> clone() const
        {
            return autoPtr<PatchInteractionModel<CloudType>>
            (
                new LocalInteraction<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~LocalInteraction() = default;


    // Member Functions

        //- Return access to the massEscape field
        volScalarField& massEscape();

        //- Return access to the massStick field
        volScalarField& massStick();

        //- Apply velocity correction
        //  Returns true if particle remains in same cell
        virtual bool correct
        (
            typename CloudType::parcelType& p,
            const polyPatch& pp,
            bool& keepParticle
        );


        //- Write patch interaction info
        virtual void info();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "LocalInteraction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
