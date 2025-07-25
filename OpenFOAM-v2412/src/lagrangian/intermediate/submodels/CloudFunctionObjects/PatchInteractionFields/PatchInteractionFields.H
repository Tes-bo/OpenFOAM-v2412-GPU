/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::PatchInteractionFields

Group
    grpLagrangianIntermediateFunctionObjects

Description
    Creates volume fields whose boundaries are used to store patch interaction
    statistics.

    Current field output per patch face:
    - \<cloud\>\<model\>:count - cumulative particle hits
    - \<cloud\>\<model\>:mass - cumulative mass of hitting particles

    Fields can be reset according to:
    - none: fields are not reset
    - timeStep: reset at each time step
    - writeTime: reset at each write time

Usage

    \verbatim
    patchInteractionFields1
    {
        type            patchInteractionFields;
        resetMode       writeTime;
    }
    \endverbatim

SourceFiles
    PatchInteractionFields.C

\*---------------------------------------------------------------------------*/

#ifndef PatchInteractionFields_H
#define PatchInteractionFields_H

#include "CloudFunctionObject.H"
#include "volFields.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class PatchInteractionFields Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class PatchInteractionFields
:
    public CloudFunctionObject<CloudType>
{
public:

    // Public enumerations

        //- Reset mode
        enum class resetMode
        {
            none,
            timeStep,
            writeTime
        };

        static const Enum<resetMode> resetModeNames_;


private:

    // Private Data

        //- Convenience typedef for parcel type
        typedef typename CloudType::parcelType parcelType;

        //- Parcel mass field
        autoPtr<volScalarField> massPtr_;

        //- Parcel count field
        autoPtr<volScalarField> countPtr_;

        //- Reset mode
        resetMode resetMode_;


protected:

    // Protected Member Functions

        //- Create|read|reset the fields
        void reset();

        //- Write post-processing info
        virtual void write();

        //- Helper function to clear or reset fields
        void clearOrReset
        (
            autoPtr<volScalarField>& fieldPtr,
            const word& fieldName,
            const dimensionSet& dims
        ) const;


public:

    //- Runtime type information
    TypeName("patchInteractionFields");


    // Constructors

        //- Construct from dictionary
        PatchInteractionFields
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        PatchInteractionFields(const PatchInteractionFields<CloudType>& pii);

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObject<CloudType>> clone() const
        {
            return autoPtr<CloudFunctionObject<CloudType>>
            (
                new PatchInteractionFields<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~PatchInteractionFields() = default;


    // Member Functions

        // Evaluation

            //- Pre-evolve hook
            virtual void preEvolve
            (
                const typename parcelType::trackingData& td
            );

            //- Post-patch hook
            virtual bool postPatch
            (
                const parcelType& p,
                const polyPatch& pp,
                const typename parcelType::trackingData& td
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PatchInteractionFields.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
