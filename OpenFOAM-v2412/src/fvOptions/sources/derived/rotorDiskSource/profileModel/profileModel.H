/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::profileModel

Description
    Base class for profile models for handling
    aerofoil lift and drag polar diagrams.

Usage
    Minimal example by using \c constant/fvOptions:
    rotorDiskSource1
    {
        // Mandatory/Optional (inherited) entries
        ...

        // Mandatory entries (runtime modifiable)
        profiles
        {
            // Mandatory entries (runtime modifiable)
            <profile1>
            {
                type    <profileModel>;

                // Mandatory/Optional (inherited) entries
                ...
            }
            ...
            <profileN>
            {
                ...
            }

            // Optional entries (runtime modifiable)
            file        <fileName>;
        }
    }

    where the entries mean:
    \table
      Property   | Description                          | Type  | Reqd | Dflt
      type       | Profile model - see below            | word  | yes  | -
      file  | Name of file containing profile characteristics | word   | no  | -
    \endtable

    Options for the \c type entry:
    \verbatim
      lookup    | Polar data is linearly interpolated based on angle of attack
      series    | Polar data is computed as sum of cosine series
    \endverbatim

See also
  - Foam::fv::rotorDiskSource
  - Foam::lookupProfile
  - Foam::seriesProfile

SourceFiles
    profileModel.C
    profileModelList.C

\*---------------------------------------------------------------------------*/

#ifndef profileModel_H
#define profileModel_H

#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class profileModel Declaration
\*---------------------------------------------------------------------------*/

class profileModel
{
protected:

    // Protected Data

        //- Coefficients dictionary
        const dictionary dict_;

        //- Name of profile model
        const word name_;

        //- File name (optional)
        fileName fName_;


    // Protected Member Functions

        //- Return true if file name is set
        bool readFromFile() const;


public:

    //- Runtime type information
    TypeName("profileModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            profileModel,
            dictionary,
            (
                const dictionary& dict,
                const word& modelName
            ),
            (dict, modelName)
        );


    // Selectors

        //- Return a reference to the selected fvOption model
        static autoPtr<profileModel> New(const dictionary& dict);


    //- Constructor
    profileModel(const dictionary& dict, const word& modelName);


    //- Destructor
    virtual ~profileModel() = default;


    // Member Functions

        // Access

            //- Return const access to the source name
            const word& name() const;


        // Evaluation

            //- Return the Cd and Cl for a given angle-of-attack
            virtual void Cdl
            (
                const scalar alpha,
                scalar& Cd,
                scalar& Cl
            ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
