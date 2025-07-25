/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::NoHeatTransfer

Group
    grpLagrangianIntermediateHeatTransferSubModels

Description
    Nusselt-number model providing an interface to the
    properties of \c HeatTransferModel without any active model.

Usage
    Minimal example by using \c constant/\<CloudProperties\>:
    \verbatim
    subModels
    {
        heatTransferModel    none;
    }
    \endverbatim

    where the entries mean:
    \table
      Property           | Description              | Type | Reqd | Deflt
      heatTransferModel  | Type name: none          | word | yes  | -
    \endtable

SourceFiles
    NoHeatTransfer.C

\*---------------------------------------------------------------------------*/

#ifndef NoHeatTransfer_H
#define NoHeatTransfer_H

#include "HeatTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class NoHeatTransfer Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NoHeatTransfer
:
    public HeatTransferModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("none");


    // Generated Methods

        //- No copy assignment
        void operator=(const NoHeatTransfer&) = delete;


    // Constructors

        //- Construct from dictionary
        NoHeatTransfer(const dictionary&, CloudType& owner);

        //- Copy construct
        NoHeatTransfer(const NoHeatTransfer<CloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<HeatTransferModel<CloudType>> clone() const
        {
            return autoPtr<HeatTransferModel<CloudType>>
            (
                new NoHeatTransfer<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~NoHeatTransfer() = default;


    // Member Functions

        //- Flag to indicate whether model activates heat transfer model
        virtual bool active() const
        {
            return false;
        }


        // Evaluation

            //- Return Nusselt number
            virtual scalar Nu(const scalar Re, const scalar Pr) const
            {
                return 0.0;
            }

            //- Return Prandtl number
            virtual scalar Pr() const
            {
                return 1.0;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "NoHeatTransfer.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
