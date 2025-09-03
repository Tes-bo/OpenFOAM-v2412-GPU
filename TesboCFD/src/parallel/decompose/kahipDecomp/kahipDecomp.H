/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2023 OpenCFD Ltd.
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
    Foam::kahipDecomp

Description
    Domain decomposition using KaHIP
    http://algo2.iti.kit.edu/documents/kahip/

    When run in parallel will collect the entire graph on to the master,
    decompose and send back.

    Coefficients dictionary: \a kahipCoeffs, \a coeffs.

    \verbatim
    numberOfSubdomains   N;
    method               kahip;

    kahipCoeffs
    {
        config          fast;
        imbalance       0.01;
    }
    \endverbatim

    Method coefficients:
    \table
        Property  | Description                          | Required | Default
        config    | fast / eco / strong                  | no       | fast
        imbalance | imbalance on cells between domains   | no       | 0.01
        seed      | initial value for random number generator | no  | 0
    \endtable

SourceFiles
    kahipDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_kahipDecomp_H
#define Foam_kahipDecomp_H

#include "metisLikeDecomp.H"
#include "Enum.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class kahipDecomp Declaration
\*---------------------------------------------------------------------------*/

class kahipDecomp
:
    public metisLikeDecomp
{
protected:

    // Protected Member Functions

        //- Call kahip with options from dictionary.
        virtual label decomposeSerial
        (
            const labelList& adjncy,
            const labelList& xadj,
            const List<scalar>& cellWeights,
            labelList& decomp
        ) const;


public:

    //- The predefined KaHIP configuration types
    enum class configs
    {
        FAST = 0, //!< default
        ECO = 1,
        STRONG = 2,
        FASTSOCIAL = 3,
        ECOSOCIAL = 4,
        STRONGSOCIAL = 5,
    };


    //- The selection names for predefined KaHIP configurations
    static const Enum<configs> configNames;


    // Generated Methods

        //- No copy construct
        kahipDecomp(const kahipDecomp&) = delete;

        //- No copy assignment
        void operator=(const kahipDecomp&) = delete;


    //- Runtime type information
    TypeName("kahip");


    // Constructors

        //- Construct with number of domains (no coefficients or constraints)
        explicit kahipDecomp(const label numDomains);

        //- Construct given decomposition dictionary and optional region name
        explicit kahipDecomp
        (
            const dictionary& decompDict,
            const word& regionName = ""
        );


    //- Destructor
    virtual ~kahipDecomp() = default;


    // Member Functions

        //- Method is parallel aware
        virtual bool parallelAware() const
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
