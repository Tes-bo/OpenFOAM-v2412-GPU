/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::fv::interRegionHeatTransferModel

Group
    grpFvOptionsSources

Description
    Intermediate class for handling inter-region heat exchanges.

    The derived classes must provide the heat transfer
    coefficients (htc) which is used as follows in the energy equation.

    \f[
        -htc*Tmapped + Sp(htc, T)
    \f]

Usage
    Minimal example by using \c constant/fvOptions:
    \verbatim
    <userDefinedName1>
    {
        // Mandatory/Optional (inherited) entries
        ...

        // Mandatory entries (unmodifiable)
        fields          (<field1> <field2> ... <fieldN>);
        nbrModel        <nbrModelName>;
        semiImplicit    true;

        // Optional entries (unmodifiable)
        T               <Tname>;
        Tnbr            <TnbrName>;

        // Mandatory/Optional (derived) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                             | Type  | Reqd | Dflt
      nbrModel  | Name of the model in the neighbour mesh | word  | yes  | -
      fields    | Names of operand fields              | wordList | yes  | -
      semiImplicit | Flag to activate semi-implicit coupling | bool | yes | -
      T         | Name of operand temperature field       | word  | no   | T
      Tnbr      | Name of operand neighbour temperature field | word | no | T
    \endtable

    The inherited entries are elaborated in:
     - \link interRegionOption.H \endlink

See also
  - Foam::fv::constantHeatTransfer
  - Foam::fv::tabulatedHeatTransfer
  - Foam::fv::tabulatedNTUHeatTransfer
  - Foam::fv::variableHeatTransfer

SourceFiles
    interRegionHeatTransferModel.C
    interRegionHeatTransferModelTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef interRegionHeatTransferModel_H
#define interRegionHeatTransferModel_H

#include "interRegionOption.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                Class interRegionHeatTransferModel Declaration
\*---------------------------------------------------------------------------*/

class interRegionHeatTransferModel
:
    public interRegionOption
{
protected:

    // Protected Data

        //- Name of the model in the neighbour mesh
        word nbrModelName_;

        //- Pointer to neighbour interRegionHeatTransferModel
        interRegionHeatTransferModel* nbrModel_;

        //- Flag to determine the first iteration
        bool firstIter_;

        //- Flag to activate semi-implicit coupling
        bool semiImplicit_;

        //- Time index - used for updating htc
        label timeIndex_;

        //- Heat transfer coefficient [W/m2/k] times area/volume [1/m]
        volScalarField htc_;

        //- Name of operand temperature field
        word TName_;

        //- Name of operand neighbour temperature field
        word TNbrName_;


    // Protected Member Functions

        //- Set the neighbour interRegionHeatTransferModel
        void setNbrModel();

        //- Inherit correct from interRegionOption
        using interRegionOption::correct;

        //- Correct to calculate the inter-region heat transfer coefficient
        void correct();

        //- Interpolate field with nbrModel specified
        template<class Type>
        tmp<Field<Type>> interpolate
        (
            const interRegionHeatTransferModel& nbrModel,
            const Field<Type>& field
        ) const;

        //- Interpolate field without nbrModel specified
        template<class Type>
        tmp<Field<Type>> interpolate
        (
            const Field<Type>& field
        ) const;

        //- Interpolate field with nbrModel specified
        template<class Type>
        void interpolate
        (
            const interRegionHeatTransferModel& nbrModel,
            const Field<Type>& field,
            Field<Type>& result
        ) const;

        //- Interpolate field without nbrModel specified
        template<class Type>
        void interpolate
        (
            const Field<Type>& field,
            Field<Type>& result
        ) const;


public:

    //- Runtime type information
    TypeName("interRegionHeatTransferModel");


    // Constructors

        //- Construct from dictionary
        interRegionHeatTransferModel
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~interRegionHeatTransferModel() = default;


    // Member Functions

        // Access

            //- Return const access to the neighbour region name
            inline const word& nbrRegionName() const;

            //- Return const access to the mapToMap pointer
            inline const meshToMesh& meshInterp() const;

            //- Return the heat transfer coefficient
            inline const volScalarField& htc() const;

            //- Return const access to the neighbour model
            inline const interRegionHeatTransferModel& nbrModel() const;

            //- Return access to the neighbour model
            inline interRegionHeatTransferModel& nbrModel();

            //- Source term to energy equation
            virtual void addSup
            (
                fvMatrix<scalar>& eqn,
                const label fieldi
            );

            //- Source term to compressible energy equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const label fieldi
            );

            //- Calculate heat transfer coefficient
            virtual void calculateHtc() = 0;


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "interRegionHeatTransferModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "interRegionHeatTransferModelTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
