/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenCFD Ltd.
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
    Foam::functionObjects::reactionsSensitivityAnalysis

Group
    grpFieldFunctionObjects grpThermophysicalFunctionObjects

Description
    Computes indicators for reaction rates of creation or destruction
    of species in each reaction.

    This function object creates four data files named:

    - \c consumption   :   consumption rate
    - \c production    :   destruction rate
    - \c productionInt  :   integral between dumps of the production rate
    - \c consumptionInt :   integral between dumps of the consumption rate

    Operands:
    \table
      Operand        | Type | Location
      input          | -    | -
      output file    | dat  | $FOAM_CASE/postProcessing/\<FO\>/\<time\>/\<file\>
      output field   | -    | -
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    reactionSensitivityAnalysis1
    {
        // Mandatory entries (unmodifiable)
        type            reactionSensitivityAnalysis;
        libs            (fieldFunctionObjects);

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                             | Type | Req'd | Dflt
      type      | Type name: reactionSensitivityAnalysis  | word |  yes  | -
      libs      | Library name: fieldFunctionObjects      | word |  yes  | -
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link writeFile.H \endlink

    Usage by the \c postProcess utility is not available.

Note
    - Function object only applicable to single cell cases.
    - Needs a \c chemistryModel chosen.

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - Foam::functionObjects::writeFile
    - ExtendedCodeGuide::functionObjects::field::reactionsSensitivityAnalysis

SourceFiles
    reactionsSensitivityAnalysis.C
    reactionsSensitivityAnalysisObjects.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjescts_reactionsSensitivityAnalysis_H
#define functionObjescts_reactionsSensitivityAnalysis_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFieldsFwd.H"
#include "basicChemistryModel.H"
#include "autoPtr.H"
#include "basicMultiComponentMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                Class reactionsSensitivityAnalysis Declaration
\*---------------------------------------------------------------------------*/

template<class chemistryType>
class reactionsSensitivityAnalysis
:
    public fvMeshFunctionObject,
    public writeFile
{
    // Private Data

        //- Number of reactions
        label nReactions_;

        //- Start time of integration
        scalar startTime_;

        //- End time of integration
        scalar endTime_;

        //- List list for species production
        scalarListList production_;

        //- List list for species consumption
        scalarListList consumption_;

        //- List list for species production integral
        scalarListList productionInt_;

        //- List list for species consumption integral
        scalarListList consumptionInt_;

        //- Word list of species
        wordList speciesNames_;

         // File streams

            //- Integrated coefficients
            autoPtr<OFstream> prodFilePtr_;

            //- Moment coefficient
            autoPtr<OFstream> consFilePtr_;

            //- Drag coefficient
            autoPtr<OFstream> prodIntFilePtr_;

            //- Lift coefficient
            autoPtr<OFstream> consIntFilePtr_;


    // Private Member Functions

        //- Create file names for forces and bins
        void createFileNames();

        //- Output file header information
        void writeFileHeader(OFstream& os);

        //- Calculate production and destruction of each species
        void calculateSpeciesRR(const basicChemistryModel&);

        //- Write species production/consumption rates
        void writeSpeciesRR();


public:

    //- Runtime type information
    TypeName("reactionsSensitivityAnalysis");


    // Constructors

        //- Construct from Time and dictionary
        reactionsSensitivityAnalysis
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        reactionsSensitivityAnalysis
        (
            const reactionsSensitivityAnalysis&
        ) = delete;

        //- No copy assignment
        void operator=(const reactionsSensitivityAnalysis&) = delete;


    //- Destructor
    virtual ~reactionsSensitivityAnalysis() = default;


    // Member Functions

        //- Read the reactionsSensitivityAnalysis data
        virtual bool read(const dictionary&);

        //- Execute
        virtual bool execute();

        //- Calculate the reactionsSensitivityAnalysis and write
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "reactionsSensitivityAnalysis.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
