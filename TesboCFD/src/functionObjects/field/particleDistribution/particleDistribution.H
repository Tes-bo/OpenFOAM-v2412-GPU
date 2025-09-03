/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2022 OpenCFD Ltd.
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
    Foam::functionObjects::particleDistribution

Group
    grpFieldFunctionObjects

Description
    Generates a particle distribution for lagrangian data at a given time.

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
    particleDistribution1
    {
        // Mandatory entries (unmodifiable)
        type        particleDistribution;
        libs        (fieldFunctionObjects);

        // Mandatory entries (runtime modifiable)
        cloud       <cloudName>;
        nameVsBinWidth
        (
            (d 0.1)
            (U 10)
        );
        setFormat   raw;

        // Optional entries (runtime modifiable)
        tagField    none;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Req'd | Dflt
      type         | Type name: particleDistribution    | word |  yes  | -
      libs         | Library name: fieldFunctionObjects | word |  yes  | -
      cloud        | Name of cloud to process           | word |  yes  | -
      nameVsBinWidth | List of cloud field-bin width | wordHashTable | yes | -
      setFormat      | Output format                 | word          |  yes  | -
      tagField | Name of cloud field to use group particles | word | no | none
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link writeFile.H \endlink

    Usage by the \c postProcess utility is not available.

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - Foam::functionObjects::writeFile
    - ExtendedCodeGuide::functionObjects::field::particleDistribution

SourceFiles
    particleDistribution.C
    particleDistributionTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_particleDistribution_H
#define functionObjects_particleDistribution_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "scalarField.H"
#include "Random.H"
#include "Tuple2.H"
#include "coordSetWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                    Class particleDistribution Declaration
\*---------------------------------------------------------------------------*/

class particleDistribution
:
    public functionObjects::fvMeshFunctionObject,
    public writeFile
{
protected:

    // Protected Data

        //- Cloud name
        word cloudName_;

        //- Tag field name - used to filter the particles into groups
        word tagFieldName_;

        //- Random number generator - used by distribution models
        Random rndGen_;

        //- List of field name vs. bin width
        List<Tuple2<word, scalar>> nameVsBinWidth_;

        //- Writer
        mutable autoPtr<coordSetWriter> writerPtr_;


    // Protected Member Functions

        //- Generate the distribution
        void generateDistribution
        (
            const word& fieldName,
            const scalarField& field,
            const scalar binWidth,
            const label tag = -1
        );

        // Helper function to retrieve the particle data
        template<class Type>
        bool processField
        (
            const objectRegistry& obr,
            const label fieldi,
            const List<DynamicList<label>>& addr
        );


public:

    //- Runtime type information
    TypeName("particleDistribution");


    // Constructors

        //- Construct from Time and dictionary
        particleDistribution
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        particleDistribution(const particleDistribution&) = delete;

        //- No copy assignment
        void operator=(const particleDistribution&) = delete;


    //- Destructor
    virtual ~particleDistribution() = default;


    // Member Functions

        //- Read the particleDistribution data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the particleDistribution
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "particleDistributionTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
