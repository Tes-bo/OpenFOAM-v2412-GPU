/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
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
    Foam::functionObjects::valueAverageBase

Description
    Base class that computes the ensemble- or time-based singular-value average
    values, with optional windowing, from the output of function objects
    that generate non-field type values, e.g. \c Cd of \c forceCoeffs or
    \c momentum_x in \c momentum function objects.

See also
    - Foam::functionObjects::writeFile

SourceFiles
    valueAverageBase.C
    valueAverageBaseTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_valueAverageBase_H
#define functionObjects_valueAverageBase_H

#include "writeFile.H"
#include "Enum.H"
#include "stateFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                      Class valueAverageBase Declaration
\*---------------------------------------------------------------------------*/

class valueAverageBase
:
    public functionObjects::writeFile
{
public:

    // Public Enumerations

        //- Averaging window types
        enum class windowType : uint8_t
        {
            NONE,           //!< "none"
            APPROXIMATE,    //!< "approximate"
            EXACT           //!< "exact"
        };

        //- Names for windowType enumeration
        static const Enum<windowType> windowTypeNames;

private:

    // Private Data

        //- Log flag
        bool log;


protected:

    // Protected Data

        //- Reset the averaging process on restart
        bool resetOnRestart_;

        //- Window type
        windowType windowType_;

        //- Reference to the state functionObject
        stateFunctionObject& state_;

        //- Name of function object to retrieve data from
        word functionObjectName_;

        //- List of fields on which to operate
        wordList fieldNames_;

        //- Optional tolerance to check for converged results
        scalar tolerance_;

        //- Averaging window
        scalar window_;

        //- Average time per field
        List<scalar> totalTime_;


    // Protected Member Functions

        //- Templated function to calculate the average
        //  Return true if processed
        template<class Type, class Type2=Type>
        bool calc(const label fieldi, bool& converged, dictionary& dict);

        //- Output file header information
        virtual void writeFileHeader(Ostream& os) const;

        //- Read state from dictionary
        virtual void readState(dictionary& dict);

        //- Write state to dictionary for restarts
        virtual void writeState(dictionary& dict);

        //- Reset state
        virtual void resetState(dictionary& dict);


public:

    // Constructors

        //- Construct from Time and dictionary
        valueAverageBase
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict,
            stateFunctionObject& state,
            const bool writeToFile = true
        );

        //- No copy construct
        valueAverageBase(const valueAverageBase&) = delete;

        //- No copy assignment
        void operator=(const valueAverageBase&) = delete;


    //- Destructor
    virtual ~valueAverageBase() = default;


    // Public Member Functions

        //- Read the field value average data
        virtual bool read(const dictionary& dict);

        //- Calculate the averages
        virtual bool calculate(dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "valueAverageBaseTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
