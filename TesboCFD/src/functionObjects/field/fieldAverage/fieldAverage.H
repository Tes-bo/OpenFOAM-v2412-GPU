/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
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
    Foam::functionObjects::fieldAverage

Group
    grpFieldFunctionObjects

Description
    Computes ensemble- and/or time-based field averages, with optional
    windowing, for a user-specified selection of volumetric and/or surface
    fields.

    Fields are entered as a list of sub-dictionaries, which indicate the type of
    averages to perform, and can be updated during the calculation.  The current
    options include:
    - \c mean: arithmetic mean
        \f[
            \overline{x} = \frac{1}{N}\displaystyle\sum\limits_{i=0}^N x_i
        \f]
    - \c prime2Mean: prime-squared mean
        \f[
            \overline{x'}^2 = \frac{1}{N}\displaystyle\sum\limits_{i=0}^N
            (x_i - \overline{x})^2
        \f]
    - \c base: average over 'time', or 'iteration' (\c N in the above)
    - \c window: optional averaging window, specified in 'base' units

    Average field names are constructed by concatenating the base field with
    the averaging type, e.g. when averaging field 'U', the name of resultant
    fields becomes:
    - arithmetic mean field, \c UMean
    - prime-squared field, \c UPrime2Mean

    Information regarding the number of averaging steps, and total averaging
    time are written on a per-field basis to the
    \c "<functionObject name>Properties" dictionary,
    located in \c \<time\>/uniform.

    When restarting form a previous calculation, the averaging is continuous or
    may be restarted using the \c restartOnRestart option.

    The averaging process may be restarted after each calculation output time
    using the \c restartOnOutput option or restarted periodically using the \c
    periodicRestart option and setting \c restartPeriod to the required
    averaging period.

    With the \c subRegion option, also supports fields on function object
    surface output (e.g., \c sampledSurfaces).

    Operands:
    \table
      Operand      | Type                          | Location
      input        | {vol,surface}\<Type\>Field(s) <!--
               --> |$FOAM_CASE/\<time\>/\<inpField\>s
      output file  | -                             | -
      output field | {vol,surface}\<Type\>Field(s) <!--
               --> | $FOAM_CASE/\<time\>/\<outField\>s
    \endtable

    where \c \<Type\>=Scalar/Vector/SphericalTensor/SymmTensor/Tensor.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    fieldAverage1
    {
        // Mandatory entries (unmodifiable)
        type                fieldAverage;
        libs                (fieldFunctionObjects);

        // Mandatory entries (runtime modifiable)
        fields
        (
            <field1>
            {
                // Optional (inherited) entries
                ...
            }

            ...

            <fieldN>
            {
                ...
            }
        );

        // Optional entries (runtime modifiable)
        restartOnRestart    false;
        restartOnOutput     false;
        periodicRestart     false;
        restartPeriod       0.002;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Req'd | Dflt
      type         | Type name: fieldAverage            | word |  yes  | -
      libs         | Library name: fieldFunctionObjects | word |  yes  | -
      fields       | Names of the operand fields and averaging options <!--
               --> | dict |  yes  | -
      restartOnRestart| Restart the averaging on restart | bool | no     | false
      restartOnOutput | Restart the averaging on output  | bool | no     | false
      periodicRestart | Periodically restart the averaging | bool | no   | false
      restartPeriod   | Periodic restart period   | scalar | conditional | -
      restartTime     | One-shot reset of the averaging  | scalar | no   | GREAT
      subRegion    | Name for alternative objectRegistry | word   | no   | ""
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link fieldAverageItem.H \endlink

    Usage by the \c postProcess utility is not available.

See also
    - Foam::functionObject
    - Foam::functionObjects::fieldAverageItem
    - Foam::functionObjects::fvMeshFunctionObject
    - ExtendedCodeGuide::functionObjects::field::fieldAverage

SourceFiles
    fieldAverage.C
    fieldAverageTemplates.C
    fieldAverageItem.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_fieldAverage_H
#define functionObjects_fieldAverage_H

#include "fvMeshFunctionObject.H"
#include "FIFOStack.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

// Forward declaration of classes
class fieldAverageItem;

/*---------------------------------------------------------------------------*\
                        Class fieldAverage Declaration
\*---------------------------------------------------------------------------*/

class fieldAverage
:
    public fvMeshFunctionObject
{
protected:

    // Protected Data

        //- Time at last call, prevents repeated averaging
        label prevTimeIndex_;

        //- Initialised flag
        bool initialised_;

        //- Restart the averaging process on restart
        Switch restartOnRestart_;

        //- Restart the averaging process on output
        Switch restartOnOutput_;

        //- Periodically restart the averaging process
        Switch periodicRestart_;

        //- Restart period
        scalar restartPeriod_;

        //- Specific restart time
        scalar restartTime_;

        //- List of field average items, describing what averages to be
        //  calculated and output
        List<fieldAverageItem> faItems_;

        // Counters

            //- Iteration steps counter
            List<label> totalIter_;

            //- Total time counter
            List<scalar> totalTime_;

            //- Index for periodic restart
            label periodIndex_;


    // Protected Member Functions

        // Initialisation routines

            //- Reset lists (clear existing values) and initialize averaging.
            //  Check requested field averages are valid, populate field lists
            void initialize();

            //- Restart averaging for restartOnOutput
            void restart();

            //- Add mean average field to database
            template<class Type>
            bool addMeanFieldType(fieldAverageItem& item);

            //- Add mean average field to database
            template<class Type>
            bool addMeanField(fieldAverageItem& item);

            //- Add prime-squared average field to database
            template<class Type1, class Type2>
            bool addPrime2MeanFieldType(fieldAverageItem& item);

            //- Add prime-squared average field to database
            template<class Type1, class Type2>
            bool addPrime2MeanField(fieldAverageItem& item);


        // Calculation functions

            //- Main calculation routine
            virtual void calcAverages();

            //- Calculate mean average fields
            template<class Type>
            void calculateMeanFields() const;

            //- Calculate prime-squared average fields
            template<class Type1, class Type2>
            void calculatePrime2MeanFields() const;

            //- Add mean-squared field value to prime-squared mean field
            template<class Type1, class Type2>
            bool addMeanSqrToPrime2MeanType(const fieldAverageItem& item) const;

            //- Add mean-squared field value to prime-squared mean field
            template<class Type1, class Type2>
            void addMeanSqrToPrime2Mean() const;

            template<class Type>
            bool storeWindowFieldType(fieldAverageItem& item);

            template<class Type>
            void storeWindowFields();

            template<class Type>
            bool restoreWindowFieldsType(const fieldAverageItem& item);

            template<class Type>
            void restoreWindowFields(const fieldAverageItem& item);

        // I-O

            //- Write averages
            virtual void writeAverages() const;

            //- Write fields
            template<class Type>
            bool writeFieldType(const word& fieldName) const;

            //- Write fields
            template<class Type>
            void writeFields() const;

            //- Write averaging properties - steps and time
            void writeAveragingProperties();

            //- Read averaging properties - steps and time
            void readAveragingProperties();


public:

    //- Runtime type information
    TypeName("fieldAverage");


    // Constructors

        //- Construct from Time and dictionary
        fieldAverage
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- No copy construct
        fieldAverage(const fieldAverage&) = delete;

        //- No copy assignment
        void operator=(const fieldAverage&) = delete;


    //- Destructor
    virtual ~fieldAverage() = default;


    // Member Functions

        //- Read the field average data
        virtual bool read(const dictionary&);

        //- Calculate the field averages
        virtual bool execute();

        //- Write the field averages
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fieldAverageTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
