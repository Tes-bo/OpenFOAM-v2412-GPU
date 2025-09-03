/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2020 OpenCFD Ltd.
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
    Foam::functionObjects::fieldAverageItem

Description
    Helper class to describe what form of averaging to apply.  A set will be
    applied to each base field in Foam::fieldAverage, of the following form.

Usage
    \verbatim
    <field1>
    {
        mean            on;
        prime2Mean      on;
        base            time; // iteration
        window          200;  // optional averaging window
        windowName      w1;   // optional window name (default = "")
        windowType      approximate; // window type

        allowRestart    yes;  // optional, used for windowType 'exact'
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Req'd | Dflt
      mean         | Flag to calculate average          | bool |  yes  | -
      prime2Mean   | Flag to calculate prime-square average | bool |  yes  | -
      base         | Type of averaging interval         | word | yes   | -
      window       | Averaging window                   | scalar | no  |
      windowName   | Name of the averaging window       | word   | no  | ""
      windowType   | Type of averaging window           | word   | no  |
      allowRestart | Flag to allow restart for windowType=exact | bool | no |
    \endtable

    Options for the \c base entry:
    \verbatim
      time | Averaging interval is based on time
      iter | Averaging interval is based on iterations
    \endverbatim

    Options for the \c windowType entry:
    \verbatim
      none  | no windowing
      exact | allow additional files will be stored and written
      approximate | disallow additional files will be stored and written
    \endverbatim

Note
    To employ the \c prime2Mean option, the \c mean option must be enabled.

SourceFiles
    fieldAverageItem.C
    fieldAverageItemIO.C

\*---------------------------------------------------------------------------*/

#ifndef fieldAverageItem_H
#define fieldAverageItem_H

#include "Enum.H"
#include "Switch.H"
#include "FIFOStack.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class objectRegistry;

namespace functionObjects
{

// Forward Declarations
class fieldAverageItem;
Istream& operator>>(Istream&, fieldAverageItem&);
Ostream& operator<<(Ostream&, const fieldAverageItem&);

/*---------------------------------------------------------------------------*\
                      Class fieldAverageItem Declaration
\*---------------------------------------------------------------------------*/

class fieldAverageItem
{
public:

    // Public Data

        // File and field name extensions

            //- Mean average
            static const word EXT_MEAN;

            //- Prime-squared average
            static const word EXT_PRIME2MEAN;


        //- Enumeration defining the averaging base type
        enum class baseType
        {
            ITER,
            TIME
        };

        //- Enumeration defining the averaging window type
        enum class windowType
        {
            NONE,
            APPROXIMATE,
            EXACT
        };


private:

    // Private Data

        //- Active flag
        bool active_;

        //- Field name
        word fieldName_;

        //- Compute mean flag
        bool mean_;

        //- Name of mean field
        word meanFieldName_;

        //- Compute prime-squared mean flag
        bool prime2Mean_;

        //- Name of prime-squared mean field
        word prime2MeanFieldName_;

        //- Averaging base type names
        static const Enum<baseType> baseTypeNames_;

        //- Averaging base type
        baseType base_;

        //- Total number of iterations item has been evolved
        label totalIter_;

        //- Total time item has been evolved
        scalar totalTime_;

        //- Averaging window - defaults to -1 for 'all iters/time'
        scalar window_;

        //- Averaging window name - defaults to 'window'
        word windowName_;

        //- Averaging window type
        windowType windowType_;

        //- Averaging window type names
        static const Enum<windowType> windowTypeNames_;

        //- List of window times (windowType = EXACT)
        FIFOStack<scalar> windowTimes_;

        //- List of window field names (windowType = EXACT)
        FIFOStack<word> windowFieldNames_;

        //- Switch to write all necessary files for clean restart
        bool allowRestart_;


public:

    // Constructors

        //- Construct null
        fieldAverageItem();

        //- Construct from Istream
        fieldAverageItem(Istream&);

        //- Construct as copy
        fieldAverageItem(const fieldAverageItem&);


    //- Destructor
    ~fieldAverageItem();


    // Member Functions

        // Access

            //- Return const access to the active flag
            inline bool active() const;

            //- Return non-const access to the active flag
            inline bool& active();

            //- Return const access to the field name
            inline const word& fieldName() const;

            //- Return const access to the mean flag
            inline bool mean() const;

            //- Return non-const access to the mean flag
            inline bool& mean();

            //- Return const access to the mean field name
            inline const word& meanFieldName() const;

            //- Set the mean field name
            inline void setMeanFieldName(const word& name);

            //- Return const access to the prime-squared mean flag
            inline bool prime2Mean() const;

            //- Return non-const access to the prime-squared mean flag
            inline bool& prime2Mean();

            //- Return const access to the prime-squared mean field name
            inline const word& prime2MeanFieldName() const;

            //- Set the prime-squared mean field name
            inline void setPrime2MeanFieldName(const word& name);

            //- Return averaging base type name
            inline const word& base() const;

            //- Return the total number of iterations item has been evolved
            inline label totalIter() const;

            //- Return the total time item has been evolved
            inline scalar totalTime() const;

            //- Return the window length (iterations or seconds)
            inline scalar window() const;

            //- Return the (optional) window name
            inline const word& windowName() const;

            //- Return the list of window times (windowType = EXACT)
            inline const FIFOStack<scalar>& windowTimes() const;

            //- Return the list of window field names (windowType = EXACT)
            inline const FIFOStack<word>& windowFieldNames() const;

            //- Return the allow restart flag
            inline bool allowRestart() const;

            //- Return the current time interval
            inline scalar dt(const scalar deltaT) const;

            //- Return the total time interval
            inline scalar Dt() const;

            //- Helper function to construct a window field name
            inline word windowFieldName(const word& prefix) const;

            //- Return true if time is inside window (including boundaries)
            inline bool inWindow(const scalar t) const;

            //- Return true if we wish to store window fields
            inline bool storeWindowFields() const;

            //- Return true if we wish to write window fields
            inline bool writeWindowFields() const;

            //- Add field to window
            void addToWindow(const word& fieldName, const scalar deltaT);

            //- Evolve and update
            void evolve(const objectRegistry& obr);

            //- Clear out all mean fields and (optionally) supporting data
            void clear(const objectRegistry& obr, const bool fullClean);

            //- Read state and re-initialise values
            bool readState(const dictionary& dict);

            //- Write state for restart
            void writeState(dictionary& dict) const;

            //- Calculate the mean field value
            template<class Type>
            bool calculateMeanField(const objectRegistry& obr) const;

            //- Calculate prime-squared average fields
            template<class Type1, class Type2>
            bool calculatePrime2MeanField(const objectRegistry& obr) const;


    // Member Operators

        void operator=(const fieldAverageItem&);


    // Friend Operators

        friend bool operator==
        (
            const fieldAverageItem& a,
            const fieldAverageItem& b
        )
        {
            return
                a.active_ == b.active_
             && a.fieldName_ == b.fieldName_
             && a.mean_ == b.mean_
             && a.meanFieldName_ == b.meanFieldName_
             && a.prime2Mean_ == b.prime2Mean_
             && a.prime2MeanFieldName_ == b.prime2MeanFieldName_
             && a.base_ == b.base_
             && a.totalIter_ == b.totalIter_
             && a.totalTime_ == b.totalTime_
             && a.window_ == b.window_
             && a.windowName_ == b.windowName_
             && a.windowType_ == b.windowType_
             && a.allowRestart_ == b.allowRestart_;
        }

        friend bool operator!=
        (
            const fieldAverageItem& a,
            const fieldAverageItem& b
        )
        {
            return !(a == b);
        }


    // IOstream Operators

        friend Istream& operator>>(Istream&, fieldAverageItem&);
        friend Ostream& operator<<(Ostream&, const fieldAverageItem&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fieldAverageItemI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fieldAverageItemTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
