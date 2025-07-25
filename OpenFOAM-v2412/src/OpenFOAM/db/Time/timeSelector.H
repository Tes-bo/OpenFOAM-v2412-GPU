/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2018-2024 OpenCFD Ltd.
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
    Foam::timeSelector

Description
    The timeSelector provides a convenient means of selecting multiple
    times.

    A typical use would be the following:
    \verbatim
    timeSelector::addOptions();
    ...
    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    ...
    forAll(timeDirs, timei)
    {
        ...
    }
    \endverbatim

    With the \c addOptions() method, application receives
    \b -time, \b -latestTime, \b -constant and \b -noZero options.
    The \b -constant option explicitly includes the
    \c constant/ directory in the time list and the \b -noZero option
    explicitly excludes the \c 0/ directory from the time list.
    There may however also be many cases in which neither the \c constant/
    directory nor the \c 0/ directory contain particularly relevant
    information. This might occur, for example, when post-processing
    results. In this case, addOptions is called with optional boolean
    arguments.

    \verbatim
    timeSelector::addOptions(false, true);
    \endverbatim

    The first argument avoids adding the \b -constant option. The second
    argument adds an additional \b -withZero option and also prevents the
    \c 0/ directory from being included in the default time range and in the
    \b -latestTime selection.

    It is also possible to use the timeSelector for setting a single time.
    Typical use would be the following:

    \verbatim
    timeSelector::addOptions_singleTime();
    ...
    #include "setRootCase.H"
    #include "createTime.H"
    timeSelector::setTimeIfPresent(runTime, args);
    \endverbatim

    With the \c addOptions_singleTime() method, application receives
    \b -time, \b -latestTime, \b -constant options. In the case,
    the \b -time option is intended to be a single value and not include
    any ranges.

    The subsequent call to \c setTimeIfPresent() will scan the arguments
    for relevant time options and use them to set the time.

SourceFiles
    timeSelector.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_timeSelector_H
#define Foam_timeSelector_H

#include "scalarRanges.H"
#include "instantList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class argList;
class Time;

/*---------------------------------------------------------------------------*\
                        Class timeSelector Declaration
\*---------------------------------------------------------------------------*/

class timeSelector
{
    // Private Data

        //- The selectable time ranges
        scalarRanges ranges_;


public:

    // Constructors

        //- Default construct
        timeSelector() noexcept = default;

        //- Construct by parsing string for time ranges
        explicit timeSelector(const std::string& str);


    // Member Functions

        //- True if value is within any of the ranges
        bool contains(const scalar value) const;

        //- True if value of the instant is within any of the ranges
        bool contains(const instant& t) const;

        //- True if value of the instant is within any of the ranges
        bool selected(const instant& t) const { return contains(t); }

        //- Return the set of selected instants in the given list that are
        //- within the ranges
        List<bool> selected(const instantList& times) const;

        //- Select a list of Time values that are within the ranges
        instantList select(const instantList& times) const;

        //- Select a list of Time values that are within the ranges
        void inplaceSelect(instantList& times) const;


    // Static Member Functions

        //- Add timeSelector options to argList::validOptions
        //
        //  \par Options added:
        //  - \c -constant
        //  - \c -time
        //  - \c -latestTime
        //  - \c -noZero
        //  - \c -withZero
        //  .
        //
        // \param constant
        //   Add the \b -constant option to include the \c constant/ directory
        //
        // \param withZero
        //   Enable the \b -withZero option and alter the normal time selection
        //   behaviour (and \b -latestTime behaviour) to exclude the \c 0/
        //   directory. The \c 0/ directory will only be included when
        //   \b -withZero is specified.
        //   The \b -noZero option has precedence over the \b -withZero option.
        static void addOptions
        (
            const bool constant=true,
            const bool withZero=false
        );

        //- Add single-time timeSelector options to argList::validOptions()
        //
        //  \par Options added:
        //  - \c -constant
        //  - \c -time
        //  - \c -latestTime
        //  - \c -noZero  (ignored)
        //  .
        static void addOptions_singleTime();

        //- Return the set of times selected based on the argList options
        static instantList select
        (
            const instantList& times,
            const argList& args,
            const word& constantName = "constant"
        );

        //- Return the set of times selected based on the argList options
        //- and also set the runTime to the first instance or the
        //- \c constant/ directory if no instances are specified or available
        static instantList select0
        (
            Time& runTime,
            const argList& args
        );

        //- If any time option provided return the set of times -
        //- as per select0() - otherwise return just the current time.
        //  Also set the runTime to the first instance
        static instantList selectIfPresent
        (
            Time& runTime,
            const argList& args
        );

        //- Set the runTime based on \c -constant (if present),
        //- \c -time (value), or \c -latestTime.
        //  This method is a no-op if no relevant options have been specified.
        static bool setTimeIfPresent
        (
            Time& runTime,
            const argList& args,
            //! Force initial time (default: 0) even if no options specified
            const bool forceInitial = false
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
