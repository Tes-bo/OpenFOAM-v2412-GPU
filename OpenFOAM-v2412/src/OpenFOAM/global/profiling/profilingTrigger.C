/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2009-2016 Bernhard Gschaider
    Copyright (C) 2016-2024 OpenCFD Ltd.
    Copyright (C) 2023 Josep Pocurull Serra, Barcelona Supercomputing Center
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

\*---------------------------------------------------------------------------*/

#include "profiling.H"
#include "profilingTrigger.H"
#include "profilingInformation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// Extrae profiling hooks
// ----------------------
// https://tools.bsc.es/extrae
// ----------------------

#ifdef HAVE_EXTRAE
#include <map>
#include <utility>

// Weak functions for Extrae C api
extern "C"
{
    typedef unsigned extrae_type_t;
    typedef unsigned long long extrae_value_t;

    // Adds to the Paraver Configuration File human readable information
    // regarding type and its values.
    void Extrae_define_event_type
    (
        extrae_type_t *type,
        char *type_description,
        unsigned *nvalues,
        extrae_value_t *values,
        char **values_description
    ) __attribute__((weak));

    // Adds a single timestamped event into the tracefile
    void Extrae_event
    (
        extrae_type_t type,
        extrae_value_t value
    ) __attribute__((weak));

}  // End extern "C"


// Descriptor for the events
static char myExtrae_description[] = "OpenFOAM Extrae Profiling";

static void open_extrae_region(const std::string& name)
{
    // Event history (by name) of profiling triggers
    static std::map<std::string, extrae_value_t> event_history;

    // Scratch space for transcribing map -> flat lists
    static Foam::DynamicList<char*> eventNames;
    static Foam::DynamicList<extrae_value_t> eventValues;


    if (event_history.empty())
    {
        event_history.insert(std::make_pair("End", 0));
    }

    extrae_type_t event_type = 7000;
    extrae_value_t event_name;

    // Check if there is already an event with that name
    auto iter = event_history.find(name);
    if (iter != event_history.end())
    {
        event_name = iter->second;
    }
    else
    {
        // Update extrae defined events

        event_name = static_cast<extrae_value_t>(event_history.size());
        event_history.insert(std::make_pair(name, event_name));

        unsigned numEvents = event_history.size();

        const Foam::label len(numEvents);

        eventNames.resize_nocopy(len);
        eventValues.resize_nocopy(len);

        Foam::label i = 0;
        for (const auto& iter : event_history)
        {
            eventNames[i]  = const_cast<char*>(iter.first.data());
            eventValues[i] = iter.second;
            ++i;
        }

        Extrae_define_event_type
        (
            &event_type,
            myExtrae_description,
            &numEvents,
            eventValues.data(),
            eventNames.data()
        );
    }

    Extrae_event(event_type, event_name);
}


static void close_extrae_region()
{
    Extrae_event(7000, 0);
}

#endif  /* HAVE_EXTRAE */


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

bool Foam::profilingTrigger::possible() noexcept
{
    return
    (
        profiling::active()
        #ifdef HAVE_EXTRAE
     || bool(Extrae_event)
        #endif
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::profilingTrigger::enter(const std::string& name)
{
    // debug:  std::cerr<< "enter profiling: " << name << '\n';

    ptr_ = profiling::New(name);

    #ifdef HAVE_EXTRAE
    if (Extrae_event)
    {
        open_extrae_region(name);
    }
    #endif
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::profilingTrigger::~profilingTrigger()
{
    stop();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::profilingTrigger::running() const noexcept
{
    return ptr_;
}


void Foam::profilingTrigger::stop()
{
    #ifdef HAVE_EXTRAE
    if (Extrae_event)
    {
        close_extrae_region();
    }
    #endif

    if (ptr_)
    {
        // profiling info pointer managed by pool storage, so no delete here
        profiling::unstack(ptr_);
    }

    ptr_ = nullptr;
}


// ************************************************************************* //
