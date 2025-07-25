/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2009-2016 Bernhard Gschaider
    Copyright (C) 2016-2023 OpenCFD Ltd.
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

#include "argList.H"
#include "profiling.H"
#include "profilingInformation.H"
#include "profilingSysInfo.H"
#include "cpuInfo.H"
#include "memInfo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::profiling::allowed(Foam::debug::infoSwitch("allowProfiling", 1));

std::unique_ptr<Foam::profiling> Foam::profiling::singleton_(nullptr);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::profilingInformation* Foam::profiling::create()
{
    // Top-level entry: reset everything
    pool_.clear();
    children_.clear();
    stack_.clear();
    times_.clear();

    Information* info = new Information;

    pool_.push_back(info);
    children_.resize(pool_.size());
    children_.back().clear();  // safety

    return info;
}


Foam::profilingInformation* Foam::profiling::create
(
    profilingInformation *parent,
    const std::string& descr
)
{
    const label parentId = parent->id();

    for (Information* child : children_[parentId])
    {
        if (descr == child->description())
        {
            return child;  // Found existing
        }
    }

    Information* info = new Information(parent, descr, pool_.size());

    pool_.push_back(info);
    children_.resize(pool_.size());
    children_.back().clear();  // safety
    children_[parentId].push_back(info);

    return info;
}


void Foam::profiling::beginTimer(profilingInformation *info)
{
    stack_.push_back(info);
    times_.push_back(clockValue::now());
    info->setActive(true);              // Mark as on stack
}


Foam::profilingInformation* Foam::profiling::endTimer()
{
    Information *info = stack_.back();
    clockValue clockval = times_.back();
    stack_.pop_back();
    times_.pop_back();

    info->update(clockval.elapsed());   // Update elapsed time
    info->setActive(false);             // Mark as off stack

    return info;
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

bool Foam::profiling::active() noexcept
{
    return allowed && singleton_;
}


void Foam::profiling::disable() noexcept
{
    allowed = 0;
}


bool Foam::profiling::print(Ostream& os)
{
    if (allowed && singleton_)
    {
        return singleton_->writeData(os);
    }

    return false;
}


bool Foam::profiling::writeNow()
{
    if (allowed && singleton_)
    {
        return singleton_->regIOobject::write();
    }

    return false;
}


void Foam::profiling::initialize
(
    const IOobject& ioObj,
    const Time& owner
)
{
    if (allowed && !singleton_)
    {
        singleton_.reset(new profiling(ioObj, owner));
    }
}


void Foam::profiling::initialize
(
    const dictionary& dict,
    const IOobject& ioObj,
    const Time& owner
)
{
    if (allowed && !singleton_)
    {
        singleton_.reset(new profiling(dict, ioObj, owner));
    }
}


void Foam::profiling::stop(const Time& owner)
{
    if (singleton_ && &owner == &(singleton_->owner_))
    {
        singleton_.reset(nullptr);
    }
}


Foam::profilingInformation* Foam::profiling::New(const std::string& descr)
{
    Information *info = nullptr;

    if (active())
    {
        Information *parent = singleton_->stack_.back();

        info = singleton_->create(parent, descr);
        singleton_->beginTimer(info);

        if (singleton_->memInfo_)
        {
            singleton_->memInfo_->update();

            info->maxMem_ = Foam::max
            (
                info->maxMem_,
                singleton_->memInfo_->size()
            );
        }
    }

    return info;
}


void Foam::profiling::unstack(const profilingInformation *info)
{
    if (active() && info)
    {
        Information *top = singleton_->endTimer();

        if (info->id() != top->id())
        {
            FatalErrorInFunction
                << "Profiling information to unstack has different id than"
                << " the top of the profiling stack" << nl
                << "  info: " << info->id() << " (" << info->description()
                << ")\n"
                << "  top:  " << top->id()  << " (" << top->description()
                << ")\n" << endl
                << abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::profiling::profiling
(
    const IOobject& io,
    const Time& owner,
    const bool allEnabled
)
:
    IOdictionary(io),
    owner_(owner)
{
    if (allEnabled)
    {
        sysInfo_.reset(new profilingSysInfo);
        cpuInfo_.reset(new cpuInfo);
        memInfo_.reset(new memInfo);
    }

    Information *info = this->create();
    this->beginTimer(info);

    DetailInfo << "profiling initialized" << nl;
}


Foam::profiling::profiling
(
    const dictionary& dict,
    const IOobject& io,
    const Time& owner
)
:
    profiling(io, owner, false)
{
    {
        bool on = false;

        if (dict.readIfPresent("sysInfo", on) && on)
        {
            sysInfo_.reset(new profilingSysInfo);
        }
        if (dict.readIfPresent("cpuInfo", on) && on)
        {
            cpuInfo_.reset(new cpuInfo);
        }
        if (dict.readIfPresent("memInfo", on) && on)
        {
            memInfo_.reset(new memInfo);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::profiling::~profiling()
{
    if (this == singleton_.get())
    {
        singleton_.reset(nullptr);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::Time& Foam::profiling::owner() const noexcept
{
    return owner_;
}


Foam::label Foam::profiling::size() const noexcept
{
    return stack_.size();
}


bool Foam::profiling::writeData(Ostream& os) const
{
    static DynamicList<scalar> elapsed;

    const clockValue now(clockValue::now());

    const label nstack = stack_.size();

    elapsed.resize(nstack+1);   // extend for last entry, which has no child.

    for (label stacki=0; stacki < nstack; ++stacki)
    {
        elapsed[stacki] = (now - times_[stacki]);
    }
    elapsed.back() = 0;

    os.beginBlock("profiling");

    // Active items
    for (label stacki=0; stacki < nstack; ++stacki)
    {
        if (stacki) os << nl;   // Extra line between entries

        stack_[stacki]->write
        (
            os,
            true,
            elapsed[stacki],    // elapsedTime
            elapsed[stacki+1]   // childTimes
        );
    }

    // Non-active items
    for (const Information& info : pool_)
    {
        if (!info.active())
        {
            os << nl;
            info.write(os);
        }
    }

    os.endBlock();

    if (sysInfo_)
    {
        os << nl;
        sysInfo_->writeEntry("sysInfo", os);
    }

    if (cpuInfo_)
    {
        os << nl;
        cpuInfo_->writeEntry("cpuInfo", os);
    }

    if (memInfo_)
    {
        memInfo_->update();
        os << nl;
        memInfo_->writeEntry("memInfo", os);
    }

    return os.good();
}


bool Foam::profiling::writeObject
(
    IOstreamOption,
    const bool writeOnProc
) const
{
    return
        regIOobject::writeObject
        (
            IOstreamOption(IOstreamOption::ASCII),
            true  // always writeOnProc
        );
}


// ************************************************************************* //
