/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2023 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Application
    Test-refPtr

Description
    Tests some basic functionality of refPtr

\*---------------------------------------------------------------------------*/

#include "primitiveFields.H"
#include "autoPtr.H"
#include "refPtr.H"
#include "tmp.H"
#include "Switch.H"

using namespace Foam;

struct myScalarField : public scalarField
{
    using scalarField::scalarField;
};


template<class T>
void constructInfo()
{
    Info<< " move-constructible:"
        << std::is_move_constructible<T>::value
        << " move-assignable:"
        << std::is_move_assignable<T>::value
        << " nothrow:"
        << std::is_nothrow_move_assignable<T>::value
        << " trivially:"
        << std::is_trivially_move_assignable<T>::value
        << nl;
}


template<class T>
void printInfo(const autoPtr<T>& item, const bool verbose = false)
{
    Info<< "autoPtr good:" << Switch::name(item.good())
        << " addr: " << Foam::name(item.get());

    constructInfo<autoPtr<T>>();

    if (verbose && item)
    {
        Info<< "content: " << item() << nl;
    }
}


template<class T>
void printInfo(const refPtr<T>& item, const bool verbose = false)
{
    Info<< "refPtr good:" << Switch::name(item.good())
        << " pointer:" << Switch::name(item.is_pointer())
        << " addr: " << Foam::name(item.get())
        << " movable:" << Switch(item.movable());

    constructInfo<refPtr<T>>();

    if (verbose && item)
    {
        Info<< "content: " << item() << nl;
    }
}


template<class T>
void printInfo(const tmp<T>& item, const bool verbose = false)
{
    Info<< "tmp good:" << Switch::name(item.good())
        << " pointer:" << Switch::name(item.is_pointer())
        << " addr: " << Foam::name(item.get())
        << " movable:" << Switch(item.movable());

    constructInfo<tmp<T>>();

    if (verbose && item)
    {
        Info<< "content: " << item() << nl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main()
{
    {
        Info<< nl << "Construct from reference" << nl;

        scalarField f2(10, Foam::sqrt(2.0));
        refPtr<scalarField> tfld2(f2);
        printInfo(tfld2, true);

        Info<< nl << "emplaced:"<< nl;
        tfld2.emplace(25, scalar(1));
        printInfo(tfld2, true);
    }

    {
        Info<< nl << "Construct from New (is_pointer)" << nl;
        auto tfld1 = refPtr<scalarField>::New(10, scalar(-1));
        printInfo(tfld1, true);

        Info<< nl << "emplaced:"<< nl;
        tfld1.emplace(15, scalar(1));
        printInfo(tfld1, true);

        Info<< nl << "Dereferenced: " << *tfld1 << nl;

        Info<< nl << "Construct from autoPtr" << nl;
        refPtr<scalarField> tfld2(autoPtr<scalarField>::New(10, scalar(2)));
        printInfo(tfld2, true);


        Info<< nl << "Construct from unique_ptr" << nl;
        std::unique_ptr<scalarField> ptr(new scalarField(10, scalar(3)));
        refPtr<scalarField> tfld3(std::move(ptr));
        printInfo(tfld3, true);


        Info<< nl << "Reset from autoPtr" << nl;
        tfld2.reset(autoPtr<scalarField>::New(3, scalar(13)));
        printInfo(tfld2, true);


        Info<< nl << "Reset from unique_ptr" << nl;
        ptr.reset(new scalarField(5, scalar(15)));
        tfld3.reset(std::move(ptr));
        printInfo(tfld3, true);


        ptr.reset(new scalarField(2, scalar(1)));
        Info<< nl << "const-ref from pointer: " << name(ptr.get()) << nl;
        tfld3.cref(ptr.get());
        printInfo(tfld3, true);
    }

    {
        refPtr<scalarField> tfld1;
        auto aptr = autoPtr<scalarField>::New(2, scalar(2));

        tmp<scalarField> tfld2;
        printInfo(tfld2, true);

        tfld2 = new scalarField(10, Zero);

    /*
        tfld2 = aptr.get();

        // tfld1.reset(aptr);
        // tfld1 = std::move(aptr);
        // tfld1 = aptr;

        Info<< nl << "From autoPtr" << nl;
        printInfo(aptr, true);
        //& Info<< nl << "Construct from autoPtr" << nl;
        //& // refPtr<scalarField> tfld2(autoPtr<scalarField>::New(10, scalar(2)));
        //& printInfo(tfld2, true);
*/
    }

    {
        auto aptr1 = autoPtr<labelField>::New(2, Zero);
        //auto aptr1 = autoPtr<scalarField>::New(2, scalar(2));
        auto aptr2 = autoPtr<scalarField>::New(2, scalar(2));

        refPtr<scalarField> tfld2(std::move(aptr2));

        // aptr2 = std::move(aptr1);
    }

    {
        auto tptr1 = tmp<labelField>::New(2, Zero);
        auto aptr1 = autoPtr<labelField>::New(2, Zero);
        auto tfld2 = refPtr<labelField>::New(2, Zero);

        // Deleted: refPtr<labelField> tfld1(aptr1);
        refPtr<labelField> tfld1;

        // refPtr<labelField> tfld1(std::move(tptr1));
        // refPtr<labelField> tfld1(tptr1);

        tfld1 = std::move(aptr1);

        // tfld1.reset(aptr1);
        // tfld1.reset(tfld2);

        // tfld1 = std::move(tptr1);
        // Does not compile: tfld1.ref(tptr1);
        // Deleted: tfld1.cref(tptr1);
        // Deleted: tfld1.ref(aptr1);
    }

    Info<< "\nEnd" << endl;

    return 0;
}


// ************************************************************************* //
