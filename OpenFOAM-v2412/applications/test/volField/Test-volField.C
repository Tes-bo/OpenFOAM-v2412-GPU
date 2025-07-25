/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

Application
    Test-volField

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "GeometricFields.H"
#include "transformGeometricField.H"

// #undef TEST_UINT8_FIELD

#ifdef TEST_UINT8_FIELD
namespace Foam
{
    // Something like a state field. Probably only dimensionless
    // - still needs some basic boundary conditions!!
    typedef GeometricField<uint8_t, fvPatchField, volMesh> volUint8Field;

    defineTemplateTypeNameAndDebug(volUint8Field, 0);

} // End namespace Foam
#endif


template<class GeoField>
void setDims(UPtrList<GeoField>& list, const dimensionSet& ds)
{
    forAll(list, i)
    {
        if (list.set(i))
        {
            list[i].dimensions().reset(ds);
        }
    }
}


template<class GeoField>
void rename(UPtrList<GeoField>& list, const word& baseName)
{
    forAll(list, i)
    {
        if (list.set(i))
        {
            list[i].rename(IOobject::groupName(baseName, name(i)));
        }
    }
}


template<class GeoField, class GeoField2>
void setDims(UPtrList<GeoField>& list, const GeoField2& refField)
{
    setDims(list, refField.dimensions());
}


template<class GeoField>
void write
(
    UPtrList<GeoField>& list,
    label len
)
{
    len = min(len, list.size());

    for (label i=0; i < len; ++i)
    {
        if (list.set(i))
        {
            list[i].write();
        }
    }
}


template<class GeoField>
PtrList<GeoField> create
(
    const fvMesh& mesh,
    const label len
)
{
    PtrList<GeoField> list(len);

    forAll(list, i)
    {
        list.set
        (
            i,
            GeoField::New
            (
                word("cmpt." + name(i)),
                mesh,
                Foam::zero{},  // value
                dimless
            ).ptr()
        );
    }

    return list;
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addBoolOption("zip", "run zip/unzip tests");

    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"

    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< nl
        << "p.v().size(): "
        << p.v().size() << endl;

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    #include "createPhi.H"

    volSymmTensorField st
    (
        IOobject
        (
            "st",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensioned<symmTensor>(dimless, symmTensor::one),
        fvPatchFieldBase::zeroGradientType()
    );

    volTensorField tensf
    (
        IOobject
        (
            "tf",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensioned<tensor>(dimless, tensor(1,2,3,4,5,6,7,8,9)),
        fvPatchFieldBase::zeroGradientType()
    );

    SolverPerformance<symmTensor> sP =
    (
        solve
        (
            fvm::ddt(st)
          + fvm::div(phi, st)
          - fvm::laplacian
            (
                dimensionedScalar("D", sqr(dimLength)/dimTime, 1),
                st
            )
         ==
            dimensioned<symmTensor>
            (
                "source",
                dimless/dimTime,
                symmTensor(0, 2, 0, 1, 1.5, 0)
            )
        )
    );

    Info<< nl
        << "Detailed SolverPerformance<symmTensor>: " << nl
        << "  " << sP << endl;

    Info<< nl
        << "solverPerformanceDict: "
        << mesh.solverPerformanceDict() << endl;


    if (args.found("zip"))
    {
        Info<< nl << "Running some zip/unzip tests" << nl;

        auto cmpts = create<volScalarField>(mesh, 9);
        auto slice = create<volVectorField>(mesh, 3);

        // vector
        {
            setDims(cmpts, U);
            setDims(slice, U);
            rename(cmpts, U.name());
            rename(slice, U.name() + "-slice");

            Info<< "unzip: " << U.name() << nl;

            unzip(U, cmpts[0], cmpts[1], cmpts[2]);

            zip(U, cmpts[0], cmpts[1], cmpts[2]);

            U.rename(IOobject::groupName(U.name(), "rezip"));

            write(cmpts, 3);
            U.write();
        }

        // tensor
        {
            setDims(cmpts, tensf);
            setDims(slice, tensf);
            rename(cmpts, tensf.name());
            rename(slice, tensf.name() + "-slice");

            Info<< "unzip: " << tensf.name() << nl;

            unzip
            (
                tensf,
                cmpts[0], cmpts[1], cmpts[2],
                cmpts[3], cmpts[4], cmpts[5],
                cmpts[6], cmpts[7], cmpts[8]
            );

            unzipRows(tensf, slice[0], slice[1], slice[2]);

            // reszip transposed
            zipCols(tensf, slice[0], slice[1], slice[2]);

            tensf.rename(tensf.name() + "-T");

            write(cmpts, 9);
            write(slice, 3);
            tensf.write();
        }
    }


    // Note: this is definitely not working, but added to help diagnose
    // what is missing if we wish to proceed with this idea
    #ifdef TEST_UINT8_FIELD
    {
        Info<< "uint8 field\n" << endl;
        GeometricField<uint8_t, fvPatchField, volMesh> statefld
        (
            IOobject
            (
                "state",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<uint8_t>(dimless, uint8_t{100})
        );

        Info().beginBlock("stateField")
            << statefld << nl;
        Info().endBlock();
    }
    #endif

    Info<< "\nEnd\n" << nl;

    return 0;
}


// ************************************************************************* //
