#!/usr/bin/env bash

set -e
set -x
shopt -s dotglob

readonly name="visit-databases"
readonly ownership="VisIt Upstream <kwrobot@kitware.com>"
readonly subtree="databases/readers"
readonly repo="https://gitlab.kitware.com/third-party/visit.git"
readonly tag="for/paraview-20230801-g756b0aa4ee"

readonly paths="
.gitattributes
README.kitware.md

src/databases/ANSYS
src/databases/AUXFile
src/databases/BOV
src/databases/Boxlib3D
src/databases/CEAucd
src/databases/CMAT
src/databases/CTRL
src/databases/Chombo
src/databases/Claw
src/databases/Curve2D
src/databases/DDCMD
src/databases/Dyna3D
src/databases/Enzo
src/databases/ExtrudedVol
src/databases/FLASH
src/databases/Fluent
src/databases/GGCM
src/databases/GTC
src/databases/GULP
src/databases/Gadget
src/databases/H5Nimrod
src/databases/Image
src/databases/LAMMPS
src/databases/Lines
src/databases/M3D
src/databases/M3DC1
src/databases/MFIX
src/databases/MFIXCDF
src/databases/MM5
src/databases/Mili
src/databases/Miranda
src/databases/NASTRAN
src/databases/NETCDF
src/databases/Nek5000
src/databases/OVERFLOW
src/databases/OpenFOAM
src/databases/PATRAN
src/databases/PFLOTRAN
src/databases/PLOT3D
src/databases/Pixie
src/databases/Point3D
src/databases/ProteinDataBank
src/databases/RAW
src/databases/SAMRAI
src/databases/SAR
src/databases/SAS
src/databases/Silo
src/databases/Spheral
src/databases/TFT
src/databases/TSurf
src/databases/Tecplot
src/databases/Tetrad
src/databases/UNIC
src/databases/VASP
src/databases/Velodyne
src/databases/Vs
src/databases/XYZ
src/databases/Xmdv
src/databases/paraDIS
src/databases/paraDIS_tecplot
"

database_cleanup () {
    local reader="$1"
    readonly reader
    shift

    local kind
    case "$reader" in
        Chombo|FLASH|GGCM|M3DC1|MFIXCDF|MFIX|Miranda|NASTRAN|Nek5000|OpenFOAM|paraDIS|Pixie|PLOT3D|ProteinDataBank|Vs)
            kind="reader_options"
            ;;
        BOV|Curve2D|ExtrudedVol|Image|Silo|Tecplot|Xmdv)
            kind="writer"
            ;;
        *)
            kind="reader"
            ;;
    esac
    readonly kind

    rm -v *MDServerPluginInfo.C

    case "$kind" in
        reader)
            rm -vf *PluginInfo.* avt*Options.* avt*Writer.*
            ;;
        reader_options)
            rm -vf *Writer.*
            ;;
        writer)
            ;;
    esac
}

cleanup_LAMMPS () {
    # This just abstracts over the other sub-formats.
    rm -v avtLAMMPSFileFormat.*
}

cleanup_H5Nimrod () {
    # Remove the PDF file describing the format.
    rm -v nimrod.pdf
}

cleanup_PLOT3D () {
    # Documentation tarball.
    rm -v plot3d_html.tar
}

cleanup_VASP () {
    # Unwrapped VASP readers.
    rm -v avtVASPFileFormat.*
    rm -v avtXDATFileFormat.*
}

cleanup_Vs () {
    # Unused files.
    rm -v VsStaggeredField.*
}

cleanup_paraDIS () {
    rm -v RC_cpp_lib/RC_c_lib/gethostname.c
    rm -v RC_cpp_lib/RC_c_lib/gprof-helper.c
    rm -v RC_cpp_lib/RC_c_lib/signals.c
    rm -v RC_cpp_lib/ComboBoxManager.*
    rm -v RC_cpp_lib/RCLibConstants.h
    rm -v RC_cpp_lib/cppFileUtils.h
    rm -v RC_cpp_lib/debugutil-qt.*
    rm -v RC_cpp_lib/test.C
    rm -v paraDIS_tecplot.*
}

cleanup_paraDIS_tecplot () {
    rm -vrf RC_cpp_lib
    rm -vrf RC_c_lib
}

extract_source () {
    git_archive
    pushd "$extractdir/$name-reduced"
    mv --target-directory=. src/databases/*
    find . -name "*.code" -delete
    find . -name "*.xml" -delete
    for reader in *; do
        [ -d "$reader" ] || continue
        pushd "$reader"
        mv -v "CMakeLists.paraview.txt" "CMakeLists.txt"
        database_cleanup "$reader"
        if type "cleanup_$reader" >/dev/null 2>/dev/null; then
            "cleanup_$reader"
        fi
        popd
    done
    popd
}

. "${BASH_SOURCE%/*}/update-common.sh"
