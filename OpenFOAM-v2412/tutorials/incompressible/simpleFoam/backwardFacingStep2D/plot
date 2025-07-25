#!/bin/bash
cd "${0%/*}" || exit                                # Run from this directory
. ${WM_PROJECT_DIR:?}/bin/tools/RunFunctions        # Tutorial run functions
#------------------------------------------------------------------------------

plot_U() {

    sampleDir="$1"
    Uref="$2"
    href="0.0127"

    image="backwardStep2D_U.png"

    gnuplot<<EOF1
    set terminal pngcairo font "helvetica,20" size 1000, 1000
    set grid
    set key top right
    set xrange [-0.4:1.2]
    set yrange [0:3]
    set xlabel "U/U_0"
    set ylabel "y/h"
    set output "$image"
    set format x "%.1f"
    set format y "%.1f"

    # OpenFOAM
        samples="$sampleDir"
        Uref="$Uref"
        href="$href"
        fields="p_U_turbulenceProperties:devReff.xy"

    plot \
        samples."/x_by_h_01_".fields u (\$3/Uref):(\$1/href) w l lw 2 lc rgb "red" t "x/h = 1", \
        samples."/x_by_h_04_".fields u (\$3/Uref):(\$1/href) w l lw 2 lc rgb "green" t "x/h = 4", \
        samples."/x_by_h_06_".fields u (\$3/Uref):(\$1/href) w l lw 2 lc rgb "blue" t "x/h = 6", \
        samples."/x_by_h_10_".fields u (\$3/Uref):(\$1/href) w l lw 2 lc rgb "black" t "x/h = 10"
EOF1
}


plot_tau() {

    timeDir="$1"
    Uref="$2"
    href="0.0127"

    echo "    # ccx tau_xx tau_yy tau_zz" > tauw.dat
    foamDictionary -entry boundaryField/lowerWall/value -value "$timeDir"/Cx | \
        sed -n '/(/,/)/p' | sed -e 's/[()]//g;/^\s*$/d' > Cx.$$
    foamDictionary -entry boundaryField/lowerWall/value -value "$timeDir"/wallShearStress | \
        sed -n '/(/,/)/p' | sed -e 's/[()]//g;/^\s*$/d' > tau.$$
    paste -d ' ' Cx.$$ tau.$$ >> tauw.dat
    rm -f Cx.$$ tau.$$

    image="backwardStep2D_tau.png"

    gnuplot<<EOF2
    set terminal pngcairo font "helvetica,20" size 1000, 1000
    set grid
    set key bottom right
    set xrange [-5:30]
    set yrange [-0.002:0.004]
    set xlabel "x/h"
    set ylabel "C_f"
    set output "$image"

    # OpenFOAM
        Uref="$Uref"
        href="$href"

    plot \
        "tauw.dat" u (\$1/href):(-\$2/(0.5*Uref*Uref)) t "simpleFoam" w l lw 2 lc rgb "black"
EOF2
}


#------------------------------------------------------------------------------

# Requires gnuplot
command -v gnuplot >/dev/null || {
    echo "FOAM FATAL ERROR: gnuplot not found - skipping graph creation" 1>&2
    exit 1
}

# Requires awk
command -v awk >/dev/null || {
    echo "FOAM FATAL ERROR: awk not found - skipping graph creation" 1>&2
    exit 1
}

# Check "results" directory
timeDir=$(foamListTimes -latestTime)
sampleDir=postProcessing/sample/"$timeDir"

[ -d "$sampleDir" ] || {
    echo "FOAM FATAL ERROR: sample/"$timeDir" directory not found - skipping graph creation" 1>&2
    exit 1
}


#------------------------------------------------------------------------------

Uref=$(awk '{print $2}' $sampleDir/Uref_p_U_turbulenceProperties:devReff.xy)

echo ""
echo "# Plots the U profiles"
echo ""

plot_U "$sampleDir" "$Uref"

echo ""
echo "# Plots the wall-shear stress profiles"
echo ""

plot_tau "$timeDir" "$Uref"


#------------------------------------------------------------------------------
