/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2412                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    `format'      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scale   0.001;

changecom(//)changequote([,]) dnl>
define(calc, [esyscmd(perl -e 'printf ($1)')]) dnl>
define(VCOUNT, 0)
define(vlabel, [[// ]Vertex $1 = VCOUNT define($1, VCOUNT)define([VCOUNT], incr(VCOUNT))])

define(Dt, 220)
define(DdDt, 0.33)
define(HtDt, 3)
define(PI, 3.14159265)

define(Dd, calc(Dt * DdDt))
define(Ht, calc(HtDt * Dt))
define(psB, calc(Ht/10.0))
define(psH, calc(psB * 1.5))

define(Hl, calc(Ht-psH-psB))
define(Hl2, calc((Hl/2.0)+psB))
define(Hl3, calc(Hl+psB))

define(eDd, calc(Dt/100))

define(Rd, calc(Dd/2))
define(Rd2, calc(Dd/4))
define(Rd3, calc(Rd + eDd))

define(Rt, calc(Dt/2))
define(qRt, calc(Rt * 0.75))

define(cx, calc(Rd2*cos((PI/180)*45)))
define(cy, calc(Rd2*sin((PI/180)*45)))

define(ccx, calc(Rd*cos((PI/180)*45)))
define(ccy, calc(Rd*sin((PI/180)*45)))

define(ecx, calc(Rd3*cos((PI/180)*45)))
define(ecy, calc(Rd3*sin((PI/180)*45)))

define(tcx, calc(Rt*cos((PI/180)*45)))
define(tcy, calc(Rt*sin((PI/180)*45)))

define(qcx, calc(Rt*0.75*cos((PI/180)*45)))
define(qcy, calc(Rt*0.75*sin((PI/180)*45)))

define(NPS, 7)   //square section
define(NPS1, 25)
define(NPS2, 5) //vertical in the lower cylinder
define(NPS3, 10) //vertical in the upper cylindre
define(NPD, 8)   //square section to perimeter
define(NPC, 1)   //center
define(NPE, 8)  //external


vertices
(
    ( -cx -cy  0.0)
    (  cx -cy  0.0)
    (  cx  cy  0.0)
    ( -cx  cy  0.0)

    ( -cx -cy  psB)
    (  cx -cy  psB)
    (  cx  cy  psB)
    ( -cx  cy  psB)

    ( ccx -ccy  0.0)
    ( ccx  ccy  0.0)
    ( ccx -ccy  psB)
    ( ccx  ccy  psB)

    ( -ccx -ccy  0.0)
    ( -ccx  ccy  0.0)
    ( -ccx -ccy  psB)
    ( -ccx  ccy  psB)

    ( ecx -ecy  0.0)   //16
    ( ecx  ecy  0.0)
    ( ecx -ecy  psB)
    ( ecx  ecy  psB)

    ( -ecx -ecy  0.0)
    ( -ecx  ecy  0.0)
    ( -ecx -ecy  psB)
    ( -ecx  ecy  psB)

    ( qcx -qcy  0.0) //24
    ( qcx  qcy  0.0)
    ( qcx -qcy  psB)
    ( qcx  qcy  psB)

    ( -qcx -qcy  0.0)
    ( -qcx  qcy  0.0)
    ( -qcx -qcy  psB)
    ( -qcx  qcy  psB)

    ( tcx -tcy  0.0) //32
    ( tcx  tcy  0.0)
    ( tcx -tcy  psB)
    ( tcx  tcy  psB)

    ( -tcx -tcy  0.0)
    ( -tcx  tcy  0.0)
    ( -tcx -tcy  psB)
    ( -tcx  tcy  psB)

    ( qcx -qcy  Hl2) //40
    ( qcx  qcy  Hl2)
    ( -qcx -qcy  Hl2)
    ( -qcx  qcy  Hl2)

    ( tcx -tcy  Hl2) //44
    ( tcx  tcy  Hl2)
    ( -tcx -tcy  Hl2)
    ( -tcx  tcy  Hl2)

    ( ecx -ecy  Hl2) //48
    ( ecx  ecy  Hl2)
    ( -ecx -ecy  Hl2)
    ( -ecx  ecy  Hl2)

    ( -cx -cy  Hl2) //52
    (  cx -cy  Hl2)
    (  cx  cy  Hl2)
    ( -cx  cy  Hl2)

    ( ccx -ccy  Hl2) //56
    ( ccx  ccy  Hl2)
    ( -ccx -ccy  Hl2)
    ( -ccx  ccy  Hl2)

    ( -cx -cy  Hl3) //60
    (  cx -cy  Hl3)
    (  cx  cy  Hl3)
    ( -cx  cy  Hl3)

    ( ccx -ccy  Hl3) //64
    ( ccx  ccy  Hl3)
    ( -ccx -ccy  Hl3)
    ( -ccx  ccy  Hl3)

    ( qcx -qcy  Hl3) //68
    ( qcx  qcy  Hl3)
    ( -qcx -qcy  Hl3)
    ( -qcx  qcy  Hl3)

    ( ecx -ecy  Hl3) //72
    ( ecx  ecy  Hl3)
    ( -ecx -ecy  Hl3)
    ( -ecx  ecy  Hl3)

    ( tcx -tcy  Hl3) //76
    ( tcx  tcy  Hl3)
    ( -tcx -tcy  Hl3)
    ( -tcx  tcy  Hl3)

    ( -cx -cy  Ht) //80
    (  cx -cy  Ht)
    (  cx  cy  Ht)
    ( -cx  cy  Ht)

    ( ccx -ccy  Ht) //84
    ( ccx  ccy  Ht)
    ( -ccx -ccy  Ht)
    ( -ccx  ccy  Ht)

    ( qcx -qcy  Ht) //88
    ( qcx  qcy  Ht)
    ( -qcx -qcy  Ht)
    ( -qcx  qcy  Ht)

    ( ecx -ecy  Ht) //92
    ( ecx  ecy  Ht)
    ( -ecx -ecy  Ht)
    ( -ecx  ecy  Ht)

    ( tcx -tcy  Ht) //96
    ( tcx  tcy  Ht)
    ( -tcx -tcy  Ht)
    ( -tcx  tcy  Ht)
);

blocks
(
    hex
    (
        0 1 2 3
        4 5 6 7
    )
    (NPS NPS NPS2)
    simpleGrading (1 1 1)

    //quartier est
    hex
    (
        1 8 9 2
        5 10 11 6
    )
    (NPD NPS NPS2)
    simpleGrading (0.9 1 1)

    //quartier ouest
    hex
    (
        3 13 12 0
        7 15 14 4
    )
    (NPD NPS NPS2)
    simpleGrading (0.9 1 1)

    //quartier sud
    hex
    (
        0 12 8 1
        4 14 10 5
    )
    (NPD NPS NPS2)
    simpleGrading (0.9 1 1)

    //quartier nord
    hex
    (
        2 9  13 3
        6 11 15 7
    )
    (NPD NPS NPS2)
    simpleGrading (0.9 1 1)

// **************************************

// Fabrication couronne inférieure 4 et 5

    //ceinture est
    hex
    (
        8 16 17 9
        10 18 19 11
    )
    (NPC NPS NPS2)
    simpleGrading (1 1 1)

    //ceinture ouest
    hex
    (
        13 21 20 12
        15 23 22 14
    )
    (NPC NPS NPS2)
    simpleGrading (1 1 1)

    //ceinture sud
    hex
    (
        12 20 16 8
        14 22 18 10
    )
    (NPC NPS NPS2)
    simpleGrading (1 1 1)

    //ceinture nord
    hex
    (
        9  17 21 13
        11 19 23 15
    )
    (NPC NPS NPS2)
    simpleGrading (1 1 1)


// Couronne inférieure externe
    //mi-couronne est
    hex
    (
        16 24 25 17
        18 26 27 19
    )
    (NPE NPS NPS2)
    simpleGrading (1.1 1 1)

    //mi-couronne ouest
    hex
    (
        21 29 28 20
        23 31 30 22
    )
    (NPE NPS NPS2)
    simpleGrading (1.1 1 1)

    //mi-couronne sud
    hex
    (
        20 28 24 16
        22 30 26 18
    )
    (NPE NPS NPS2)
    simpleGrading (1.1 1 1)

    //mi-couronne nord
    hex
    (
        17 25 29 21
        19 27 31 23
    )
    (NPE NPS NPS2)
    simpleGrading (1.1 1 1)

    //mi-couronne est2
    hex
    (
        24 32 33 25
        26 34 35 27
    )
    (NPE NPS NPS2)
    simpleGrading (0.9 1 1)

    //mi-couronne ouest2
    hex
    (
        29 37 36 28
        31 39 38 30
    )
    (NPE NPS NPS2)
    simpleGrading (0.9 1 1)

    //mi-couronne sud2
    hex
    (
        28 36 32 24
        30 38 34 26
    )
    (NPE NPS NPS2)
    simpleGrading (0.9 1 1)

    //mi-couronne nord2
    hex
    (
        25 33 37 29
        27 35 39 31
    )
    (NPE NPS NPS2)
    simpleGrading (0.9 1 1)

// LongBas

    //est
    hex
    (
        18 26 27 19
        48 40 41 49
    )
    (NPE NPS NPS1)
    simpleGrading (1.1 1 1)

    //ouest
    hex
    (
        23 31 30 22
        51 43 42 50
    )
    (NPE NPS NPS1)
    simpleGrading (1.1 1 1)

    //sud
    hex
    (
        22 30 26 18
        50 42 40 48
    )
    (NPE NPS NPS1)
    simpleGrading (1.1 1 1)

    //Nord
    hex
    (
        19 27 31 23
        49 41 43 51
    )
    (NPE NPS NPS1)
    simpleGrading (1.1 1 1)

    //square
    hex
    (
        4  5  6  7
        52 53 54 55
    )
    (NPS NPS NPS1)
    simpleGrading (1 1 1)

    // est-in
    hex
    (
        5  10 11 6
        53 56 57 54
    )
    (NPD NPS NPS1)
    simpleGrading (0.9 1 1)

    // ouest-in
    hex
    (
        7  15 14 4
        55 59 58 52
    )
    (NPD NPS NPS1)
    simpleGrading (0.9 1 1)

    // sud-in
    hex
    (
        4  14 10 5
        52 58 56 53
    )
    (NPD NPS NPS1)
    simpleGrading (0.9 1 1)

    // nord-in
    hex
    (
        6  11 15 7
        54 57 59 55
    )
    (NPD NPS NPS1)
    simpleGrading (0.9 1 1)

// Couronne exterieure longue

    // est
    hex
    (
        26 34 35 27
        40 44 45 41
    )
    (NPE NPS NPS1)
    simpleGrading (0.9 1 1)

    // ouest
    hex
    (
        31 39 38 30
        43 47 46 42
    )
    (NPE NPS NPS1)
    simpleGrading (0.9 1 1)

    // sud
    hex
    (
        30 38 34 26
        42 46 44 40
    )
    (NPE NPS NPS1)
    simpleGrading (0.9 1 1)

    // nord
    hex
    (
        27 35 39 31
        41 45 47 43
    )
    (NPE NPS NPS1)
    simpleGrading (0.9 1 1)

// longHaut

    //square
    hex
    (
        52 53 54 55
        60 61 62 63
    )
    (NPS NPS NPS1)
    simpleGrading (1 1 1)

    // est-in
    hex
    (
        53 56 57 54
        61 64 65 62
    )
    (NPD NPS NPS1)
    simpleGrading (0.9 1 1)

    // ouest-in
    hex
    (
        55 59 58 52
        63 67 66 60
    )
    (NPD NPS NPS1)
    simpleGrading (0.9 1 1)

    // sud-in
    hex
    (
        52 58 56 53
        60 66 64 61
    )
    (NPD NPS NPS1)
    simpleGrading (0.9 1 1)

    // nord-in
    hex
    (
        54 57 59 55
        62 65 67 63
    )
    (NPD NPS NPS1)
    simpleGrading (0.9 1 1)

    //est
    hex
    (
        48 40 41 49
        72 68 69 73
    )
    (NPE NPS NPS1)
    simpleGrading (1.1 1 1)

    //ouest
    hex
    (
        51 43 42 50
        75 71 70 74
    )
    (NPE NPS NPS1)
    simpleGrading (1.1 1 1)

    //sud
    hex
    (
        50 42 40 48
        74 70 68 72
    )
    (NPE NPS NPS1)
    simpleGrading (1.1 1 1)

    //Nord
    hex
    (
        49 41 43 51
        73 69 71 75
    )
    (NPE NPS NPS1)
    simpleGrading (1.1 1 1)

// couronne externe haute

    // est
    hex
    (
        40 44 45 41
        68 76 77 69
    )
    (NPE NPS NPS1)
    simpleGrading (0.9 1 1)

    // ouest
    hex
    (
        43 47 46 42
        71 79 78 70
    )
    (NPE NPS NPS1)
    simpleGrading (0.9 1 1)

    // sud
    hex
    (
        42 46 44 40
        70 78 76 68
    )
    (NPE NPS NPS1)
    simpleGrading (0.9 1 1)

    // nord
    hex
    (
        41 45 47 43
        69 77 79 71
    )
    (NPE NPS NPS1)
    simpleGrading (0.9 1 1)

// Bloc supérieur

    //square
    hex
    (
        60 61 62 63
        80 81 82 83
    )
    (NPS NPS NPS3)
    simpleGrading (1 1 1)

    // est-in
    hex
    (
        61 64 65 62
        81 84 85 82
    )
    (NPD NPS NPS3)
    simpleGrading (0.9 1 1)

    // ouest-in
    hex
    (
        63 67 66 60
        83 87 86 80
    )
    (NPD NPS NPS3)
    simpleGrading (0.9 1 1)

    // sud-in
    hex
    (
        60 66 64 61
        80 86 84 81
    )
    (NPD NPS NPS3)
    simpleGrading (0.9 1 1)

    // nord-in
    hex
    (
        62 65 67 63
        82 85 87 83
    )
    (NPD NPS NPS3)
    simpleGrading (0.9 1 1)

// Fabrication couronne supérieure 6 et 7

    //ceinture est
    hex
    (
        64 72 73 65
        84 92 93 85
    )
    (NPC NPS NPS3)
    simpleGrading (1 1 1)

    //ceinture ouest
    hex
    (
        67 75 74 66
        87 95 94 86
    )
    (NPC NPS NPS3)
    simpleGrading (1 1 1)

    //ceinture sud
    hex
    (
        66 74 72 64
        86 94 92 84
    )
    (NPC NPS NPS3)
    simpleGrading (1 1 1)

    //ceinture nord
    hex
    (
        65 73 75 67
        85 93 95 87
    )
    (NPC NPS NPS3)
    simpleGrading (1 1 1)

    //est
    hex
    (
        72 68 69 73
        92 88 89 93
    )
    (NPE NPS NPS3)
    simpleGrading (1.1 1 1)

    //ouest
    hex
    (
        75 71 70 74
        95 91 90 94
    )
    (NPE NPS NPS3)
    simpleGrading (1.1 1 1)

    //sud
    hex
    (
        74 70 68 72
        94 90 88 92
    )
    (NPE NPS NPS3)
    simpleGrading (1.1 1 1)

    //Nord
    hex
    (
        73 69 71 75
        93 89 91 95
    )
    (NPE NPS NPS3)
    simpleGrading (1.1 1 1)

// couronne externe

    // est
    hex
    (
        68 76 77 69
        88 96 97 89
    )
    (NPE NPS NPS3)
    simpleGrading (0.9 1 1)

    // ouest
    hex
    (
        71 79 78 70
        91 99 98 90
    )
    (NPE NPS NPS3)
    simpleGrading (0.9 1 1)

    // sud
    hex
    (
        70 78 76 68
        90 98 96 88
    )
    (NPE NPS NPS3)
    simpleGrading (0.9 1 1)

    // nord
    hex
    (
        69 77 79 71
        89 97 99 91
    )
    (NPE NPS NPS3)
    simpleGrading (0.9 1 1)
);

edges
(
    arc 8 9 (Rd 0.0 0.0)
    arc 10 11 (Rd 0.0 psB)

    arc 13 12 (-Rd 0.0 0.0)
    arc 15 14 (-Rd 0.0 psB)

    arc 12 8  (0.0 -Rd 0.0)
    arc 14 10 (0.0 -Rd psB)

    arc 9  13 (0.0 Rd 0.0)
    arc 11 15 (0.0 Rd psB)

    arc 16 17 (Rd3 0.0 0.0)
    arc 18 19 (Rd3 0.0 psB)

    arc 21 20 (-Rd3 0.0 0.0)
    arc 23 22 (-Rd3 0.0 psB)

    arc 20 16  (0.0 -Rd3 0.0)
    arc 22 18 (0.0 -Rd3 psB)

    arc 17 21 (0.0 Rd3 0.0)
    arc 19 23 (0.0 Rd3 psB)

    arc 24 25 (qRt 0.0 0.0)
    arc 26 27 (qRt 0.0 psB)

    arc 28 29 (-qRt 0.0 0.0)
    arc 30 31 (-qRt 0.0 psB)

    arc 28 24  (0.0 -qRt 0.0)
    arc 30 26 (0.0 -qRt psB)

    arc 25 29 (0.0 qRt 0.0)
    arc 27 31 (0.0 qRt psB)

    arc 32 33 (Rt 0.0 0.0)
    arc 34 35 (Rt 0.0 psB)

    arc 37 36 (-Rt 0.0 0.0)
    arc 39 38 (-Rt 0.0 psB)

    arc 36 32  (0.0 -Rt 0.0)
    arc 38 34 (0.0 -Rt psB)

    arc 33 37 (0.0 Rt 0.0)
    arc 35 39 (0.0 Rt psB)

    arc 48 49 (Rd3 0.0 Hl2)
    arc 40 41 (qRt 0.0 Hl2)

    arc 51 50 (-Rd3 0.0 Hl2)
    arc 43 42 (-qRt 0.0 Hl2)

    arc 42 40 (0.0 -qRt Hl2)
    arc 50 48 (0.0 -Rd3 Hl2)

    arc 41 43 (0.0 qRt Hl2)
    arc 49 51 (0.0 Rd3 Hl2)

    arc 56 57 (Rd 0.0 Hl2)
    arc 58 59 (-Rd 0.0 Hl2)

    arc 58 56 (0.0 -Rd Hl2)
    arc 57 59 (0.0 Rd Hl2)

    arc 44 45 (Rt 0.0 Hl2)
    arc 46 47 (-Rt 0.0 Hl2)

    arc 46 44 (0.0 -Rt Hl2)
    arc 45 47 (0.0 Rt Hl2)

    arc 64 65 (Rd 0.0 Hl3)
    arc 67 66 (-Rd 0.0 Hl3)

    arc 64 66 (0.0 -Rd Hl3)
    arc 65 67 (0.0 Rd Hl3)

    arc 72 73 (Rd3 0.0 Hl3)
    arc 68 69 (qRt 0.0 Hl3)

    arc 70 71 (-qRt 0.0 Hl3)
    arc 74 75 (-Rd3 0.0 Hl3)

    arc 72 74 (0.0 -Rd3 Hl3)
    arc 68 70 (0.0 -qRt Hl3)

    arc 73 75 (0.0 Rd3 Hl3)
    arc 69 71 (0.0 qRt Hl3)

    arc 76 77 (Rt 0.0 Hl3)
    arc 78 79 (-Rt 0.0 Hl3)

    arc 76 78 (0.0 -Rt Hl3)
    arc 77 79 (0.0 Rt Hl3)

    arc 84 85 (Rd 0.0 Ht)
    arc 86 87 (-Rd 0.0 Ht)

    arc 84 86 (0.0 -Rd Ht)
    arc 85 87 (0.0 Rd Ht)

    arc 92 93 (Rd3 0.0 Ht)
    arc 94 95 (-Rd3 0.0 Ht)

    arc 94 92 (0.0 -Rd3 Ht)
    arc 93 95 (0.0 Rd3 Ht)

    arc 88 89 (qRt 0.0 Ht)
    arc 90 91 (-qRt 0.0 Ht)

    arc 88 90 (0.0 -qRt Ht)
    arc 89 91 (0.0 qRt Ht)

    arc 96 97 (Rt 0.0 Ht)
    arc 98 99 (-Rt 0.0 Ht)

    arc 98 96 (0.0 -Rt Ht)
    arc 97 99 (0.0 Rt Ht)
);


defaultPatch
{
    name    walls;
    type    wall;
}

boundary
(
    base
    {
        type wall;
        faces
        (
            (0 1 2 3)
            (1 8 9 2)
            (2 9 13 3)
            (3 13 12 0)
            (0 12 8 1)
        );
    }

    outlet
    {
        type patch;
        faces
        (
            (80 81 82 83)
            (81 84 85 82)
            (83 87 86 80)
            (80 86 84 81)
            (82 85 87 83)

            (84 92 93 85)
            (87 95 94 86)
            (86 94 92 84)
            (85 93 95 87)

            (92 88 89 93)
            (95 91 90 94)
            (94 90 88 92)
            (93 89 91 95)

            (88 96 97 89)
            (91 99 98 90)
            (90 98 96 88)
            (89 97 99 91)
        );
    }
);

// ************************************************************************* //
