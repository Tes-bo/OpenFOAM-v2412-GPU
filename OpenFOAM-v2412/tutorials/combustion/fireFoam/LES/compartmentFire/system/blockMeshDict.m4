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
// General macros to create 2D/extruded-2D meshes


changecom(//)changequote([,])
define(calc, [esyscmd(perl -e 'printf ($1)')])
define(VCOUNT, 0)
define(vlabel, [[// ]Vertex $1 = VCOUNT define($1, VCOUNT)define([VCOUNT], incr(VCOUNT))])
define(pi, 3.14159265)
define(angle, 45)

scale   1;

/*  User parameter for geometry  */
define(L_square, 0.06)
define(R_burner, 0.0475)
define(H_box, 0.4)
define(L_airBox, 0.4)
define(H_airBox, 0.6)

/*  User parameter for mesh resolution  */
define(dx, 6)
define(dz, 6)
define(dy, 40)
define(d_inlet, 2)//nb of pts between L_square and R_burner
define(d_wall, 17)//nb of pts between R_burner and H_box
define(dx_airBox, 25)
define(dy_airBox, 15)
/********************/

define(X_s, calc(L_square/2))
define(Z_s, calc(L_square/2))

define(angleRad, calc(angle*pi/180))
define(X_b, calc(cos(angleRad)*R_burner))
define(Z_b, calc(sin(angleRad)*R_burner))

define(X_h, calc(H_box/2))
define(Z_h, calc(H_box/2))

vertices
(
// (X Y Z)

( -X_s 0 -Z_s)    //vertex 0
(  X_s 0 -Z_s)    //vertex 1
(  X_s 0  Z_s)    //vertex 2
( -X_s 0  Z_s)    //vertex 3

( -X_b 0 -Z_b)    //vertex 4
(  X_b 0 -Z_b)    //vertex 5
(  X_b 0  Z_b)    //vertex 6
( -X_b 0  Z_b)    //vertex 7

( -X_h 0 -Z_h)    //vertex 8
( -X_b 0 -Z_h)    //vertex 9
(  X_b 0 -Z_h)    //vertex 10
(  X_h 0 -Z_h)    //vertex 11
(  X_h 0 -Z_b)    //vertex 12
(  X_h 0  Z_b)    //vertex 13
(  X_h 0  Z_h)    //vertex 14
(  X_b 0  Z_h)    //vertex 15
( -X_b 0  Z_h)    //vertex 16
( -X_h 0  Z_h)    //vertex 17
( -X_h 0  Z_b)    //vertex 18
( -X_h 0 -Z_b)    //vertex 19



( -X_s H_box -Z_s)    //vertex 20
(  X_s H_box -Z_s)    //vertex 21
(  X_s H_box  Z_s)    //vertex 22
( -X_s H_box  Z_s)    //vertex 23

( -X_b H_box -Z_b)    //vertex 24
(  X_b H_box -Z_b)    //vertex 25
(  X_b H_box  Z_b)    //vertex 26
( -X_b H_box  Z_b)    //vertex 27

( -X_h H_box -Z_h)    //vertex 28
( -X_b H_box -Z_h)    //vertex 29
(  X_b H_box -Z_h)    //vertex 30
(  X_h H_box -Z_h)    //vertex 31
(  X_h H_box -Z_b)    //vertex 32
(  X_h H_box  Z_b)    //vertex 33
(  X_h H_box  Z_h)    //vertex 34
(  X_b H_box  Z_h)    //vertex 35
( -X_b H_box  Z_h)    //vertex 36
( -X_h H_box  Z_h)    //vertex 37
( -X_h H_box  Z_b)    //vertex 38
( -X_h H_box -Z_b)    //vertex 39

( calc(-X_h-L_airBox) 0 -Z_h)    //vertex 40
( calc(-X_h-L_airBox) 0 -Z_b)    //vertex 41
( calc(-X_h-L_airBox) 0  Z_b)    //vertex 42
( calc(-X_h-L_airBox) 0  Z_h)    //vertex 43

( calc(-X_h-L_airBox) H_box -Z_h)     //vertex 44
( calc(-X_h-L_airBox) H_box -Z_b)    //vertex 45
( calc(-X_h-L_airBox) H_box  Z_b)    //vertex 46
( calc(-X_h-L_airBox) H_box  Z_h)    //vertex 47

( -X_h H_airBox -Z_h)    //vertex 48
( -X_h H_airBox -Z_b)    //vertex 49
( -X_h H_airBox  Z_b)    //vertex 50
( -X_h H_airBox  Z_h)    //vertex 51

( calc(-X_h-L_airBox) H_airBox -Z_h)    //vertex 52
( calc(-X_h-L_airBox) H_airBox -Z_b)    //vertex 53
( calc(-X_h-L_airBox) H_airBox  Z_b)    //vertex 54
( calc(-X_h-L_airBox) H_airBox  Z_h)    //vertex 55

);

blocks
(
    hex (0 1 21 20 3 2 22 23) (dx dy dz) simpleGrading (1 1 1) //Block 0
    hex (4 5 25 24 0 1 21 20) (dx dy d_inlet) simpleGrading (1 1 1) //Block 1
    hex (1 5 25 21 2 6 26 22) (d_inlet dy dz) simpleGrading (1 1 1) //Block 2
    hex (3 2 22 23 7 6 26 27) (dx dy d_inlet) simpleGrading (1 1 1) //Block 3
    hex (4 0 20 24 7 3 23 27) (d_inlet dy dz) simpleGrading (1 1 1) //Block 4
    hex (8 9 29 28 19 4 24 39) (d_wall dy d_wall) simpleGrading (1 1 1) //Block 5
    hex (9 10 30 29 4 5 25 24) (dx dy d_wall) simpleGrading (1 1 1) //Block 6
    hex (10 11 31 30 5 12 32 25) (d_wall dy d_wall) simpleGrading (1 1 1) //Block 7
    hex (5 12 32 25 6 13 33 26) (d_wall dy dz) simpleGrading (1 1 1) //Block 8
    hex (6 13 33 26 15 14 34 35) (d_wall dy d_wall) simpleGrading (1 1 1) //Block 9
    hex (7 6 26 27 16 15 35 36) (dx dy d_wall) simpleGrading (1 1 1) //Block 10
    hex (18 7 27 38 17 16 36 37) (d_wall dy d_wall) simpleGrading (1 1 1) //Block 11
    hex (19 4 24 39 18 7 27 38) (d_wall dy dz) simpleGrading (1 1 1) //Block 12

    hex (40 8 28 44 41 19 39 45) (dx_airBox dy d_wall) simpleGrading (0.416667  1 1) //Block 13
    hex (41 19 39 45 42 18 38 46) (dx_airBox dy dz) simpleGrading (0.416667  1 1) //Block 14
    hex (42 18 38 46 43 17 37 47) (dx_airBox dy d_wall) simpleGrading (0.416667  1 1) //Block 15

    hex (44 28 48 52 45 39 49 53) (dx_airBox dy_airBox d_wall) simpleGrading (0.416667 1.72 1) //Block 16
    hex (45 39 49 53 46 38 50 54) (dx_airBox dy_airBox dz) simpleGrading (0.416667  1.72 1) //Block 17
    hex (46 38 50 54 47 37 51 55) (dx_airBox dy_airBox d_wall) simpleGrading (0.416667  1.72 1) //Block 18
);

edges
(
    arc 4 5 origin (0 0 0)
    arc 5 6 origin (0 0 0)
    arc 6 7 origin (0 0 0)
    arc 7 4 origin (0 0 0)

    arc 24 25 origin (0 H_box 0)
    arc 25 26 origin (0 H_box 0)
    arc 26 27 origin (0 H_box 0)
    arc 27 24 origin (0 H_box 0)
);

boundary
(
    inlet
    {
        type patch;
        faces
        (
            (0 1 2 3)
            (1 5 6 2)
            (3 2 6 7)
            (4 0 3 7)
            (4 5 1 0)
        );
    }

    entrainment
    {
        type patch;
        faces
        (
            (40 8 19 41)
            (41 19 18 42)
            (42 18 17 43)

            (40 8 28 44)
            (44 28 48 52)

            (41 40 44 45)
            (42 41 45 46)
            (43 42 46 47)

            (45 44 52 53)
            (46 45 53 54)
            (47 46 54 55)

            (17 43 47 37)
            (37 47 55 51)
        );
    }

    outlet
    {
        type patch;
        faces
        (
            (52 48 49 53)
            (53 49 50 54)
            (54 50 51 55)

            (28 39 49 48)
            (39 38 50 49)
            (38 37 51 50)
        );

    }
);

mergePatchPairs
(
);

// ************************************************************************* //
