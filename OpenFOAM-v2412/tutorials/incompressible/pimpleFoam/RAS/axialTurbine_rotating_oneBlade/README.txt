- case provided by Hakan Nilsson
- demonstrates:
    - meshing in cylindrical coordinates
    - use of cyclicPeriodicAMI
    - setting of non-default matchTolerance (0.01) to get correct number
      of sectors
    - using paraview AngularPeriodicFilter

- touch axialTurbine_rotating_oneBlade.foam
- paraview
    - File -> Load State -> select 'angPer.pvsm'
    - Choose FileName and select the .foam file
