- blockMesh
- decomposePar
- in processor1:
    - swap points 2 (1 1 0) and 6 (1 1 2) in polyMesh/points
    - swap indices 2 and 6 in polyMesh/faces

so now we have the same mesh but different edge orientation on processor1
