- 2x2x1 mesh
- remove one cell, exposing two faces
- move exposed faces into two patches
- decompose onto 3
- run surfaceMeshExtract -featureAngle 180
- should also mark points on the processor that has no
  faces but is coupled
