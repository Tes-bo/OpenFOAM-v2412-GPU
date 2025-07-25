testcase to demonstrate surface-curvature-based refinement.
- starts from a single cell
- has uniform surface refinement 0
- but has curvature detection to refine up to 10
- also uses limitRegions to not refine half the domain (in x direction)
- note that there is still refinement bleeding due to the 2:1 limitation
