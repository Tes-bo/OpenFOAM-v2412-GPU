set multiplot layout 2,2
set grid

set xl 'Optimization cycle'

set yl 'C_d/C_{d,init}'
p 'optimisation/objective/0/dragas1' u 1:3 w lp notitle

set yl 'C_l'
set yr [0.05266:0.05366]
p 'optimisation/objective/0/liftlift' u 1:2 w lp notitle, 0.05340147181439196 t 'upper bound', 0.05320147181439196 t 'lower bound'

set autoscale y
set yl '(V-Vinit)/Vinit'
p 'optimisation/objective/0/volvol' u 1:2 w lp notitle, -0.15 t 'lower bound'
