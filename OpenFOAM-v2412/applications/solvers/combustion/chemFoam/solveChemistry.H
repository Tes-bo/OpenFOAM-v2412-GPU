dtChem = chemistry.solve(runTime.deltaTValue());
scalar Qdot = chemistry.Qdot()()[0]/rho[0];
integratedHeat += Qdot*runTime.deltaTValue();
