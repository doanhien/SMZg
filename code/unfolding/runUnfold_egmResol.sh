#!/bin/bash                                                                                                                                                               

for ((i=1; i <= 100;i+=1)); do
    root -l -b -q unfolding_egmResol.C\(1\,1\,1\,$i\,\"rho_up\"\)
    root -l -b -q unfolding_egmResol.C\(1\,1\,0\,$i\,\"rho_up\"\)
    root -l -b -q unfolding_egmResol.C\(1\,1\,1\,$i\,\"rho_dn\"\)
    root -l -b -q unfolding_egmResol.C\(1\,1\,0\,$i\,\"rho_dn\"\)

    root -l -b -q unfolding_egmResol.C\(1\,1\,1\,$i\,\"phi_up\"\)
    root -l -b -q unfolding_egmResol.C\(1\,1\,0\,$i\,\"phi_up\"\)
    root -l -b -q unfolding_egmResol.C\(1\,1\,1\,$i\,\"phi_dn\"\)
    root -l -b -q unfolding_egmResol.C\(1\,1\,0\,$i\,\"phi_dn\"\)

done
