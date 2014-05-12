  
 Enclosed please find my simulation result and code.

* I have implemented optimal quadratic non-coherent demodulation for CP-BFSK signal.

* I have extensively tested my code on different frequency deviation, bit period etc.

* For this simulation, I have used the parameters used in plot_ber_cpfsk.m

* For each SINR, I have used 100 runs. That takes a lot of time to run. So, I have set num_run=10 in plot_ber_cpfsk.m

* plot_ber_cpfsk.m is the main routine that calls other signal generation, demodulation and filtering routines.

* The code can be tested by typing


>> [AV_NF_BER AV_SMF_BER1 AV_SMF_BER2 AV_MF_BER]   = plot_ber_cpfsk()


