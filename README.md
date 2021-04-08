# TES_signal_sim

This project includes signal processing simulation programs(most in MATLAB)，serving to ALI TES-READOUT system in China.

## Poly-phase Filter Bank (PFB)
- pfb_fir.m:     PFB implemention with maximumly decimation by M:1
- pfb_psd_mean.m:  simple power computation of each channel, optional integration. 
- os2_pfb2v1_fir.m:  PFB implemention with oversampling by M/2:1 decimation
- os2_pfb4v3_fir.m:  PFB implemention with oversampling by 4M/3:1 decimation
- os2_pfb8v5_fir.m:  PFB implemention with oversampling by 8M/5:1 decimation
