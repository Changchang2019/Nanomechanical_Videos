# Nanomechanical_Videos

**GNU_Radio folder**

System requirements: GNU Radio Companion 3.10.1.1 (Python 3.10.12)

Typical install time on a "normal" desktop computer: immediate

Expected run time for demo on a "normal" desktop computer: a few seconds

*dvbs_tx.grc* GNU Radio flowgraph used to transmit a DVB-S2 video stream.

*dvbs_rx.grc* GNU Radio flowgraph used to receive a DVB-S2 video stream.

*bit_stream_tx.grc* GNU Radio flowgraph used to transmit a bitstream to measure the Bit Error Ratio (BER).

*bit_stream_rx.grc* GNU Radio flowgraph used to receive a bitstream to measure BER.

*labvid1M_qpsk_5to6aOFF.ts* TS file used in the video transmission experiment (Supplementary Video 1).

*random_data_byte.txt* Random data used in BER measurements.



**Matlab folder**

System requirements: MATLAB R2022a

Expected run time for demo on a "normal" desktop computer: from a few minutes to a few hours, depending on the number of points in the drive waveform and the number of drive frequencies.

*Calculate_L_epsilon_BER.m* is the main MATLAB code to calculate L, epsilon, and BER. The code is annotated.
It calls the functions *decim_joel.m*, *eq_motion_JM_QPSK_v2B.m*, *rotate_IQ.m*, and *EVM_joel_v6.m*.
