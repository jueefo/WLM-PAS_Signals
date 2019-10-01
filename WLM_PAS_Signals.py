#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 09:25:50 2019

@author: juha

Illustration of the frequency distribution of the wavelength (or wavenumber) modulation
    
"""

def main():
    import numpy as np
    import matplotlib.pyplot as plt
    import scipy.signal as signal
    
    import ModulationWN # Module to calculate modulation waveform
    import PAResponse # Photoacoustic response
    import SpectralLineShape # Module to calculate spectral lineshape
    
    PeakHeight = 1  # Spectral intensity
    FWHM = 1E-1 # 2*sqrt(2*ln(2))*c => Full width at half maximum for Gaussian lineshape
    CenterWavenumber = 1278 
    ModFreq = 80 # Modulation frequency (typical 80 Hz)
    MeasTime = 1 # Signal integration time (e.g. 1-5 s)
    Sampling = 50000 # Sampling frequency (e.g. 50000/s)
    Time = np.linspace(0, MeasTime, Sampling*MeasTime) # Sampling times

    ModulationWavenumberE = 0.05*FWHM # Modulation error with in the modulation wavenumber position
    
    CantileverGap = 5 # The gap of the cantilever acts as a high pass filter (e.g. 5 Hz)
    Heating = 10 # The heat transfer ratio acts as low pass filter (e.g. 10 Hz])

#   Modulation limit calculated from the middle of peak. Optimal modulation is about +/- 1.1 x FWHM
    ModulationStartWN = CenterWavenumber-1.1*FWHM
    ModulationEndWN = CenterWavenumber+1.1*FWHM

#  Triangle modulation waveform
    ModulationWavenumber = ModulationWN.Triangle(ModulationStartWN, ModulationEndWN, ModFreq, Sampling, MeasTime)

#   Modulation with positon error
    ModulationWavenumber = ModulationWavenumber + ModulationWavenumberE
    LaserPower = np.ones(ModulationWavenumber.size) #   Laser power is assumed to be constant over modulation    

#   Signal: Laser power x Gaussian function over Time
    GaussianSignal = SpectralLineShape.GaussianShape(PeakHeight, ModulationWavenumber, CenterWavenumber, FWHM)
    Signal = LaserPower*GaussianSignal
    
#   Siganl response with lowpass and highpass  
    PASignal = PAResponse.butter_lowpass_filter(Signal, Heating, Sampling, 1)
    PASignal = PAResponse.butter_highpass_filter(PASignal, CantileverGap, Sampling, 1)

#   Spectra calculations
    Spectrum = np.fft.fft(Signal) # spectrum of raw signal
    PASpectrum = np.fft.fft(PASignal) # PAS spectrum
    freq = np.fft.fftfreq(Time.size)
    freq = freq*Sampling # actual frequency

#   Plots
    plt.subplot(2, 1, 1)
    plt.plot(Time, Signal, 'b-', label='Signal (flat response)')
    plt.xlabel('Time [sec]')
    plt.grid()
    plt.legend()
    
    plt.subplot(2, 1, 2)
    plt.plot(Time, PASignal, 'g-', linewidth=2, label='PA-signal')
    plt.xlabel('Time [sec]')
    plt.grid()
    plt.legend()

    plt.subplots_adjust(hspace=0.35)
    plt.show()
    plt.close()

    PlotFreq = ModFreq*5
    
    plt.subplot(2, 1, 1)
    plt.plot(freq, abs(Spectrum), 'b-', label='Spectrum (flat response)')
    plt.xlim(0, PlotFreq)
    plt.xlabel('Frequency [Hz]')
    plt.grid()
    plt.legend()

    plt.subplot(2, 1, 2)
    plt.plot(freq, abs(PASpectrum), 'g-', label='PA-spectrum')
    plt.xlim(0, PlotFreq)
    plt.xlabel('Frequency [Hz]')
    plt.grid()
    plt.legend()

    plt.subplots_adjust(hspace=0.35)
    plt.show()
    plt.close()
    
#   Plot the frequency response.
#   Get the filter coefficients so we can check its frequency response.
    blow, alow = PAResponse.butter_lowpass(Heating, Sampling, 1)
    bhigh, ahigh = PAResponse.butter_highpass(CantileverGap, Sampling, 1)

    wlow, hlow = signal.freqz(blow, alow, worN=100001)
    whigh, hhigh = signal.freqz(bhigh, ahigh, worN=100001)

    plt.subplot(2, 1, 1)
    plt.plot(0.5*Sampling*wlow/np.pi, abs(hlow), 'b-', label='Heating effect')
    plt.plot(Heating, 0.5*np.sqrt(2), 'ko')
    plt.axvline(Heating, color='k')
    plt.plot(0.5*Sampling*whigh/np.pi, np.abs(hhigh), 'g-', label='Cantilever gap effect')
    plt.plot(CantileverGap, 0.5*np.sqrt(2), 'ko')
    plt.axvline(CantileverGap, color='k')
    plt.xlim(0, 4*ModFreq)
    plt.xlabel('Frequency [Hz]')
    plt.grid()
    plt.legend()
    
    plt.subplot(2, 1, 2)
    plt.plot(0.5*Sampling*whigh/np.pi, np.abs(hhigh)*np.abs(hlow), 'r-', label='Total response below resonance')
    plt.xlim(0, PlotFreq)
    plt.xlabel('Frequency [Hz]')
    plt.grid()
    plt.legend()
    
    plt.subplots_adjust(hspace=0.35)
    plt.show()
    plt.close()

if __name__=='__main__':
    main()