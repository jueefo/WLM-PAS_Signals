import numpy as np

def GaussianShape(PeakHeight, ModulationWavenumber, CenterWavenumber, FWHM):
#   Gaussian function value at points ModulationWavenumber, when center point is CenterWavenumber for FWHM
#   Gaussian shape approximation for low pressure
    Gaussian = PeakHeight*np.exp((-((ModulationWavenumber-CenterWavenumber)**2)*4*np.log(2))/(FWHM**2))
    
    return Gaussian