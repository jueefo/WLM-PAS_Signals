def Triangle(ModulationStart, ModulationEnd, Freq, Sampling, MeasTime):
#   Triange modulation of lenght Sampling*Time
    import numpy as np
    
    Time = np.linspace(0, MeasTime, Sampling*MeasTime)
#   ModulationWavenumberTemp one modulation period
    ModulationWavenumberTemp = np.append(np.linspace(ModulationStart, (ModulationEnd-(ModulationEnd-ModulationStart)/Sampling*2), int(Sampling/Freq/2)),
                                         np.linspace(ModulationEnd, (ModulationStart+(ModulationEnd-ModulationStart)/Sampling*2), int(Sampling/Freq/2)))
    ModulationWavenumber=ModulationWavenumberTemp

#   Modulation periods added to whole measurement time
    for Findex in range(MeasTime*Freq):
        ModulationWavenumber = np.append(ModulationWavenumberTemp, ModulationWavenumber)

#   Cut the modulation lenght to measurement time
    ModulationWavenumber = ModulationWavenumber[0:Time.size]
    return (ModulationWavenumber)