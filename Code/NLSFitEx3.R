library(minpack.lm)
library(graphics)

Data <- read.csv(file="../Data/ClimateData.csv", header=TRUE)

Years <- Data[, 1]

Temperatures = Data[, 4]

x = seq(length(Years)) #~ pseudo-timepoints to feed to SinMod 

# define objective function: returns the array to be minimized
SinMod = function(x, data){
    Shift = params['Shift'].value
    Amp = params['Amp'].value
    Len = params['Len'].value
    Phase = params['Phase'].value

    return (Shift + Amp * sc.sin((2 * sc.pi * x / Len) + Phase) }

# create a set of Parameters
params = Parameters()
params.add('Amp',   value= 5,  min=0)
params.add('Shift', value= 10.0, min=5, max=15) #you can add bounds
params.add('Len', value= 12.0) # Sensitive to this parameter!
params.add('Phase', value= 1.0)

# do fit (you can try different algorithms, such as lestsq, nelder, etc)

result = minimize(SinMod, params, args=(x, Temperatures),method="leastsq")

# calculate final result
final = Temperatures + result.residual

# write error report
report_fit(result.params)

# Plot results
try: # poor use of try!
    pl.close('all')
    pl.ion()
    fig = pl.plot(x, Temperatures, '+k--')
    pl.plot(x, final, 'r')
    pl.title('lmfit to Climatic fluctuations')
    pl.xlabel('Time')
    pl.ylabel('Temperature ($^\circ$C)')
except:
    pass
