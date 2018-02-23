from parcels import *
import numpy as np
from datetime import timedelta

lon = np.arange(0, 13.6-13.6/100, 13.6/100)
lat = np.arange(0, 13.6-13.6/100, 13.6/100)
#grid = RectilinearZGrid(long, lat, mesh = 'flat')  # from_data creates its own grid (like this one) I think?


field_data = np.load('/media/alexander/DATA/Ubuntu/Miniproject/Parcels Learning/fields/offline_grid_H_low.npy')
fieldset = FieldSet.from_data(data={'U': field_data[:, :, 0], 'V': field_data[:, :, 1]}, dimensions={'lon': lon, 'lat': lat}, transpose=False, mesh='flat')

# pset = ParticleSet.from_list(fieldset=fieldset,
#                              pclass=JITParticle,
#                              lon=np.hstack((np.arange(2.6, 2.8, 0.01), np.arange(10.81, 11.01, 0.01))).tolist(),
#                              lat=np.full((40, 1), 0.5).tolist())
pset = ParticleSet.from_line(fieldset=fieldset,
                             pclass=JITParticle,
                             start=[2.8, 6.8],
                             finish=[10.8, 6.8],
                             size=20)

def DeleteParticle(particle, fieldset, time, dt):
    particle.delete()

pset.execute(AdvectionRK4,                 # the kernel (which defines how particles move)
             runtime=timedelta(seconds=1),    # the total length of the run
             dt=timedelta(seconds=-0.1),      # the timestep of the kernel
             recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})
             #moviedt=timedelta(seconds=0.1),
             #movie_background_field=fieldset.U)
             #output_file=pset.ParticleFile(name="../parcels_outputs/H/testrun_backwards", outputdt=timedelta(seconds=0.1)))
pset.show(savefile='/media/alexander/DATA/Ubuntu/Miniproject/Parcels Learning/outputs/H/animations/particles', field=fieldset.U, vmin=-1.1, vmax=1.1, domain=[8.8, 10.8, 2.8, 10.8])

#plotTrajectoriesFile('../outputs/H/testrun_backwards.nc')

print()