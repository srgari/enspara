#Step 2 Fitting - Loads the assigments matrix. 
from enspara import msm
from enspara import ra

assigs = ra.load('fs-khybrid-clusters0020-assignments.h5')

#Step 3 implied time scale - determine the lagtime for your MSM.	
import numpy as np
from functools import partial

# make 20 different lag times (integers) evenly spaced between 10 and 750
lag_times = np.linspace(10, 750, num=20).astype(int)

implied_timescales = []
for time in lag_times:
    m = msm.MSM(
        lag_time=time,
        method=msm.builders.transpose)
    m.fit(assigs)

    implied_timescales.append(
        -time / np.log(msm.eigenspectrum(m.tprobs_, n_eigs=4)[0][1:3])
    )


#Step 4 Plotting (line missing in the tutorial). 
import matplotlib.pyplot as plt

implied_timescales = np.vstack(implied_timescales)

for i in range(implied_timescales.shape[1]):
    plt.plot(lag_times, implied_timescales[:, i],
             label='$\lambda_{%s}$' % (i+1))
    plt.show()	# Missing in the tutorial

#Step 5 Creating the MSM (remove extra ')' in the tutorial, check with Justin!).

m = msm.MSM(
    lag_time=10,
    method=msm.builders.transpose) # removed the extra ')'.
m.fit(assigs)

#Step 6 saving the MSM. STUCK!
with open('my-file-name.pkl', 'wb') as f:
	pickle.dump(m, f)

with open('my-file-name.pkl', 'rb') as f:
	m = pickle.load(f)

def save_msm_data(outputName, msm_object):
   pickleSave(outputName, msm_object)

   tcounts = msm_object.tcounts_
   tprobs = msm_object.tprobs_
   eqProbs = msm_object.eq_probs_

   np.savetxt("eq_probs.dat", eqProbs)





