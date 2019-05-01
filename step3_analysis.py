#Step 7 Analysis (change the input file name)
import pickle
import mdtraj as md

with open('fs-khybrid-clusters0020-centers.pickle', 'rb') as f:
    ctr_structs = md.join(pickle.load(f))

#Step 8
hbonds = md.kabsch_sander(ctr_structs)

#Step 9
import matplotlib.pylab as plt
weighted_hbond_mtx = sum(p*h for p, h in zip(m.eq_probs_, hbonds)).todense()
plt.imshow(weighted_hbond_mtx, cmap='viridis_r')
plt.colorbar()



all_hbonds = set()

#Step11
# accumulate all the possible pairs of residues involved in hbonds
for i in range(len(ctr_structs)):
    donors, acceptors = np.where(hbonds[i].todense() != 0)
    all_hbonds.update([(d, a) for d, a in zip(donors, acceptors)])

# make a list so that it's ordered
all_hbonds = list(all_hbonds)

# this matrix of length n_states will have each binary feature vector
hbond_presence = np.zeros((m.n_states_, len(all_hbonds)),
                          dtype='uint8')

# set each value i, j to one if state i has hbond j.
for i in range(len(ctr_structs)):
    donors, acceptors = np.where(hbonds[i].todense() != 0)

    for a, d in zip(donors, acceptors):
        hbond_id = all_hbonds.index((a, d))
        hbond_presence[i, hbond_id] = 1

#Step 12
p_hbond = np.dot(m.eq_probs_, hbond_presence)

plt.bar(np.arange(len(all_hbonds)), height=np.log(p_hbond))
plt.ylabel("Free Energy")
plt.xlabel("HBond ID")
plt.savefig('./hbond-free-energy.svg')

#Step 12
from enspara.info_theory import weighted_mi

hbond_mi = weighted_mi(features=hbond_presence, weights=m.eq_probs_)
hbond_mi = hbond_mi - np.diag(np.diag(hbond_mi))

plt.imshow(hbond_mi - np.diag(np.diag(hbond_mi)))
plt.colorbar()

#Step 13

hbond1, hbond2 = np.unravel_index(hbond_mi.argmax(), hbond_mi.shape)

def hbond2str(pair, top):
   return '⟶'.join([str(top.residue(i)) for i in pair])

hbond2str(all_hbonds[hbond1], ctr_structs.top), hbond2str(all_hbonds[hbond2], ctr_structs.top)

