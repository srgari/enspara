python /home/banff/enspara/enspara/apps/cluster.py --trajectories trajectory-*.xtc --topology fs-peptide.pdb --algorithm khybrid --cluster-number 20 --subsample 10 --atoms '(name N or name C or name CA)' --distances /home/banff/Desktop/MSM/fs_peptide/fs-khybrid-clusters0020-distances.h5 --center-features /home/banff/Desktop/MSM/fs_peptide/fs-khybrid-clusters0020-centers.pickle --assignments /home/banff/Desktop/MSM/fs_peptide/fs-khybrid-clusters0020-assignments.h5
############################################################################################################
#otherparameters I used
#--lag-times 5:100:2 

#All the parameters in Enspara documentation.
#
#--n-eigenvalues integer
# This is the number of eigenvalues that will be computed for each lag time.
# The default is five.
#
#--lag-times min:max:step
# This is the list of lag times (in frames).
# The default is 5:100:2. help="List of lagtimes (in frames) to compute eigenspectra for." 
# "Format is min:max:step.")
#
#--symmetrization method name
# This is the method to use to enforce detailed balance in the counts matrix.
# The default is transpose.
#
#--trj-ids trajs
# This will only use given trajectories to compute the implied timescales.
# This is useful for handling assignments for shared state space clusterings.
# The deafult is none.
#
#--processes integer
# This will set the number of cores to use.
# Because eigenvector decompositions are thread-parallelized, this should
# usually be several times smaller than the number of cores availiable on
# your machine.
# The deafult is max(1, cpu_count()/4).
#
#--trim truth statement
# This will turn on ergodic trimming
# The default is False.
#
#--plot path/to/directory/file_name.png
# This is how the plot will save.
#
#--logscale
# This will put the y-axis of the plot on a log scale.
#
#--timestep", default=None, type=float,
#Value to convert number of frames to ns (time). Ex for my case: 50 ps =1 frame, thus 20 * 50 ps = 1 ns

#################################################################################################################
#Lagtime plot
python /home/banff/enspara/enspara/apps/implied_timescales.py --assignments fs-khybrid-clusters0020-assignments.h5 --n-eigenvalues 5  --processes 2 --plot implied_timescales.png --logscale
# check implied_timescale.py for more option!
