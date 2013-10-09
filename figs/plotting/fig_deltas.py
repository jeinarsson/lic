import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as colormaps
import matplotlib.backends.backend_pdf
import time
import scipy.io
from os.path import exists
from os import mkdir
import h5py
from analyze_helpers import *

def SSS_theory_deltas():
	return [1,0,0,-0.5,0,-0.5,0.375,0,0,0,0,0.375,0,0,0,-0.5,0,1,-0.75,0,-0.5,0, \
   0,0.375,0,-0.75,0,0,0.5625,0,0,0,0,0,0,0,0,-0.75,0,0.375,0,0,0,0,0,0,1, \
   0,-0.5,0.375,0,-0.5,0,0.375,0,1]


def main():

	matplotlib.rc('text', usetex=True)
	matplotlib.rc('ps', usedistiller='xpdf')


	symmetric_indices = [0, 1, 2, 4, 5, 8]
	asymmetric_indices = [1, 2, 5]

	outputfile_name = 'turbulence/Acorr3_onetime.hdf5'
	output_file = h5py.File(outputfile_name, 'r')


	kolmogorov_time = .0446
	dataset = output_file['/A_0_0_0']
	delta_t = dataset.attrs['delta_t']
	time_downsampling = dataset.attrs['time_downsampling']
	max_lag = dataset.attrs['max_lag']
	times = np.arange(-max_lag, max_lag+1) * delta_t * time_downsampling / kolmogorov_time
	num_groups = dataset.shape[1]	

	time_slice = slice(max_lag, max_lag+1)

	norm = np.mean(get_Acorr3_onetime(output_file, 0,0,0,time_slice))
	meanvalues = []
	deviations = []
	count = 0
	ticks = []
	ticklabels = []
	for i1 in xrange(6):
		for i2 in xrange(i1, 6):
			for i3 in xrange(i2, 6):
				ij = symmetric_indices[i1]
				kl = symmetric_indices[i2]
				mn = symmetric_indices[i3]

				ji = transpose_index(ij)
				lk = transpose_index(kl)
				nm = transpose_index(mn)

				data1 = get_Acorr3_onetime(output_file, ij, kl, mn, time_slice)
				data2 = get_Acorr3_onetime(output_file, ij, kl, nm, time_slice)
				data3 = get_Acorr3_onetime(output_file, ij, lk, mn, time_slice)
				data4 = get_Acorr3_onetime(output_file, ij, lk, nm, time_slice)
				data5 = get_Acorr3_onetime(output_file, ji, kl, mn, time_slice)
				data6 = get_Acorr3_onetime(output_file, ji, kl, nm, time_slice)
				data7 = get_Acorr3_onetime(output_file, ji, lk, mn, time_slice)
				data8 = get_Acorr3_onetime(output_file, ji, lk, nm, time_slice)

				# SSS
				data = data1 + data2 + data3 + data4 + data5 + data6 + data7 + data8
				data /= 8.*norm

				meanvalues.append(np.mean(data))
				deviations.append(np.std(data))


	
				if count % 5 == 0:
					ticks.append(count + 0.5)
					ticklabels.append('$S_{{{}}}S_{{{}}}S_{{{}}}$'.format(index_name(ij), index_name(kl), index_name(mn)))

				print '$S_{{{}}}S_{{{}}}S_{{{}}}$'.format(index_name(ij), index_name(kl), index_name(mn))
				count += 1


	theory_deltas = SSS_theory_deltas()
	diffs = [a-b for (a,b) in zip(meanvalues, theory_deltas)]

	fig1, axes = plt.subplots(2,1, sharex=True,sharey=True,figsize=(11/2.5, 8/2.5))
	axes[0].bar(np.arange(len(meanvalues)), meanvalues, 1, color='b', yerr=deviations)
	axes[1].bar(np.arange(len(meanvalues)), theory_deltas, 1, color='b')
	#axes[2].bar(np.arange(len(meanvalues)), diffs, 1, color='b')
	#axes[2].set_xlabel('Component')
	axes[0].set_xlim(0,len(meanvalues))
	axes[0].set_title(r'\normalsize Turbulence data')
	axes[1].set_title(r'\normalsize Theory')
	axes[0].set_ylabel(r'$C^{SSS}_{ijklmn}/C^{SSS}_{111111}$')
	axes[1].set_ylabel(r'$C^{SSS}_{ijklmn}/C^{SSS}_{111111}$')
	#axes[2].set_ylabel('Difference')
	
	axes[1].set_xticks(ticks)
	axes[1].set_xticklabels(ticklabels,rotation=70)

	fig1.tight_layout()


	print "Starting export..."
	if not exists("output/"):
		mkdir("output")
	pp = matplotlib.backends.backend_pdf.PdfPages('output/deltas.pdf')
	pp.savefig(fig1)
	pp.close()


	plt.show()
	print "Done!"



if __name__ == '__main__':
	main()