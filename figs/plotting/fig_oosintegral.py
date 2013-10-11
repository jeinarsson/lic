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
from scipy.integrate import cumtrapz
from analyze_helpers import *

def OOS_theory_deltas():
	return [2,0,0,2,0,-4,0,0,0,0,3,0,0,0,-3,0,0,0,2,0,0,-4,0,2,0,3,0,0,0,0,-4,0,0,2,0,2]


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

	theory_deltas = OOS_theory_deltas()
	fig1, axes = plt.subplots(2,1, sharex=True,figsize=(11/2.5, 10/2.5))

	
	corr_function = np.zeros( (2*max_lag+1,))

	time_slice = slice(max_lag, max_lag+1)

	# Kristians norm: A11A11A11 = 16 at t=0
	#norm = np.mean(get_Acorr3_onetime(output_file, 0,0,0,time_slice))/16./16.
	norm = 1.

	num_functions = 0
	counter = -1

	for i1 in xrange(3):
		for i2 in xrange(3):
			for i3 in xrange(3):

				i1i2 = (i1+1, i2+1)
				i2i3 = (i2+1, i3+1)
				i3i1 = (i3+1, i1+1)

				# print i1i2
				# print i2i3
				# print i3i1
				# print ""

				ij = square_to_linear( (i1+1, i2+1) )
				jk = square_to_linear( (i2+1, i3+1) )
				ki = square_to_linear( (i3+1, i1+1) )

				print ij
				print jk
				print ki
				print ""

				kl = jk
				mn = ki

				ji = transpose_index(ij)
				lk = transpose_index(kl)
				nm = transpose_index(mn)

				data1 = get_Acorr3_onetime(output_file, ij, kl, mn)
				data2 = get_Acorr3_onetime(output_file, ij, kl, nm)
				data3 = get_Acorr3_onetime(output_file, ij, lk, mn)
				data4 = get_Acorr3_onetime(output_file, ij, lk, nm)
				data5 = get_Acorr3_onetime(output_file, ji, kl, mn)
				data6 = get_Acorr3_onetime(output_file, ji, kl, nm)
				data7 = get_Acorr3_onetime(output_file, ji, lk, mn)
				data8 = get_Acorr3_onetime(output_file, ji, lk, nm)

				# OOS
				data = data1 + data2 - data3 - data4 - data5 - data6 + data7 + data8

				# SSS
				#data = data1 + data2 + data3 + data4 + data5 + data6 + data7 + data8

				data /= 8.*norm

				meandata = np.mean(data, axis=1)
				corr_function += meandata
				
				num_functions += 1

				#plot_mean_and_deviation(axes[0], times, -data)
#				plot_mean_and_deviation_nonlog(axes[0], times, data)
				#axes[0].plot(times, meandata)




	
	print num_functions
	print norm

	corr_function *= kolmogorov_time**3
	#corr_function /= 1743.

	axes[0].plot(times, corr_function)

	integration_maxtime = 27 * kolmogorov_time
	integration_maxindex = int(np.ceil(integration_maxtime / (delta_t * time_downsampling)))

	#integration_slice = slice(max_lag,max_lag+integration_maxindex+1)
	#integration_slice = slice(max_lag, max_lag-integration_maxindex-1, -1)
	integration_slice = slice(max_lag-integration_maxindex,max_lag+integration_maxindex+1, 1)

	integration_times = times[integration_slice]
	integration_f = corr_function[integration_slice]
	
	#axes[0].plot(integration_times, integration_f, ls='dashed')


	integration = cumtrapz(integration_f, integration_times, initial=0.0)
	
	#value = np.mean(integration[-20:-1])
	value = integration[integration_maxindex]
	plt.axhline(axes=axes[1], y=value, ls='dotted', color='#c0c0c0')
	plt.axvline(axes=axes[1], x=0, ls='dotted', color='#c0c0c0')
	axes[1].plot(integration_times, integration)
	axes[1].text(5, value, '{:.3f}'.format(value))
	#axes[1].set_yscale('log')

	axes[0].set_ylabel(r'$\langle \mathrm{Tr} \mathbb{O}_0\mathbb{O}_0\mathbb{S}_t\rangle\tau^3_{K}$')
	axes[1].set_ylabel(r"$\int_{-\infty}^t\!\!\!\mathrm{d} t' \langle \mathrm{Tr} \mathbb{O}_0\mathbb{O}_0\mathbb{S}_{t'}\rangle\tau^3_{K} $")
	axes[0].set_xlim(-30,30)
	#axes[0].set_ylim(-30,12)
	axes[0].set_title(r'\normalsize Correlation function')
	axes[1].set_title(r'\normalsize Cumulative integral')
	axes[1].set_xlabel('$t/\\tau_{K}$')

	axes[0].set_yticks(np.arange(-.02, .07,.02))
	axes[1].set_yticks(np.arange(0., 0.4,.1))

	fig1.tight_layout()


	print "Starting export..."
	if not exists("output/"):
		mkdir("output")
	pp = matplotlib.backends.backend_pdf.PdfPages('output/oosintegral.pdf')
	pp.savefig(fig1)
	pp.close()


	plt.show()
	print "Done!"



if __name__ == '__main__':
	main()