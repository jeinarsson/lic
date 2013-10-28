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

def main():

	matplotlib.rc('text', usetex=True)
	matplotlib.rc('ps', usedistiller='xpdf')

	fig1 = plt.figure(1, figsize=(10/2.5, 7/2.5))
	ax = fig1.add_subplot(1,1,1)
	


	outputfile_name = 'turbulence/particle_statistics_condAA_AAA_std1.hdf5'
	output_file = h5py.File(outputfile_name, 'r')


	kolmogorov_time = .0446
	dataset_ndotsqr = output_file['/ndotsqr']
	delta_t = output_file.attrs['delta_t']

	lambdas = output_file['/lambdas'][:]
	

	# white noise theory
	# outputfile_name = 'Acorr2.hdf5'
	# output_file = h5py.File(outputfile_name)
	# kolmogorov_time = .0446
	# dataset = output_file['/A_0_0']
	# delta_t = dataset.attrs['delta_t']
	# time_downsampling = dataset.attrs['time_downsampling']
	# max_lag = dataset.attrs['max_lag']
	# time_slice = slice(max_lag,max_lag+1)
	
	
	norm = 476.267330091 # trAtA()
	
	data = dataset_ndotsqr[:,:]/norm





	b = (lambdas**2-1.)/(lambdas**2+1.0)
	whitenoise = (5.0 + 3*b**2)/30.

	SOS_int = 0.02
	SSS_int = -.3
	OOS_int = +.36	

	thirdorder = whitenoise + 2.*b/5.*(3.*b**2/7.*SSS_int-OOS_int+2.*b*SOS_int)

	[ax._get_lines.color_cycle.next() for k in xrange(1)]

	ax.plot(lambdas, thirdorder,ls='solid', label=r'$\mathrm{Ku}^3$')
	ax.plot(lambdas, whitenoise,ls='solid', label=r'$\mathrm{Ku}^2$')
	

	meanvalues = np.mean(data, axis=1)
	deviations = np.std(data, axis=1)
	#ax.errorbar(lambdas, meanvalues, yerr=deviations, ls='none',marker='o',ms=4,markerfacecolor='#ffffff',color='#004B89' ,label='Turbulence data')
	ax.plot(lambdas, meanvalues, ls='none',marker='o',ms=5,markerfacecolor='#ffffff',color='#004B89' ,label='Data')
	ax.legend(loc='upper right')
	ax.set_ylim([0.05,.50])
	ax.set_xscale('log')

	ax.set_ylabel('$\\langle \\dot{n}^2 \\tau^2_{K} \\rangle$')
	ax.set_xlabel('$\\lambda$')
	fig1.tight_layout()

	print "Starting export..."
	if not exists("output/"):
		mkdir("output")
	pp = matplotlib.backends.backend_pdf.PdfPages('output/ndotsqr_mean.pdf')
	pp.savefig(fig1)
	pp.close()


	plt.show()
	print "Done!"



if __name__ == '__main__':
	main()