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
	


	outputfile_name = 'turbulence/alignment_pdfs.hdf5'
	output_file = h5py.File(outputfile_name, 'r')


	
	lambdas = output_file['/lambdas'][:]
	

	xvalues = np.linspace(0,3,50)

	for dataset in [output_file['ndotsqr']]:
		pdf = dataset[0,:]
		pdf = pdf/np.trapz(pdf,xvalues)
		ax.semilogy(xvalues, pdf, label=r'$\lambda = 1/100$')
		pdf = dataset[10,:]
		pdf = pdf/np.trapz(pdf,xvalues)
		ax.semilogy(xvalues, pdf, label=r'$\lambda = 100$')
		pdf = dataset[5,:]
		pdf = pdf/np.trapz(pdf,xvalues)
		ax.semilogy(xvalues, pdf, label=r'$\lambda = 1$')

	print lambdas[0]
	print lambdas[5]
	print lambdas[10]
	ax.legend()


	ax.set_xlabel(r'$\dot n^2 \tau^2_{K}$')
	ax.set_ylabel(r'$P(\dot n^2 \tau^2_{K})$',rotation=90)
	#ax.set_ylim(0,1)
	ax.set_xlim(0,2.9)

	fig1.tight_layout()

	print "Starting export..."
	if not exists("output/"):
		mkdir("output")
	pp = matplotlib.backends.backend_pdf.PdfPages('output/ndotsqr_pdf.pdf')
	pp.savefig(fig1)
	pp.close()


	plt.show()
	print "Done!"



if __name__ == '__main__':
	main()