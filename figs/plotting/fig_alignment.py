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

	outputfile_name = 'turbulence/alignment_pdfs.hdf5'
	output_file = h5py.File(outputfile_name, 'r')


	lambda_idc = [0, 10]
	lambdas = output_file['/lambdas'][:]
	lambdas = lambdas[lambda_idc]
	


	fig1, axes = plt.subplots(2,2, sharex=True,figsize=(11.5/2.5, 9/2.5))
	xvalues = np.linspace(0,1,50)

	data = [
		('n_e1', r'| \mbox{\boldmath$n$}\cdot \mbox{\boldmath$e$}_1|', axes[0,0]),
		('n_e2', r'|\mbox{\boldmath$n$}\cdot \mbox{\boldmath$e$}_2|', axes[0,1]),
		('n_e3', r'|\mbox{\boldmath$n$}\cdot \mbox{\boldmath$e$}_3|', axes[1,0]),
		('n_w', r'|\mbox{\boldmath$n$}\cdot \mbox{\boldmath$\Omega$}|', axes[1,1]),
		#('w_e1', '|e_1\\cdot \\Omega|', axes[2,0]),
		#('w_e2', '|e_2\\cdot \\Omega|', axes[2,1]),
		#('w_e3', '|e_3\\cdot \\Omega|', axes[2,2]),
	]

	for (datasetname, label, ax) in data:
		dataset = output_file[datasetname]
		for k, l in enumerate(lambdas):
			pdf = dataset[lambda_idc[k],:]
			pdf = pdf/np.trapz(pdf,xvalues)
			
			if k==0:
				lambdaname = '1/100'
			elif k==1:
				lambdaname = '100'
			else:
				lambdaname = '1'

			ax.plot(xvalues,pdf,label='$\\lambda = {}$'.format(lambdaname))
			

		ax.set_ylabel('$P({})$'.format(label))
		ax.set_xlabel('${}$'.format(label))
		ax.set_ylim([0,3.5])
		ax.set_yticks([0,1,2,3])

	axes[0,0].legend(loc=9)

	fig1.tight_layout()

	print "Starting export..."
	if not exists("output/"):
		mkdir("output")
	pp = matplotlib.backends.backend_pdf.PdfPages('output/alignment.pdf')
	pp.savefig(fig1)
	pp.close()


	plt.show()
	print "Done!"



if __name__ == '__main__':
	main()