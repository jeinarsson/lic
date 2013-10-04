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
import matplotlib.gridspec as gridspec
from itertools import cycle



def main():

	matplotlib.rc('text', usetex=True)
	matplotlib.rc('ps', usedistiller='xpdf')

	fig = plt.figure(figsize=(10/2.5, 7/2.5))
	
	#gs = gridspec.GridSpec(1, 1, width_ratios=[1,1], wspace=0.02,hspace=0.5,bottom=0.2,left=0.08,right=0.99)
	
	ax1 = fig.add_subplot(1,1,1)

	datafile_name = "data/lic_fig1.h5"
	datafile = h5py.File(datafile_name, 'r')
	print datafile.keys()
	
	data = datafile['/sinsqrtheta'][:,:]
	pes = datafile['/pes'][:]
	ls = datafile['/ls'][:]
	#dt = float(datafile.attrs['dt'])

 	colors=['#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#FF7F00', '#FFFF33']
	ax = ax1

	for k in xrange(5):
		b = (ls[k]**2-1.)/(ls[k]**2+1.)
		ax.semilogx(pes,data[k,:], color=colors[k],label='$\\lambda={}$'.format(ls[k]))
		range1 = slice(-25,-1,None)
		range2 = slice(0,25,None)

		if ls[k]>=25.:
			ax.semilogx(pes[range1],pes[range1]*0.0 + 1. - 1.792/ls[k],ls='dotted', color=colors[k])
	
		ax.semilogx(pes[range2],2./3 + 1./630 * (b*pes[range2])**2,ls='dotted',color=colors[k])
	
	ax.semilogx(pes,1. - 0.974*pes**(-1./3),ls='dotted')
	


	#ax.legend(loc='lower left')
#	ax.set_xticks(xticks)
#	ax.set_xticklabels(['${}$'.format((xtick-xticks[0])*1000.0) for xtick in xticks])
#	ax.set_xlim(xmin, xmax)
	ax.set_ylim(1./3,1.02)
	ax.set_ylabel("$\\langle \\sin^2\\theta\\rangle$")
	ax.set_xlabel("$\mathrm{Pe}$")
#	ax.set_yticks([-1,0, 1])
#	ax.set_yticklabels(["$-1$", "$0$", "$1$"])







	plt.tight_layout()

	print "Starting export..."
	if not exists("output/"):
		mkdir("output")
	pp = matplotlib.backends.backend_pdf.PdfPages('output/sintheta1.pdf')
	pp.savefig(fig, dpi=300)
	pp.close()


	plt.show()
	print "Done!"



if __name__ == '__main__':
	main()