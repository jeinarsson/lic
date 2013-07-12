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


from matplotlib.patches import Ellipse

def main():

	matplotlib.rc('text', usetex=True)
	matplotlib.rc('ps', usedistiller='xpdf')

	
	fig1 = plt.figure(1, figsize=(8/2.5, 7/2.5))
	ax = fig1.add_subplot(1,1,1)
	
	xsm = np.linspace(-15,0,100)
	xsp = np.linspace(0,15,100)
	ysm = 6.0**(1./3.) * (-xsm)**(2./3.)
	ysp = 6.0**(1./3.) * xsp**(2./3.)
	ax.plot(xsm,ysm, color='#E41A1C')
	ax.plot(xsp,ysp, color='#E41A1C')
	ax.plot([0,0],[0,-100], color='#E41A1C')
	#ax.axhline(y=0,ls='dotted',color='#777777')
	#ax.axvline(x=0,ls='dotted',color='#777777')
	ax.set_xlim(-15,15)
	ax.set_ylim(-15,15)
	ax.set_xlabel(r'$\mathrm{Tr} \mathbb{B}^3$')
	ax.set_ylabel(r'$\mathrm{Tr} \mathbb{B}^2$',rotation=0,labelpad=-5)
	ax.set_xticks([-10,0,10])
	ax.set_yticks([-10,0,10])
	
	ax.text(0,7.5,'Aligning',ha='center',va='center')
	ax.text(-7.5,-8,'Spiral out',ha='center',va='center')
	ax.text(7.5,-8,'Spiral in',ha='center',va='center')

	ax.annotate("$(\mathrm{Tr} \mathbb{B}^2)^3 = 6(\mathrm{Tr} \mathbb{B}^3)^2$",
                  xy=(7, 7), xycoords='data',
                  xytext=(5, 19), textcoords='data',
                  va="center", ha="center",
                  bbox=dict(boxstyle="round", fc="w",lw=.3),
                  arrowprops=dict(arrowstyle="-|>",
                                  connectionstyle="arc3,rad=0.2",
                                  relpos=(.5, 0.),
                                  fc="w"), 
                  )

	plt.subplots_adjust(top=0.8,left=0.2,right=0.8)

	print "Starting export..."
	if not exists("output/"):
		mkdir("output")
	pp = matplotlib.backends.backend_pdf.PdfPages('output/bmap.pdf')
	pp.savefig(fig1)
	pp.close()


	plt.show()
	print "Done!"



if __name__ == '__main__':
	main()