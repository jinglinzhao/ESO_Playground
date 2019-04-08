import matplotlib.pyplot as plt
import numpy as np
from astropy.stats import LombScargle
import glob


from PyAstronomy.pyTiming import pyPeriod


# STAR 		= 'Gl628'
# STAR 		= 'HD103720'
STAR 		= 'Gl358'
MJD     	= np.loadtxt('../' + STAR + '/MJD.dat')
RV_HARPS 	= np.loadtxt('../' + STAR + '/RV_HARPS.dat')
RV_noise 	= np.loadtxt('../' + STAR + '/RV_noise.dat')

ORDER = 22

Hermite_coeff 	= np.zeros((ORDER, len(MJD)))
FILE        	= glob.glob('../' + STAR + '/Periodogram_*.txt')

for i in range(ORDER):
	Hermite_coeff[i, :]	= np.loadtxt(FILE[i])


#####################
# Periodogram of RV #
#####################
Pbeg = 2
Pend = max(MJD) - min(MJD)

# Compute the GLS periodogram with default options.
# Choose Zechmeister-Kuerster normalization explicitly
clp = pyPeriod.Gls((MJD, RV_HARPS, RV_noise), norm = "ZK", Pbeg=Pbeg, Pend=Pend)
# Print helpful information to screen
clp.info()

# Define FAP levels of 10%, 5%, and 1%
fapLevels = np.array([0.01, 0.001, 0.0001])
# Obtain the associated power thresholds
plevels = clp.powerLevel(fapLevels)


########################
# Periodogram of coeff #
########################
left  = 0.15  # the left side of the subplots of the figure
right = 0.95    # the right side of the subplots of the figure
bottom = 0.15   # the bottom of the subplots of the figure
top = 0.9      # the top of the subplots of the figure
wspace = 0.55   # the amount of width reserved for blank space between subplots
hspace = 0.2   # the amount of height reserved for white space between subplots
alpha 	= 0.5


for i in range(ORDER):

	plt.rcParams.update({'font.size': 14})
	plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
	clp_c = pyPeriod.Gls((MJD, Hermite_coeff[i, :], RV_noise/1000), norm = "ZK", Pbeg=Pbeg, Pend=Pend)
	clp_c.info()
	plevels_c = clp_c.powerLevel(fapLevels)

	ax = plt.subplot(111)
	plt.plot(1/clp.freq, clp.power, 'k-', linewidth=3, alpha=0.2)
	if i%2 == 0:
		plt.plot(1/clp_c.freq, clp_c.power, 'r-', linewidth=1, alpha=alpha)
	else:
		plt.plot(1/clp_c.freq, clp_c.power, 'b-', linewidth=1, alpha=alpha)

	plt.title(r'$h_{%d}$' %i)
	plt.xlabel('Period [days]')
	plt.ylabel("Power")
	ax.set_xscale('log')
	if STAR == 'Gl628':
		ax.axvline(x=4.9, color='k', linestyle = '--', linewidth=1.0, alpha=0.8)
		ax.axvline(x=17.9, color='k', linestyle = '--', linewidth=1.0, alpha=0.8)
		ax.axvline(x=217.2, color='k', linestyle = '--', linewidth=1.0, alpha=0.8)
		ax.axvline(x=94, color='k', linestyle = '--', linewidth=1.0, alpha=0.8)
		ax.text(4.9, -0.02, '$P_1$', horizontalalignment='center', fontsize=12)
		ax.text(17.9, -0.02, '$P_2$', horizontalalignment='center', fontsize=12)
		ax.text(217.2, -0.02, '$P_3$', horizontalalignment='center', fontsize=12)
		ax.text(94.2, 0.25, '$P^{*}$', fontsize=12)
		plt.ylim(0, 0.27)
		ax.fill_between(1/clp.freq, 1, 0, where=(1/clp.freq>84.2)&(1/clp.freq<106.3), facecolor='red', alpha=0.2)
	plt.xlim(Pbeg, Pend)

	for j in range(len(fapLevels)):
	    plt.plot([min(1/clp.freq), max(1/clp.freq)], [plevels[j]]*2, '--',
	    label="FAP = %4.2f%%" % (fapLevels[j]*100))

	    np.savetxt('../' + STAR + '/power' + str(i) + '.out', clp_c.power)

	plt.legend()
	plt.savefig('../' + STAR + '/' + STAR + '_Order' + str(i) +'.png')
	plt.close()
	# plt.show()

np.savetxt('../' + STAR + '/frequency.out', clp_c.freq)