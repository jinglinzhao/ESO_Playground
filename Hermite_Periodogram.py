import matplotlib.pyplot as plt
import numpy as np
from astropy.stats import LombScargle
import glob


from PyAstronomy.pyTiming import pyPeriod


STAR 		= 'Gl628'
# STAR 		= 'HD103720'
# STAR 		= 'Gl358'
# STAR 		= 'Gl479'
# STAR 		= 'Gl581'
# STAR 		= 'Gl674'
# STAR 		= 'Gl176'
# STAR 		= 'Gl388'
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
Pend = (max(MJD) - min(MJD)) / 2
# Pend = 20

# Compute the GLS periodogram with default options.
# Choose Zechmeister-Kuerster normalization explicitly
clp = pyPeriod.Gls((MJD, RV_HARPS, RV_noise), norm = "ZK", Pbeg=Pbeg, Pend=Pend)
# Print helpful information to screen
clp.info()

# Define FAP levels of 10%, 5%, and 1%
fapLevels = np.array([0.5, 0.01, 0.001, 0.0001])
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
fontsize = 16

for i in range(ORDER):

	print('\n', 'Order', str(i))

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
	plt.ylim(0, 1.1*max(clp.power))
	ax.set_xscale('log')
	if STAR == 'Gl628':
		ax.axvline(x=4.9, color='k', linestyle = '--', linewidth=1.0, alpha=0.8)
		ax.axvline(x=17.9, color='k', linestyle = '--', linewidth=1.0, alpha=0.8)
		ax.axvline(x=217.2, color='k', linestyle = '--', linewidth=1.0, alpha=0.8)
		ax.axvline(x=94, color='k', linestyle = '--', linewidth=1.0, alpha=0.8)
		ax.text(4.9, -0.02, '$P_b$', horizontalalignment='center', fontsize=fontsize)
		ax.text(17.9, -0.02, '$P_c$', horizontalalignment='center', fontsize=fontsize)
		ax.text(217.2, -0.02, '$P_d$', horizontalalignment='center', fontsize=fontsize)
		ax.text(94.2, 0.25, '$P^{*}$', fontsize=fontsize)
		plt.ylim(0, 0.27)
		ax.fill_between(1/clp.freq, 1, 0, where=(1/clp.freq>84.2)&(1/clp.freq<106.3), facecolor='red', alpha=0.2)
	elif STAR == 'HD103720':
		ax.axvline(x=4.5557, color='k', linestyle = '--', linewidth=1.0, alpha=0.8)
		ax.axvline(x=17, color='k', linestyle = '--', linewidth=1.0, alpha=0.8)
		ax.axvline(x=68, color='k', linestyle = '--', linewidth=1.0, alpha=0.8)
		ax.text(4.5557, max(clp.power), '$P_b$', fontsize=fontsize)
		ax.text(17, max(clp.power), '$P^{*}$', fontsize=fontsize)
		ax.text(68, max(clp.power), '$4P^{*}$', fontsize=fontsize)
	elif STAR == 'Gl674':
		ax.axvline(x=4.694, color='k', linestyle = '--', linewidth=1.0, alpha=0.8)
		ax.axvline(x=32.9, color='k', linestyle = '--', linewidth=1.0, alpha=0.8)
		ax.text(4.694, max(clp.power), '$P_b$', fontsize=fontsize)
		ax.text(32.9, max(clp.power), '$P^{*}$', fontsize=fontsize)	
	elif STAR == 'Gl479':
		ax.axvline(x=11.3, color='k', linestyle = '--', linewidth=1.0, alpha=0.8)
		ax.axvline(x=23.1, color='k', linestyle = '--', linewidth=1.0, alpha=0.8)
		ax.axvline(x=24.8, color='k', linestyle = '--', linewidth=1.0, alpha=0.8)
		ax.text(11.3, max(clp.power), '$P_b$', fontsize=fontsize)
		ax.text(23.1, max(clp.power), '$P_c$', fontsize=fontsize, horizontalalignment='right')
		ax.text(24.8, max(clp.power), '$P^{*}$', fontsize=fontsize)		
	elif STAR == 'Gl358':
		ax.axvline(x=13.16, color='k', linestyle = '--', linewidth=1.0, alpha=0.8)
		ax.axvline(x=26.27, color='k', linestyle = '--', linewidth=1.0, alpha=0.8)
		ax.axvline(x=26.8, color='k', linestyle = '--', linewidth=1.0, alpha=0.8)
		ax.text(13.16, max(clp.power), '$P_b$', fontsize=fontsize)
		ax.text(26.27, max(clp.power), '$P_c$', fontsize=fontsize, horizontalalignment='right')
		ax.text(26.8, max(clp.power), '$P^{*}$', fontsize=fontsize)
		# ax.fill_between(1/clp.freq, 1, 0, where=(1/clp.freq>26)&(1/clp.freq<27.6), facecolor='red', alpha=0.2)
	plt.xlim(Pbeg, Pend)

	for j in range(len(fapLevels)):
		if fapLevels[j] == 0.5:
			plt.plot([min(1/clp.freq), max(1/clp.freq)], [plevels_c[j]]*2, '--', label="FAP = %d%%" % (fapLevels[j]*100))
		else:
			plt.plot([min(1/clp.freq), max(1/clp.freq)], [plevels_c[j]]*2, '--', label="FAP = %4.2f%%" % (fapLevels[j]*100))

		np.savetxt('../' + STAR + '/power' + str(i) + '.out', clp_c.power)

	plt.legend()
	plt.savefig('../' + STAR + '/' + STAR + '_Order' + str(i) +'.png')
	plt.close()
	# plt.show()

np.savetxt('../' + STAR + '/frequency.out', clp_c.freq)
np.savetxt('../' + STAR + '/FAP_50.out', plevels_c)


###################
# Window Function #
###################
if STAR == 'Gl628':
	star_name = 'Wolf 1061'
elif STAR[0:2] == 'Gl':
	star_name = 'GJ ' + STAR[2:]
else:
	star_name = STAR[0:2] + ' ' + STAR[2:]
print('\n', 'Window Function for ', star_name)
clp = pyPeriod.Gls((MJD, np.ones(len(MJD)), np.ones(len(MJD))), norm = "ZK", Pbeg=Pbeg, Pend=Pend)
clp.info()
plevels = clp.powerLevel(fapLevels)
plt.rcParams.update({'font.size': 14})
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
ax = plt.subplot(111)
plt.plot(1/clp.freq, clp.power, 'k-', linewidth=1, alpha=1)
plt.title('Window Function for ' + star_name)
plt.xlabel('Period [days]')
plt.ylabel("Power")
plt.xlim(Pbeg, Pend)
plt.ylim(0, 1.1*max(clp.power))
ax.set_xscale('log')
plt.savefig('../' + STAR + '/' + STAR + 'window.png')
plt.close()
# plt.show()