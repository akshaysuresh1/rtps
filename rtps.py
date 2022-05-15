# Plot the radio transient phase space diagram.

# Load relevant packages.
import os
import astropy.units as u
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Enable use of LaTeX labels in plots.
plt.rcParams.update({"text.usetex": True, "text.latex.preamble": r"\usepackage{amsmath}"})
#####################################################################
# INPUTS

# Name (including path) of .csv files for different source classes
# Pulsars
psr_file = 'Data/pulsars.csv'
# Rotating radio transients (RRATs)
rrat_file = 'Data/rrats.csv'
# Crab giant pulses
crab_gp_file = 'Data/crab_gp.csv'
# Crab nanoshots
crab_nano_file = 'Data/crab_nano.csv'
# SGR 1935+2154
sgr1935_file = 'Data/sgr1935.csv'
# Fast radio bursts
frb_file = 'Data/frb.csv'
# GLEAM-X J162759.5-523504.3
gleamx_file = 'Data/gleamx.csv'
# Solar radio bursts
solar_file = 'Data/solar_bursts.csv'
# Supernovae
sn_file = 'Data/sn.csv'
# AGN/Blazar/QSO
agn_file = 'Data/agn.csv'
# Gamma-ray burst afterglows
grb_file = 'Data/grb.csv'
# X-ray binaries
xrb_file = 'Data/xrb.csv'
# Novae
novae_file = 'Data/novae.csv'
# Flare stars/brown dwarfs
flarestars_file = 'Data/flarestars.csv'
# RS Canum Venaticorum variables
rscvn_file = 'Data/rscvn.csv'
# Magnetic Cataclysmic Variables
magcv_file = 'Data/magcv.csv'

# Data for miscellaneous sources (e.g., AR Sco, GCRT 1745, MKT J1704, Jupiter DAM, GW 170817)
# are typed directly into their respective plotting commands.

# Output parameters
show_plot = False # Do you want to view the output plot live? (True/False)
save_plot = True # Do you want to save the output plot to disk? (True/False)
basename = 'rtps' # Output basename
plot_formats = ['.png', '.pdf'] # List of plot formats
OUTPUT_DIR = 'Plots' # Output path (will be created if non-existent)
#####################################################################
# FUNCTION DEFINITIONS

# Compute pseudoluminosity (Jy kpc^2) power law curves for brightness temperature TB.
def calc_Lnu_B(TB, vW):
    '''
    Inputs:
    TB: Brightness temperature (K)
    vW: Product of radio frequency (GHz) and transient duration (s)

    Output:
    L: Spectral pseudo-luminosity (Jy kpc^2)

    Conversions:
    1 W/Hz = (1.05e-13) Jy kpc^2

    Formula (Equation 2 of Pietka et al., 2015 with the 4\pi solid angle factor removed):
    For radio frequency v in GHz and pulse duration W in seconds,
    L = TB * (vW)**2 * (2.761e-5) Watts/Hz
      = TB * (vW)**2 * (2.761 * 1.05026e-18) Jy kpc^2
    '''
    L = TB*(vW**2)*2.761*1.05026e-18
    return L

# Convert erg s^{-1} Hz^{-1} to Jy kpc^2.
def cgs_to_astro(L_cgs):
    L_astro = L_cgs * 1e19/(u.kpc.to(u.m))**2
    return L_astro

# Convert Jy kpc^2 to erg s^{-1} Hz^{-1}.
def astro_to_cgs(L_astro):
    L_cgs = L_astro * 1e-19 * (u.kpc.to(u.m))**2
    return L_cgs
#####################################################################
# PLOTTING
if show_plot or save_plot:
    fig = plt.figure(figsize=(10,8), tight_layout=True)

    # Lines of constant brightness temperature
    TB = np.array([1e4,1e8,1e12,1e16,1e20,1e24,1e28,1e32,1e36,1e40])
    vW=np.geomspace(1e-10,1e10,1000)
    for tb in TB:
        plt.plot(vW,calc_Lnu_B(tb,vW), color='k', lw=1, alpha=0.3, linestyle=':')

    # TB labels
    # Slope of TB curves in log-log space = 2
    rot_angle = 50.5 # Rotation angle (degree)
    plt.annotate(r'$10^{40}$ K', xy=(9e-10,1.5e5), xycoords='data', rotation=rot_angle, color='k', fontsize=13)
    plt.annotate(r'$10^{36}$ K', xy=(9e-10,15), xycoords='data', rotation=rot_angle, color='k', fontsize=13)
    plt.annotate(r'$10^{32}$ K', xy=(9e-10,1.5e-3), xycoords='data', rotation=rot_angle, color='k', fontsize=13)
    plt.annotate(r'$10^{28}$ K', xy=(9e-10,1.5e-7), xycoords='data', rotation=rot_angle, color='k', fontsize=13)
    plt.annotate(r'$10^{24}$ K', xy=(6e-9,7e-10), xycoords='data', rotation=rot_angle, color='k', fontsize=13)
    plt.annotate(r'$10^{20}$ K', xy=(6e-7,7e-10), xycoords='data', rotation=rot_angle, color='k', fontsize=13)
    plt.annotate(r'$10^{16}$ K', xy=(6e-5,7e-10), xycoords='data', rotation=rot_angle, color='k', fontsize=13)
    plt.annotate(r'$10^{12}$ K', xy=(1e-1,2e-7), xycoords='data', rotation=rot_angle, color='k', fontsize=13)
    plt.annotate(r'$10^{8}$ K', xy=(6e-1,7e-10), xycoords='data', rotation=rot_angle, color='k', fontsize=13)
    plt.annotate(r'$10^{4}$ K', xy=(6e1,7e-10), xycoords='data', rotation=rot_angle, color='k', fontsize=13)

    # Show light blue shading below incoherent emission TB upper bound of 10^{12} K.
    plt.fill_between(vW, np.ones(len(vW))*1.0e-10, calc_Lnu_B(1e12,vW), color='lightcyan')

    # Pulsars
    psrs = pd.read_csv(psr_file)
    plt.scatter(psrs['vW (GHz s)'], psrs['L (Jy kpc^2)'], c='blue', marker='o', alpha=0.8)
    plt.text(1.35e-4, 4.04e-6, 'Pulsars', color='blue', fontsize=14)

    # RRATs
    rrat = pd.read_csv(rrat_file)
    plt.scatter(rrat['vW (GHz s)'], rrat['L (Jy kpc^2)'], color='red', marker='o', alpha=0.8)
    plt.text(9e-3, 92, 'RRATs', color='red', fontsize=14)

    # Crab giant pulses
    crab_gp = pd.read_csv(crab_gp_file)
    plt.scatter(crab_gp['vW (GHz s)'], crab_gp['L (Jy kpc^2)'], color='slateblue', marker='o', alpha=0.8)
    plt.text(5e-8, 8.3, 'Crab GPs', color='slateblue', fontsize=14)

    # Crab nano-shots (Hankins+2003 and Jessner+2010)
    crab_nano = pd.read_csv(crab_nano_file)
    plt.scatter(crab_nano['vW (GHz s)'], crab_nano['L (Jy kpc^2)'], color='maroon', marker='o', alpha=0.8)
    plt.text(8e-10, 1305, 'Crab nanoshots', color='maroon', fontsize=14)

    # SGR 1935+2154
    sgr = pd.read_csv(sgr1935_file)
    plt.scatter(sgr['vW (GHz s)'], sgr['L (Jy kpc^2)'], color='purple', marker='o', alpha=0.8)
    plt.text(1.2e-3, 78805.0, r'SGR 1935$+$2154', color='purple', fontsize=14)

    # FRBs
    frb = pd.read_csv(frb_file)
    plt.scatter(frb['vW (GHz s)'], frb['L (Jy kpc^2)'], color='orangered', marker='o', alpha=0.9)
    plt.text(2.44e-4, 2.55e14, 'Localized fast radio bursts', color='orangered', fontsize=14)

    # GLEAM-X J162759.5-523504.3
    gleamx = pd.read_csv(gleamx_file)
    plt.scatter(gleamx['vW (GHz s)'], gleamx['L (Jy kpc^2)'], color='k', marker='o', alpha=0.9)
    plt.text(1.1, 1e-1, r'GLEAM-X'+'\n'+r'J162759.5', color='k', fontsize=14)

    # Solar radio bursts
    solar_bursts = pd.read_csv(solar_file)
    plt.scatter(solar_bursts['vW (GHz s)'], solar_bursts['L (Jy kpc^2)'], color='orange', marker='o', alpha=0.9)
    plt.text(8, 3e-6, r'Solar bursts', color='orange', fontsize=14)

    # Supernovae
    sn = pd.read_csv(sn_file)
    plt.scatter(sn['vW (GHz s)'], sn['L (Jy kpc^2)'], color='#808080', marker='o', alpha=0.9)
    plt.text(1e4, 1e8, r'Supernovae', color='#808080', fontsize=14)

    # AGN
    agn = pd.read_csv(agn_file)
    plt.scatter(agn['vW (GHz s)'], agn['L (Jy kpc^2)'], color='navy', marker='o', alpha=0.9)
    plt.text(3e3, 1e15, r'AGN/Blazar/QSO', color='navy', fontsize=14)

    # GRB afterglows
    grb = pd.read_csv(grb_file)
    plt.scatter(grb['vW (GHz s)'], grb['L (Jy kpc^2)'], color='rebeccapurple', marker='o', alpha=0.9)
    plt.text(5e2, 1e12, r'$\gamma$-ray burst afterglows', color='rebeccapurple', fontsize=14)

    # X-ray binaries
    xrb = pd.read_csv(xrb_file)
    plt.scatter(xrb['vW (GHz s)'], xrb['L (Jy kpc^2)'], color='darkgoldenrod', marker='o', alpha=0.9)
    plt.text(1e3, 2e4, r'X-ray binaries', color='darkgoldenrod', fontsize=14)

    # Novae
    novae = pd.read_csv(novae_file)
    plt.scatter(novae['vW (GHz s)'], novae['L (Jy kpc^2)'], color='palevioletred', marker='o', alpha=0.9)
    plt.text(1e8, 1, r'Novae', color='palevioletred', fontsize=14)

    # Flare stars/brown dwarfs
    fsbd = pd.read_csv(flarestars_file)
    plt.scatter(fsbd['vW (GHz s)'], fsbd['L (Jy kpc^2)'], color='brown', marker='o', alpha=0.9)
    plt.text(1e3, 6e-9, r'Flare stars/brown dwarfs', color='brown', fontsize=14)

    # RSCVn
    rscvn = pd.read_csv(rscvn_file)
    plt.scatter(rscvn['vW (GHz s)'], rscvn['L (Jy kpc^2)'], color='#8B0000', marker='o', alpha=0.9)
    plt.text(5e6, 6e-3, r'RSCVn', color='#8B0000', fontsize=14)

    # Mag CV
    magcv = pd.read_csv(magcv_file)
    plt.scatter(magcv['vW (GHz s)'], magcv['L (Jy kpc^2)'], color='#228B22', marker='o', alpha=0.9)
    plt.text(6e3, 2e-5, r'Magnetic CV', color='#228B22', fontsize=14)

    # Miscallenous sources
    # Ar Sco
    plt.scatter(0.016, 1053.0, color='orchid', marker='o', alpha=0.9)
    plt.text(3e-2, 2e3, r'AR Sco', color='orchid', fontsize=14)
    # GCRT J1745
    plt.scatter(200, 200, color='teal', marker='o', alpha=0.9)
    plt.text(10, 8e2, r'GCRT J1745$-$3009', color='teal', fontsize=14)
    # MKT J1704
    plt.scatter(2570400, 0.214775, color='darkslategray', marker='o', alpha=0.9)
    plt.text(6e6, 9e-2, r'MKT J170456.2', color='darkslategray', fontsize=14)
    # Jupiter DAM
    plt.scatter(4.0e-3, 5.0e-9, color='saddlebrown', marker='o', alpha=0.9)
    plt.text(5e-4, 1e-9, r'Jupiter DAM', color='saddlebrown', fontsize=14)
    # GW 170817 afterglow
    plt.scatter(7.1e+7, 1.4e+5, color='indigo', marker='o', alpha=0.9)
    plt.text(2e7, 2e3, r'GW 170817'+'\n'+r'afterglow', color='indigo', fontsize=14)


    # Uncertainty principle (vW <= 5 \times 10^{-10} for v in GHz and W in s)
    min_allowed_vW = 5.0e-10
    plt.axvline(min_allowed_vW, linestyle='-', color='k')
    plt.axvspan(np.min(vW), min_allowed_vW, ymin=0.0, ymax=1.0, facecolor='darkgray')
    plt.text(1.5e-10, 10, 'Uncertainty principle', color='k', rotation=90, fontsize=14)

    # Axes labels
    plt.xlabel(r'Radio frequency (GHz) $\times$ Transient duration (s)', fontsize=16)
    plt.ylabel(r'Spectral pseudo-luminosity, $L_{\nu}$ (Jy kpc$^{-2}$)', fontsize=16)

    # Axes tick properites
    plt.gca().set_xscale('log')
    plt.gca().set_xlim((1e-10,1e10))
    plt.xticks(np.geomspace(1e-10,1e10,21),
               labels=[r'$10^{-10}$','',r'$10^{-8}$','',r'$10^{-6}$','',r'$10^{-4}$','',r'$10^{-2}$','',r'$1$',
                       '',r'$10^2$','',r'$10^4$','',r'$10^6$','',r'$10^8$','',r'$10^{10}$'],
               fontsize=14)
    plt.gca().set_yscale('log')
    plt.gca().set_ylim((1e-10,1e16))
    plt.yticks(np.geomspace(1e-10,1e16,27), fontsize=14,
               labels=[r'$10^{-10}$','',r'$10^{-8}$','',r'$10^{-6}$','',r'$10^{-4}$','',r'$10^{-2}$','',r'$1$',
                       '',r'$10^2$','',r'$10^4$','',r'$10^6$','',r'$10^8$','',r'$10^{10}$','',r'$10^{12}$','',
                       r'$10^{14}$','',r'$10^{16}$'])

    # Set up secondary y-axis for luminosity units of erg s^{-1} Hz^{-1}.
    sec_ax = plt.gca().secondary_yaxis('right', functions=(astro_to_cgs, cgs_to_astro))
    sec_ax.set_ylabel(r'Spectral pseudo-luminosity, $L_{\nu}$ (erg s$^{-1}$ Hz$^{-1}$)', fontsize=16)
    sec_ax.set_yticks(np.geomspace(1e10,1e36,27))
    sec_ax.set_yticklabels([r'$10^{10}$','','','','',r'$10^{15}$','','','','',r'$10^{20}$','','','','',
                            r'$10^{25}$','','','','',r'$10^{30}$','','','','',r'$10^{35}$',''])
    sec_ax.tick_params(axis='y', labelsize=14)

    # Save plot to disk.
    if save_plot:
        if not os.path.isdir(OUTPUT_DIR):
            os.makedirs(OUTPUT_DIR)
        for format in plot_formats:
            plt.savefig(OUTPUT_DIR+'/'+basename+format)
    # Show plot live.
    if show_plot:
        plt.show()
    else:
        plt.close()

# END OF CODE
#####################################################################
