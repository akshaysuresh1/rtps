# rtps
Plot a phase space diagram of radio transients using Python. For a C implementation of the code, see [Evan Keane's repository](https://github.com/FRBs/Transient_Phase_Space).<br>

Plot design is motivated by Figure 1 of [Cordes et al. (2004)](https://ui.adsabs.harvard.edu/abs/2004NewAR..48.1459C/abstract) and Figure 5 of [Pietka et al. (2015)](https://ui.adsabs.harvard.edu/abs/2015MNRAS.446.3687P/abstract).

![Phase space of radio transients](https://github.com/akshaysuresh1/rtps/blob/main/Plots/rtps.png?raw=True)

In the above figure, colored markers indicate various source classes. The uncertainty principle (gray shaded box) restricts radio transient discovery below $\nu_{\rm GHz}W_{\rm s} \leq 5 \times 10^{-10}$. Diagonal dotted lines label contours of constant brightness temperature, $T_B$. The line $T_B  = 10^{12}\ {\rm K}$ marks the theoretical maximum brightness temperature of an incoherently emitting synchrotron radio source. Above this value, inverse Compton losses significantly limit the emission flux density. Note that while emissions with $T_B \gg 10^{12}\ {\rm K}$  are strictly coherent, $T_B \lesssim 10^{12}\ {\rm K}$  does not imply incoherent emission. <br>  

---

## Data
Plotted data for different source classes are available as `.csv` files under the `Data` directory of this repository. <br>

Presently, we take data points from the following papers and the references therein.
1. [Pietka et al. (2015)](https://ui.adsabs.harvard.edu/abs/2015MNRAS.446.3687P/abstract) (see Figure 5)
2. [Nimmo et al. (2022)](https://ui.adsabs.harvard.edu/abs/2022NatAs...6..393N/abstract) (see Figure 3)
3. [Hurley-Walker et al. (2022)](https://ui.adsabs.harvard.edu/abs/2022Natur.601..526H/abstract) (see Figure 4)
4. [Caleb et al. (2022)](https://ui.adsabs.harvard.edu/abs/2022NatAs.tmp..123C/abstract) (see Table 1)

## Installation and usage  
Execute the following steps in sequence, provided you have a working Python3 installation including the `numpy`, `matplotlib` and `astropy` packages.
```
git clone git@github.com:akshaysuresh1/rtps.git
cd rtps
python rtps.py
```

## Troubleshooting <a name="troubleshooting"></a>
Please submit an issue to voice any problems or requests. Assistance with keeping data up-to-date is greatly appreciated.
