"""
Calculate Antoine Equation Coefficients for the vapor pressure of
methyl valerate.

This code prints the best-fit coefficients of the Antoine Equation utilizing
the data files available in this folder. The data have been collected from:

Juan Ortega, Fernando Espiau, Jose Tojo, Jose Canosa, and Ana Rodriguez.
"Isobaric Vapor-Liquid Equilibria and Excess Properties for the Binary
Systems of Methyl Esters + Heptane." Journal of Chemical & Engineering
Data 48, no. 5 (September 2003): 1183-1190. doi:10.1021/je030117d.

Sergey P. Verevkin, and Vladimir N. Emel'yanenko. "Transpiration Method:
Vapor Pressures and Enthalpies of Vaporization of Some Low-Boiling Esters."
Fluid Phase Equilibria 266, no. 1-2 (April 2008): 64-75.
doi:10.1016/j.fluid.2008.02.001.

Aad C.G. van Genderen, J.Cees van Miltenburg, Jacobus G. Blok,
Mark J. van Bommel, Paul J. van Ekeren, Gerrit J.K. van den Berg, and
Harry A.J. Oonk. "Liquid-vapour Equilibria of the Methyl Esters of
Alkanoic Acids: Vapour Pressures as a Function of Temperature and Standard
Thermodynamic Function Changes." Fluid Phase Equilibria 202, no. 1
(October 2002): 109-120. doi:10.1016/S0378-3812(02)00097-3.

The data were transcribed from the original references by Bryan W. Weber.
This file is distributed as part of the work "Experiments and Modeling of
the Autoignition of Methyl Valerate at Low to Intermediate Temperatures
and Elevated Pressures in a Rapid Compression Machine" by Bryan W. Weber,
Justin A. Bunnell, Kamal Kumar, and Chih-Jen Sung.

To the extent possible under law, the authors waive all copyright and related
or neighboring rights to this file. This work is licensed under the CC0 public
domain license, available at http://creativecommons.org/publicdomain/zero/1.0/
"""
import numpy as np
from scipy.optimize import curve_fit


def antoine(T, A, B, C):
    return A - B/(T - C)


ortega_data = np.loadtxt('ortega.csv', delimiter=',')
vanGenderen_data = np.loadtxt('vanGenderen.csv', delimiter=',')
verevkin_data = np.loadtxt('verevkin.csv', delimiter=',')
all_data = np.vstack((ortega_data, vanGenderen_data, verevkin_data))
all_data = all_data[all_data[:, 0].argsort()]
coeffs, cov_mat = curve_fit(antoine, all_data[:, 0], np.log10(all_data[:, 1]),
                            p0=[6, 1500, 50], maxfev=10000)
print(coeffs)
conf_int = 2*np.sqrt(np.diag(cov_mat))
print(conf_int)
