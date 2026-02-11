import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("epsilon_chi2_utf8.txt", skiprows=1)

eps  = data[:, 0]
chi2 = data[:, 1]

cl68, cl90, cl95 = 2.30, 4.61, 5.99

plt.plot(eps, chi2, marker="o", linewidth=2)
plt.axhline(cl68, linestyle="--", label="68% C.L.")
plt.axhline(cl90, linestyle="--", label="90% C.L.")
plt.axhline(cl95, linestyle="--", label="95% C.L.")

plt.xlabel(r"$\epsilon_R$")
plt.ylabel(r"$\Delta\chi^2_{\rm bias}$")
plt.legend(frameon=False)
plt.tight_layout()
plt.savefig("epsilon_chi2_bias.pdf")
plt.savefig("epsilon_chi2_bias.png", dpi=200)
