
import numpy as np
from scipy.special import spherical_jn, spherical_yn
import matplotlib.pyplot as plt
import json

class EPRCalculator:
    def __init__(self, D: float, c=3e8, N: int = 50):
        self.c = c
        self.N = N
        self.r = D / 2
    
    def hankel(self, n, k, r):
        return spherical_jn(n, k * r) + 1j * spherical_yn(n, k * r)
    
    def an(self, n, k, r):
        return spherical_jn(n, k * r) / self.hankel(n, k, r)
    
    def bn(self, n, k, r):
        if n == 1:
            jn_minus1 = spherical_jn(0, k * r)
            hn_minus1 = spherical_jn(0, k * r) + 1j * spherical_yn(0, k * r)
        else:
            jn_minus1 = spherical_jn(n - 1, k * r)
            hn_minus1 = self.hankel(n - 1, k, r)
        
        jn = spherical_jn(n, k * r)
        hn = self.hankel(n, k, r)
        
        numerator = k * r * jn_minus1 - n * jn
        denominator = k * r * hn_minus1 - n * hn
        
        return numerator / denominator
    
    def sigma(self, lmbd, k, r):
        s = 0 + 0j
        for n in range(1, self.N + 1):
            term = ((-1) ** n) * (n + 0.5) * (self.bn(n, k, r) - self.an(n, k, r))
            s += term
        return (lmbd ** 2 / np.pi) * abs(s) ** 2
    
    def calculate(self, frequencies):
        sigmas = []
        for f in frequencies:
            lmbd = self.c / f
            k = 2 * np.pi / lmbd
            sigmas.append(self.sigma(lmbd, k, self.r))
        return sigmas

class Plotter:
    def save_to_json(frequencies, sigmas, D, filename="Результаты.json"):
        data_list = []
        for f, s in zip(frequencies, sigmas):
            entry = {"freq": float(f), "lambda": float(3e8/f), "rcs": float(s)}
            data_list.append(entry)
        
        result = {"data": data_list}
        with open(filename, 'w', encoding='utf-8') as f:
            json.dump(result, f, indent=4)
    
    def plot_results(frequencies, sigmas, D, filename="график.png"):
        plt.figure(figsize=(10, 6))
        frequencies = np.array(frequencies) / 1e9
        plt.plot(frequencies, sigmas, 'b-', linewidth=2)
        plt.xlabel('Частота, ГГц')
        plt.ylabel('ЭПР, м²')
        plt.title(f'Зависимость ЭПР от частоты\nD = {D*1000:.0f} мм')
        plt.grid(True)
        plt.savefig(filename, dpi=150)
        plt.show()

def main():
    D = 30e-3
    fmin = 0.01e9
    fmax = 25e9
    
    frequencies = np.linspace(fmin, fmax, 200)
    calculator = EPRCalculator(D=D, N=50)
    sigmas = calculator.calculate(frequencies)
    
    Plotter.save_to_json(frequencies, sigmas, D)
    Plotter.plot_results(frequencies, sigmas, D)

if __name__ == '__main__':
    main()
