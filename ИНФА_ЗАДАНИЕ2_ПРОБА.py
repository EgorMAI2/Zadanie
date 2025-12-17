
import numpy as np
import matplotlib.pyplot as plt

def E(theta):
    num = np.cos(k * l * np.cos(theta)) - np.cos(k * l)
    den = np.sin(theta)
    return num/den

def F(theta):
    return abs(E(theta)) / abs(E(theta).max())

def Dmax(theta):
    formula = (F(theta)**2 * np.sin(theta))
    return 2 / np.trapezoid(formula, theta)

def D(theta):
    return F(theta)**2 * Dmax(theta)

def creating_plot(d_times, d_dB, theta):
    fig, axs = plt.subplots(2, 2, figsize=(12,10), subplot_kw={'polar': False})
    fig.suptitle('D(Theta)')
    
    axs[0,0].plot(theta, d_times, color='blue')
    axs[0,0].set_title("КНД (разы, декарт)")
    axs[0,0].set_xlabel("θ (рад)")
    axs[0,0].set_ylabel("D(θ)")
    axs[0,0].grid(True)
    
    axs[0,1].plot(theta, d_dB, color='red')
    axs[0,1].set_title("КНД (дБ, декарт)")
    axs[0,1].set_xlabel("θ (рад)")
    axs[0,1].set_ylabel("D(θ) [дБ]")
    axs[0,1].grid(True)
    
    axs[1,0] = plt.subplot(2,2,3, polar=True)
    axs[1,0].plot(theta, d_times, color='blue')
    axs[1,0].set_title("КНД (разы, поляр)")
    
    axs[1,1] = plt.subplot(2,2,4, polar=True)
    axs[1,1].plot(theta, d_dB, color='red')
    axs[1,1].set_title("КНД (дБ, поляр)")
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig('График_задание_2.png')
    print("График сохранен как 'График_задание_2.png'")
    plt.show()

def main():
    global l, k
    f = 0.3 * 10 ** 9
    lmbd = 3 * 10 ** 8 / f
    l = 0.3  
    k = 2 * np.pi / lmbd
    
    theta = np.linspace(1e-9, np.pi-(1e-9), 2000)
    
    print(f'{Dmax(theta=theta):.3f} times\n{10 * np.log10(Dmax(theta=theta)):.3f} dB')
    
    d_times = D(theta)
    d_db = 10*np.log10(D(theta) + 1e-9)
    
    # Сохранение результатов
    with open('Резульаты_питон.txt', 'w', encoding='utf-8') as file:
        file.write('theta   d_times   d_db\n')
        for i in range(len(theta)):
            file.write(f'{theta[i]:.6f}   {d_times[i]:.6f}   {d_db[i]:.6f}\n')
    
    print("Результаты сохранены в файл 'Резульаты_питон.txt'")
    
    # Построение графика
    creating_plot(d_times=d_times, d_dB=d_db, theta=theta)

if __name__=="__main__":
    main()
