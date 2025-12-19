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
    
    theta_deg = np.degrees(theta)
    theta_degs = np.concatenate([theta_deg, theta_deg + 180])
    d_times_degs = np.concatenate([d_times, d_times[::-1]])
    d_db_degs = np.concatenate([d_dB, d_dB[::-1]])
    
    axs[0,0].plot(theta_degs, d_times_degs, color='blue')
    axs[0,0].set_title("КНД (разы, декарт)")
    axs[0,0].set_xlabel("θ (градусы)")
    axs[0,0].set_ylabel("D(θ)")
    axs[0,0].grid(True)
    axs[0,0].set_xlim(0, 360)
    
    axs[0,1].plot(theta_degs, d_db_degs, color='red')
    axs[0,1].set_title("КНД (дБ, декарт)")
    axs[0,1].set_xlabel("θ (градусы)")
    axs[0,1].set_ylabel("D(θ) [дБ]")
    axs[0,1].grid(True)
    axs[0,1].set_xlim(0, 360)
    
    axs[1,0] = plt.subplot(2,2,3, polar=True)
    theta_rads = np.concatenate([theta, theta + np.pi])
    d_times_rads = np.concatenate([d_times, d_times[::-1]])
    d_db_rads = np.concatenate([d_dB, d_dB[::-1]])
    axs[1,0].plot(theta_rads, d_times_rads, color='blue')
    axs[1,0].set_title("КНД (разы, поляр)")
    
    axs[1,1] = plt.subplot(2,2,4, polar=True)
    axs[1,1].plot(theta_rads, d_db_rads, color='red')
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
    
    d_times_0_180 = D(theta)
    d_db_0_180 = 10*np.log10(D(theta) + 1e-9)
    
    theta_0_360 = np.concatenate([theta, theta + np.pi])
    d_times_0_360 = np.concatenate([d_times_0_180, d_times_0_180[::-1]])
    d_db_0_360 = np.concatenate([d_db_0_180, d_db_0_180[::-1]])
    theta_deg_0_360 = np.degrees(theta_0_360)  
    
    with open('Резульаты_питон.txt', 'w', encoding='utf-8') as file:
        file.write('theta(град)   d_times   d_db\n')
        for i in range(len(theta_0_360)):
            angle_deg = theta_deg_0_360[i] % 360
            file.write(f'{angle_deg:.6f}   {d_times_0_360[i]:.6f}   {d_db_0_360[i]:.6f}\n')
    
    print("Результаты сохранены в файл 'Резульаты_питон.txt'")

    creating_plot(d_times=d_times_0_360, d_dB=d_db_0_360, theta=theta_0_360)

if __name__=="__main__":
    main()
