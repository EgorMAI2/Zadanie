import matplotlib.pyplot as plt
import numpy as np

def results_from_py():
    with open('Резульаты_питон.txt', 'r', encoding='utf-8') as file:
        file.readline()
        
        axis =[[],[],[]]

        for line in file:
            theta, d_times, d_db = map(float, line.split())
            axis[0].append(theta)
            axis[1].append(d_db)
            axis[2].append(d_times)
        
        return axis

def results_from_CST():
    with open('CST_результаты.txt', 'r', encoding='utf-8') as file:
        axis = [[],[],[]]
        for line in file:
            parts = line.split()
            if len(parts) >= 3:
                try:
                    #1 столбец 
                    theta_deg = float(parts[0])
                    #3 столбец
                    d_db_str = parts[2].replace(',', '.')
                    d_db = float(d_db_str)
                    
                    axis[0].append(theta_deg)  
                    axis[1].append(d_db)
                    axis[2].append(10**(d_db/10))
                except:
                    continue
        return axis

def creating_plot(cst, python):
        
    fig = plt.figure(figsize=(12, 10))
    fig.suptitle('Разница между\nPython и CST')

    ax1 = fig.add_subplot(2, 2, 1)  
    ax2 = fig.add_subplot(2, 2, 2)  

    ax1.plot(cst[0], cst[1], label='CST')
    ax1.plot(python[0], python[1], label='Python')
    ax1.set_title('D, dBi')
    ax1.set_xlabel("Theta, градусы")
    ax1.set_ylabel("dBi")
    ax1.set_xlim(0, 360)  
    ax1.grid(True)
    ax1.legend()

    ax2.plot(cst[0], cst[2], label='CST')
    ax2.plot(python[0], python[2], label='Python')
    ax2.set_title('D, разы')
    ax2.set_xlabel("Theta, градусы")
    ax2.set_ylabel("разы")
    ax2.set_xlim(0, 360) 
    ax2.grid(True)
    ax2.legend()

    ax3 = fig.add_subplot(2, 2, 3, polar=True)
    ax4 = fig.add_subplot(2, 2, 4, polar=True)

    cst_rad = [np.deg2rad(angle) for angle in cst[0]]
    python_rad = [np.deg2rad(angle) for angle in python[0]]

    ax3.plot(cst_rad, cst[1], label='CST')
    ax3.plot(python_rad, python[1], label='Python')
    ax3.set_title('D, dBi (полярные координаты)')
    ax3.grid(True)
    ax3.legend()

    ax4.plot(cst_rad, cst[2], label='CST')
    ax4.plot(python_rad, python[2], label='Python')
    ax4.set_title('D, разы (полярные координаты)')
    ax4.grid(True)
    ax4.legend()

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig("График_задание2.png")
    plt.show()

def main():
    creating_plot(cst=results_from_CST(), python=results_from_py())

if __name__=='__main__':
    main()
