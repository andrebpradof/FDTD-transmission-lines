import numpy as np
import math
import matplotlib.pyplot as pyplot
import matplotlib.animation as animation
import scipy.constants as sci
import matplotlib as mpl
mpl.rcParams['animation.ffmpeg_path'] = r'C:\Users\Andre Prado\Documents\ffmpeg\bin\ffmpeg.exe'

def fonte(res_fonte,t,l):
    if int(res_fonte) == 1:
        return Vs1(t)
    else:
        return Vs2(t,l)

#Função Degrau Unitário
def u(t):
    if t >= 0:
        return 1
    else:
        return 0

#Primeira fonte de tensão
def Vs1(t):
    return 2*u(t)

#Segunda fonte de tensão
def Vs2(t,l):
    return u(t) - u(t - l/(10*uf))

def animate(n):
    voltage_plot.set_ydata(V[n])
    current_plot.set_ydata(I[n])
    return voltage_plot,

res_fonte = input("Fonte: ")
res_carga = input("Carga: ")

l = int(input("Tamanho da linha (l): "))
K = int(input("Numero de interacoes: "))

#Carga
if int(res_carga) == 0:
    Rl = 0
elif int(res_carga) == 100:
    Rl = 100
else:
    Rl = math.inf

#Velocidade do sinal de tensão ou de corrente
uf = 0.9*sci.speed_of_light
tempo_t = 10*l/uf #3.70E-5
dz = l/K #Respeita a condição de estabilidade
dt = dz/uf*0.25
N = int(tempo_t/dt)
#Impedância característica
Zo = 50

#Resistência interna da fonte de tensão
Rs = 75 


#Cálculo da capacitância a da indutância
C = 1/(Zo*uf)
L = Zo/uf

#Matrizes de tensão e corrente
V = np.zeros((N, K))       #Tensão
I = np.zeros((N, K))     #Corrente

#Tamanho do eixo x dos gráficos
eixo_x = np.array(range(K))*dz

#Tensão e corrente com t --> ao infinito
V_inf = (Rl*fonte(res_fonte,N*dt,l))/(Rl+Rs)
I_inf = fonte(res_fonte,N*dt,l)/(Rl+Rs)

print("\nResultados:")
#print("V∞(z,t) = "+V_inf)
#print("I∞(z,t) = "+I_inf)


#Condições iniciais de tensão e corrente
V[0][0] = Zo*2/(Rs+Zo)
I[0][0] = V[0][0]/Zo

beta_S = 2*dt/(Rs*C*dz)
r = dt*dt/(L*C*dz*dz)

if Rl == 0:
    beta_L = math.inf
elif Rl == math.inf:
    beta_L = 0
else:
    beta_L = 2*dt/(Rl*C*dz)

for n in range(1,N):
    V[n][0] = (1 - beta_S)*V[n-1][0] -2*I[n-1][0] + (2/Rs)*fonte(res_fonte,n*dt,l)
    
    for k in range (1,K-1):
        V[n][k] = V[n-1][k] - (I[n-1][k]-I[n-1][k-1])

    if Rl == 0:
        V[n][K-1] = 0
    else:
        V[n][K-1] = (1 - beta_L)*V[n-1][K-1] + 2*I[n-1][K-2]

    for k in range (0, K-1):
        I[n][k] = I[n-1][k] - (dt**2)/(L*C*dz**2)*(V[n][k+1] - V[n][k])

V = V*(dt/(C*dz))


figure, (voltage) = pyplot.subplots(1,1)
voltage.grid(True)

voltage_plot, = voltage.plot(np.linspace(0,l,K), V[0], color='b', label='Voltage [V]')
voltage.set_ylim(-1, 3)
voltage.legend()

#Writer = animation.writers['ffmpeg']
#writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

animation1 = animation.FuncAnimation(figure, func = animate, frames=np.arange(0, N, (int)(K/10)), interval = 100, repeat = False)

#animation.save('animated_coil.mp4', writer = 'ffmpeg', fps = 30)

figure2, (current) = pyplot.subplots(1,1)
current.grid(True)

current_plot, = current.plot(np.linspace(0,l,K), I[0], color='g', label='Corrente [A]')
current.set_ylim(-0.05, 0.05)
current.legend()

#Writer = animation.writers['ffmpeg']
#writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

animation2 = animation.FuncAnimation(figure2, func = animate, frames=np.arange(0, N, (int)(K/10)), interval = 100, repeat = False)
pyplot.show()