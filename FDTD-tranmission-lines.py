# //////////////////////////////////////////////////////////////////
# //  SEL0612 - Ondas eletromagneticas  (2020)                    //
# //  Trabalho 01 - Linhas de trasmissao                          //
# //                                                              //
# //  Alunos:                                                     //
# //  André Baconcelo Prado Furlanetti - Nº USP: 10748305         //
# //  Diego da Silva Parra             - Nº USP: 10716550         //
# //  Mateus Fernandes Doimo           - Nº USP: 10691971         //
# //////////////////////////////////////////////////////////////////

import numpy as np
import math
import matplotlib.pyplot as pyplot
import matplotlib.animation as animation
import scipy.constants as sci
import matplotlib as mpl
#mpl.rcParams['animation.ffmpeg_path'] = 'ffmpeg/bin/ffmpeg.exe' # Usado para salvar as animações

# Retona qual fonte o usuário escolheu
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

# Atualiza os valores da voltagem na animação
def animacao_voltagem(n):
    voltage_plot.set_ydata(V[n])
    return voltage_plot

# Atualiza os valores da corrente na animação
def animacao_corrente(n):
    current_plot.set_ydata(I[n])
    return current_plot

# Pega as informações do usuário 
# Fonte = 1 --> *u(t)   2 --> u(t) - u(t - l/(10*uf))
res_fonte = input("Fonte: 1 ou 2 ? ")
# Carga
res_carga = input("Carga: 0, 100 ou inf? ")

l = int(input("Tamanho da linha (l): "))
K = int(input("Numero de interacoes: "))

#Carga
if str(res_carga) == "0":
    Rl = 0
    nome_anim_RL = "0"
elif str(res_carga) == "100":
    Rl = 100
    nome_anim_RL = "100"
else:
    Rl = math.inf
    nome_anim_RL = "Inf"

#Velocidade do sinal de tensão ou de corrente
uf = 0.9*sci.speed_of_light
#Tempo total
tempo_t = 10*l/uf
#Tamanho de dz
dz = l/K
#Tamanho do dt
dt = dz/uf*0.2
# Número de iterações no tempo
N = int(tempo_t/dt)
#Impedância característica
Zo = 50
#Resistência interna da fonte de tensão
Rs = 75 
#Cálculo da capacitância a da indutância
C = 1/(Zo*uf)   #Capacitância
L = Zo/uf       #Indutância

#Matrizes de tensão e corrente
V = np.zeros((N, K))        #Tensão
I = np.zeros((N, K-1))      #Corrente

#Tamanho do eixo x dos gráficos
eixo_x = np.array(range(K))*dz

#Tensão e corrente com t --> ao infinito
V_inf = (Rl*fonte(res_fonte,N*dt,l))/(Rl+Rs)
I_inf = fonte(res_fonte,N*dt,l)/(Rl+Rs)

print("\nResultados:")
print("V∞(z,t) = "+str(V_inf))
print("I∞(z,t) = "+str(I_inf))


#Condições iniciais de tensão e corrente
V[0][0] = (Zo*fonte(res_fonte, 0,K*dz))/(Rs+Zo)
I[0][0] = fonte(res_fonte,0, K*dz)/(Rs+Zo)


beta_1 = 2*dt/(Rs*C*dz)

if Rl == 0:
    beta_2 = math.inf
elif Rl == math.inf:
    beta_2 = 0
else:
    beta_2 = 2*dt/(Rl*C*dz)

# Metodo FDTD
for n in range(1,N):
    # Condicao de controrno da tensao no inicio da linha
    V[n][0] = (1 - beta_1)*V[n-1][0] -2*I[n-1][0] + (2/Rs)*fonte(res_fonte,(n-1)*dt,l)
    
    # Tensao para toda linha K
    for k in range (1,K-1):
        V[n][k] = V[n-1][k] - (I[n-1][k]-I[n-1][k-1])

    # Condicao de contorno da tensao no final da linha
    if Rl == 0:
        V[n][K-1] = 0
    else:
        V[n][K-1] = (1 - beta_2)*V[n-1][K-1] + 2*I[n-1][K-2]

    # Corrente para toda linha K-1
    for k in range (0, K-1):
        I[n][k] = I[n-1][k] - (dt**2)/(L*C*dz**2)*(V[n][k+1] - V[n][k])

V = V*(dt/(C*dz))

################### Animações ###################

# Ajustes no tamanho do eixo y
if Rl == math.inf:
    if int(res_fonte) == 1:
        volt_limite_y = 3
    else:
        volt_limite_y = 2
    corrente_limite_y = 0.04
elif Rl == 0:
    volt_limite_y = 2.2
    if int(res_fonte) == 1:
        corrente_limite_y = 0.06
    else:
        corrente_limite_y = 0.03
else:
    volt_limite_y = 2.2
    corrente_limite_y = 0.04

# Configurando gráfico da tensão
figure, voltage = pyplot.subplots(1,1)
voltage.grid(True)
figure.patch.set_facecolor('#E0E0E0')
figure.patch.set_alpha(0.7)
voltage_plot, = voltage.plot(np.linspace(0,l,K), V[0], color='b', label='Voltage [V]')
voltage.set_ylim(-1.2, volt_limite_y)
voltage.legend()
voltage.set_ylabel('V(z,t)')
voltage.set_xlabel('z (m)')
voltage.set_title('Voltagem')


#Writer = animation.writers['ffmpeg']
#writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

animation1 = animation.FuncAnimation(figure, func = animacao_voltagem, frames=np.arange(0, N, (int)(K/10)), interval = 100, repeat = False)

# Para salvar a animação em .mp4
#animation1.save("Fonte_"+res_fonte+"_voltagem_RL_"+nome_anim_RL+".mp4", writer = 'ffmpeg')

# Configurando gráfico da corrente
figure2, current = pyplot.subplots(1,1)
current.grid(True)
figure2.patch.set_facecolor('#E0E0E0')
figure2.patch.set_alpha(0.7)
current_plot, = current.plot(np.linspace(0,l,K-1), I[0], color='g', label='Corrente [A]')
current.set_ylim(-0.025, corrente_limite_y)
current.legend()
current.set_ylabel('I(z,t)')
current.set_xlabel('z (m)')
current.set_title('Corrente')

animation2 = animation.FuncAnimation(figure2, func = animacao_corrente, frames=np.arange(0, N, (int)(K/10)), interval = 100, repeat = False)

# Para salvar a animação em .mp4
#animation2.save("Fonte_"+res_fonte+"_corrente_RL_"+nome_anim_RL+".mp4", writer = 'ffmpeg')

#Mostra as animações
pyplot.show()