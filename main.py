import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm_notebook

# Aceleración
grav = 9.81

# Masa del dipolo
m = 0.01

# es la permeabilidad máxima del MetGlass
mu = 1e6

# Permeabilidad del espacio libre
mu_0 = 4 * np.pi * 1e-7

# Radio del anillo de corriente
a = 0.08

# Resistividad
R = 9e-5

# Calculo constante k
# $$k = \frac{9(\mu\mu_0)^2a^4}{4R}$$
k = (9 * (mu * mu_0) ** 2 * a**4) / (4 * R)


def p(x):
    return x**2 / ((x**2 + a**2) ** (5 / 2))


def dp_dx(x):
    return (10 * x**2) / ((x**2 + a**2) ** (7 / 2)) - (x**2) / (
        (x**2 + a**2) ** (5 / 2)
    )


def f(t, x, y):  # ecuacion de la posicion
    return y


def g(t, x, y):  # ecuacion de la velocidad
    return -grav - (k / m) * p(x) * y


def dg_dt(
    t, x, y
):  # ecuacion de la derivada de la velocidad # ecuacion de la aceleracion
    return (dp_dx(x) * -(y**2) - (k / m) * p(x) * y) - 10


# Tiempo inicial
t_0 = 0
# Tiempo final
t_f = 6

# Altura Inicial x_0
x_0 = 10
# Velocidad Inicial y_0
y_0 = 0

# Paso horizontal
h = 0.01

# Número de pasos
# Calcula el número de pasos entre t_0 y t_f cuando tienes un paso horizontal de h
N = int((t_f - t_0) / h)

# Arreglo con valores del tiempo:
# Usa np.linspace() para definir un arreglo que vaya de t_0 a t_f en N pasos
t = np.linspace(t_0, t_f, N)

all_xs = []
all_ys = []
all_as = []  # List to store acceleration values

# Usar los valores iniciales para el primer paso
x_n = x_0
y_n = y_0

# Iterar para cada paso de la lista de tiempo t
# Usamos tqdm para generar una barra de progreso
for t_n in tqdm_notebook(t):
    # Para cada iteración debes calcular x_nplus1 y y_nplus1
    kx1 = f(t_n, x_n, y_n)
    ky1 = g(t_n, x_n, y_n)

    kx2 = f(t_n + 0.5 * h, x_n + 0.5 * h * kx1, y_n + 0.5 * h * ky1)
    ky2 = g(t_n + 0.5 * h, x_n + 0.5 * h * kx1, y_n + 0.5 * h * ky1)

    kx3 = f(t_n + 0.5 * h, x_n + 0.5 * h * kx2, y_n + 0.5 * h * ky2)
    ky3 = g(t_n + 0.5 * h, x_n + 0.5 * h * kx2, y_n + 0.5 * h * ky2)

    kx4 = f(t_n + h, x_n + h * kx3, y_n + h * ky3)
    ky4 = g(t_n + h, x_n + h * kx3, y_n + h * ky3)

    x_nplus1 = x_n + (1 / 6) * h * (kx1 + 2 * kx2 + 2 * kx3 + kx4)
    y_nplus1 = y_n + (1 / 6) * h * (ky1 + 2 * ky2 + 2 * ky3 + ky4)

    # No olvides guardar cada valor x_nplus1 y y_nplus1 calculado
    # en las listas all_xs y all_ys
    all_xs.append(x_nplus1)
    all_ys.append(y_nplus1)

    # Tampoco olvides actualizar x_n y y_n para calcular algo nuevo
    # en cada iteración
    x_n = x_nplus1
    y_n = y_nplus1

    # Calculate acceleration using the current velocity
    a_n = dg_dt(t_n, x_n, y_n)

    # Store the calculated acceleration in all_as list
    all_as.append(a_n)


# plt.figure(figsize=(16, 4))
# subplot de 4
plt.figure(figsize=(16, 4))

plt.subplot(1, 3, 1)
plt.plot(t, all_xs)
plt.title("Altura de Dipolo con Frenado Magnético")
plt.xlabel("Tiempo [s]")
plt.ylabel("Altura [m]")

plt.subplot(1, 3, 2)
plt.plot(t, all_ys)
plt.title("Velocidad de Dipolo con Frenado Magnético")
plt.xlabel("Tiempo [s]")
plt.ylabel("Velocidad [m/s]")

plt.subplot(1, 3, 3)
plt.plot(t, all_as)
plt.title("Aceleración de Dipolo con Frenado Magnético")
plt.xlabel("Tiempo [s]")
plt.ylabel("Aceleración [m/s^2]")


plt.tight_layout()
plt.show()
