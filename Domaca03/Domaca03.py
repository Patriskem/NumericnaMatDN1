import math
import matplotlib.pyplot as plt
from fractions import Fraction

def solve_dopri5(f, u0, tspan, h, p=None, tol=1e-11):
    """
    DOPRI5 metoda (Dormand-Prince 5. reda) za reševanje sistema diferencialnih enačb.
    f(u, p, t) - funkcija desne strani
    u0 - začetni pogoji (seznam ali tuple)
    tspan - (t0, T)
    h - začetni korak
    p - dodatni parametri
    tol - relativna natančnost
    Vrne seznam časov ts in rešitev us (seznam vektorjev).
    """

    # Koeficienti DOPRI5 (Dormand-Prince)
    c = [Fraction(0,1), Fraction(1,5), Fraction(3,10), Fraction(4,5), 
     Fraction(8,9), Fraction(1,1), Fraction(1,1)]
    a = [
    [],
    [Fraction(1,5)],
    [Fraction(3,40), Fraction(9,40)],
    [Fraction(44,45), Fraction(-56,15), Fraction(32,9)],
    [Fraction(19372,6561), Fraction(-25360,2187), Fraction(64448,6561), Fraction(-212,729)],
    [Fraction(9017,3168), Fraction(-355,33), Fraction(46732,5247), Fraction(49,176), Fraction(-5103,18656)],
    [Fraction(35,384), Fraction(0,1), Fraction(500,1113), Fraction(125,192), Fraction(-2187,6784), Fraction(11,84)]
        ]
    b5 = [Fraction(35,384), 0, Fraction(500,1113), Fraction(125,192), Fraction(-2187,6784), Fraction(11,84), 0]
    b4 = [Fraction(5179,57600), 0, Fraction(7571,16695), Fraction(393,640), Fraction(-92097,339200), Fraction(187,2100), Fraction(1,40)]

    t0, T = tspan
    ts = [t0]
    us = [list(u0)]
    t = t0
    u = list(u0)

    while t < T - 1e-15:
        if t + h > T:
            h = T - t

        k = []
        for i in range(7):
            u_stage = [u[j] + h * float(sum(a[i][m] * k[m][j] for m in range(len(a[i])))) for j in range(len(u))]
            t_stage = t + c[i] * h
            k.append(f(u_stage, p, t_stage))

        u5 = [u[j] + h * float(sum(b5[m] * k[m][j] for m in range(7))) for j in range(len(u))]
        u4 = [u[j] + h * float(sum(b4[m] * k[m][j] for m in range(7))) for j in range(len(u))]

        # napaka
        err = max(abs(u5[j] - u4[j]) for j in range(len(u)))

        if err < tol:
            t += h
            u = u5
            ts.append(t)
            us.append(u)

        # prilagoditev koraka
        if err == 0:
            s = 2
        else:
            s = 0.9 * (tol / err) ** 0.2
        h = h * min(5, max(0.1, s))

    return ts, us


def nihalo(l, t, theta0, omega0, h, g=9.80665):
    """
    Reši matematično nihalo:
    θ''(t) + g/l * sin(θ(t)) = 0
    začetni pogoji: θ(0) = theta0, θ'(0) = omega0
    Vrne θ(t) ob času t.
    """
    def f(u, p, tau):
        theta, omega = u
        return [omega, -(g/l) * math.sin(theta)]

    ts, us = solve_dopri5(f, [theta0, omega0], (0.0, t), h)
    return us[-1][0]


def harm_nahalo(l, t, theta0, omega0, h, g=9.80665):
    """
    Reši harmonično nihalo (linearizirano):
    θ''(t) + g/l * θ(t) = 0
    začetni pogoji: θ(0) = theta0, θ'(0) = omega0
    Vrne θ(t) ob času t.
    """
    def f(u, p, tau):
        theta, omega = u
        return [omega, -(g/l) * theta]

    ts, us = solve_dopri5(f, [theta0, omega0], (0.0, t), h)
    return us[-1][0]


def narisi_energijo_proti_casu(theta0_range, l=1.0, g=9.80665, h=1e-3):
    """
    Nariše graf odvisnosti energije nihala od nihajnega časa (E(T)).

    PARAMETRI:
    theta0_range : seznam kotov v radianih
        Začetni koti nihala (v radianih), ki določajo začetno energijo.
    l : float
        Dolžina nihala v metrih.
    g : float
        Gravitacijski pospešek (privzeto 9.80665 m/s²).
    h : float
        Začetna velikost koraka za integracijo (za DOPRI5 metodo).

    OPIS:
    Matematično nihalo sledi enačbi gibanja:
        θ''(t) + (g/l) * sin(θ(t)) = 0

    Energija nihala (brez dušenja, brez zunanjih sil) se ohranja in je enaka:
        E = 1/2 * l² * ω² + g * l * (1 - cos(θ))
      kjer:
        - ω je kotna hitrost (θ'),
        - θ je kot glede na navpičnico.
    """
    energije = []
    casi = []

    for theta0 in theta0_range:
        omega0 = 0.0  # začetna hitrost

        def f(u, p, t):
            theta, omega = u
            return [omega, -(g / l) * math.sin(theta)]

        u0 = [theta0, omega0]
        t_max = 20  # maksimalni čas simulacije

        ts, us = solve_dopri5(f, u0, (0.0, t_max), h)

        # Poiščemo čas, ko theta prvič prečka 0 z negativno hitrostjo
        crossing_time = None
        for i in range(1, len(us)):
            if us[i-1][0] > 0 and us[i][0] <= 0 and us[i][1] < 0:
                t1, t2 = ts[i-1], ts[i]
                theta1, theta2 = us[i-1][0], us[i][0]
                frac = theta1 / (theta1 - theta2)
                crossing_time = t1 + frac * (t2 - t1)
                break

        if crossing_time:
            T = 2 * crossing_time  # celoten nihajni čas (simetrija)
            E = 0.5 * (l**2) * omega0**2 + g * l * (1 - math.cos(theta0))
            casi.append(T)
            energije.append(E)

    plt.figure(figsize=(8, 5))
    plt.plot(energije, casi, marker='o')
    plt.xlabel('Energija [J]')
    plt.ylabel('Nihajni čas [s]')
    plt.title('Odvisnost nihajnega časa nihala od energije')
    plt.grid(True)
    plt.show()


def primerjaj_nihali(theta0, l=1.0, g=9.80665, h=1e-3, t_max=10):
    """
    Nariše graf odmika matematičnega in harmoničnega nihala za podane začetne pogoje.
    """
    omega0 = 0.0

    # Matematično nihalo
    def f_math(u, p, t):
        theta, omega = u
        return [omega, -(g/l) * math.sin(theta)]

    ts_math, us_math = solve_dopri5(f_math, [theta0, omega0], (0.0, t_max), h)

    # Harmonično nihalo
    def f_harm(u, p, t):
        theta, omega = u
        return [omega, -(g/l) * theta]

    ts_harm, us_harm = solve_dopri5(f_harm, [theta0, omega0], (0.0, t_max), h)

    # Risanje
    plt.figure(figsize=(8, 5))
    plt.plot(ts_math, [u[0] for u in us_math], label="Matematično nihalo", color="blue")
    plt.plot(ts_harm, [u[0] for u in us_harm], label="Harmonično nihalo", linestyle="--", color="red")
    plt.xlabel("Čas [s]")
    plt.ylabel("Kot θ(t) [rad]")
    plt.title(f"Primerjava nihanja (θ0={theta0:.2f} rad)")
    plt.legend()
    plt.grid(True)
    plt.show()


__all__ = ["solve_dopri5", "nihalo", "narisi_energijo_proti_casu", "harm_nahalo", "primerjaj_nihali"]
