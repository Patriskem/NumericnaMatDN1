import math
from Domaca03 import solve_dopri5, nihalo, narisi_energijo_proti_casu, primerjaj_nihali, harm_nahalo

# funkcija za primerjanje vrednosti z zazeljeno relativno tocnostjo
def almost_equal(a, b, tol=1e-10):
    return abs(a - b) < tol * (1 + abs(b))


# tu imamo referenco e^T, ki jo lahko izračunamo z math.exp z visoko natančnostjo
def test_exp_growth():
    print("Test: eksponentna rast y' = y, y(0)=1")
    f = lambda u, p, t: [u[0]]
    for T in [0.5, 1.0, 2.0, 5.0]:
        ts, us = solve_dopri5(f, [1.0], (0.0, T), 0.1)
        approx = us[-1][0]
        exact = math.exp(T)
        if almost_equal(approx, exact):
            print(f"  ✅ PASS T={T}: approx={approx:.12f}, exact={exact:.12f}")
        else:
            print(f"  ❌ FAIL T={T}: approx={approx:.12f}, exact={exact:.12f}")


# test periode
def test_nihajoce():
    print("Test: harmonično nihalo (majhen kot)")
    l = 1.0
    g = 9.80665
    T_period = 2 * math.pi * math.sqrt(l/g)
    theta0 = 0.01
    omega0 = 0.0

    for k in range(1, 4):  # 1,2,3 periode
        val = nihalo(l, k*T_period, theta0, omega0, h=1e-2)
        if almost_equal(val, theta0):
            print(f"  ✅ PASS {k} period: theta≈{val:.12f}")
        else:
            print(f"  ❌ FAIL {k} period: theta≈{val:.12f}, expected≈{theta0:.12f}")


# sanity check da numerična rešitev ne "eksplodira"
def test_vecje_amplitude():
    print("Test: matematično nihalo (večji odmik)")
    l = 1.0
    theta0 = math.pi/2
    omega0 = 0.0

    times = [0.5, 1.0, 2.0, 3.0]
    for t in times:
        val = nihalo(l, t, theta0, omega0, h=1e-3)
        if abs(val) < math.pi:
            print(f"  ✅ PASS t={t}: theta≈{val:.12f}")
        else:
            print(f"  ❌ FAIL t={t}: theta≈{val:.12f}, out of bounds")


# narisemo graf kota matematičnega in harmoničnega nihala
def test_primerjava():
    print("1. primerjava matematičnega in harmoničnega nihala")
    theta0 = 1  # 60 stopinj priblizno
    l = 1.0
    t = 2.0
    val_math = nihalo(l, t, theta0, 0.0, h=1e-3)
    val_harm = harm_nahalo(l, t, theta0, 0.0, h=1e-3)
    print(f"  Matematično nihalo θ(t)≈{val_math:.12f}")
    print(f"  Harmonično nihalo θ(t)≈{val_harm:.12f}")

    # nariši primerjalni graf
    print("...Risanje matematičnega in harmoničnega nihala...")
    primerjaj_nihali(theta0, l=l)

    print("2. primerjava matematičnega in harmoničnega nihala")
    theta0 = 0.34  # 20 stopinj priblizno
    l = 1.0
    t = 2.0
    val_math = nihalo(l, t, theta0, 0.0, h=1e-3)
    val_harm = harm_nahalo(l, t, theta0, 0.0, h=1e-3)
    print(f"  Matematično nihalo θ(t)≈{val_math:.12f}")
    print(f"  Harmonično nihalo θ(t)≈{val_harm:.12f}")

    # nariši primerjalni graf
    print("...Risanje matematičnega in harmoničnega nihala...")
    primerjaj_nihali(theta0, l=l)

    
# narisemo graf nihajnega casa v odvisnosti od energije
def graf():
    print("...Risanje grafa nihajnega casa v odvisnosti od energije...")
    theta0_range = [i * math.pi / 180 for i in range(5, 180, 5)]  # od 5° do 175°
    narisi_energijo_proti_casu(theta0_range, l=1.0)

if __name__ == "__main__":
    test_exp_growth()
    test_nihajoce()
    test_vecje_amplitude()
    test_primerjava()
    graf()
