#' # DN2: Porazdelitvena funkcija za N(0,1) in Gauss-Legendrove kvadrature
#' Avtor: Patrik Schembri

#' ## Naloge s funkcijami: Porazdelitvena funkcija za N(0, 1)

#' V prvem delu naloge se ukvarjamo s pisanjem učinkovite funkcije za izračun vrednosti
#' porazdelitvene funkcije standardne normalne porazdelitve $N(0, 1)$.

#' Računamo torej:
#' $$
#' F(x) = \frac{1}{\sqrt{2\pi}} \int_{-\infty}^{x} e^{-\frac{t^2}{2}} dt
#' $$

using Domaca02, Distributions, Plots

#' ### Metoda
#' 
#' Funkcija `normalPor(x)` uporablja različne pristope glede na območje:
#' - za $x < 0$ izkorišča simetrijo $F(x) = 1 - F(-x)$,
#' - za $x \in [-10, 10]$ uporabi sestavljeno Simpsonovo metodo za integracijo,
#' - za $x > 10$ uporabi repno integracijo $Q(x) = \int_x^{\infty} f(t)\,dt$,
#'   kar omogoča, da relativna napaka ostane majhna tudi pri zelo velikih vrednostih $x$.
#'
#' Tako zagotovimo, da je relativna natančnost < $10^{-10}$ na celem definicijskem območju.

#' Za test najprej uporabimo funkcijo $f(x) = x$ in izračunamo vrednost integrala s sestavljeno
#' Simpsonovo metodo na intervalu [0,1].
f(x) = x
a1, b1 = 0, 1
simpson = simpsonovoPravilo(f, a1, b1, 100)

#' Nato izračunamo vrednosti porazdelitvene funkcije standardne normalne porazdelitve
#' s funkcijo `normalPor(x)` in rezultate primerjamo z vrednostmi knjižnice `Distributions.jl`.
dist = Normal(0, 1)
xs = [-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0]
razlike = Float64[]

for x in xs
    res = normalPor(x)
    pricakovano = cdf(dist, x)  # Točna vrednost iz Distributions.jl
    razlika = abs(res - pricakovano)
    push!(razlike, razlika)
end

plot(xs, razlike, label="Razlika", xlabel="x", 
    ylabel="Razlika", title="Razlike med točnim in izračunanim rezultatom")


#' ## Naloge s števili: Gauss-Legendrove kvadrature

#' V drugem delu naloge se ukvarjamo z izpeljavo Gauss-Legendreovega pravila za numerično integracijo.
#' Najprej ga izračunamo za dve točki po formuli:
#' $$
#' \int_a^b f(x) = A(f(x_1)) + B(f(x_2)) + R_f
#' $$
#' Za primer izračunajmo $\int_1^2 \sin((x))dx$ s pomočjo tega pravila.
f(x) = sin(x)
a, b = 1.0, 2.0
gaussLegendre2P(f, a, b)

#' Sedaj lahko izpeljemo še sestavljeno Gauss-Legendreovo pravilo, kjer moramo določiti
#' še število intervalov na katerega bo glavni interval [a, b] razdeljen. Izračun lahko
#' naredimo na enaki funkciji kot v prvem primeru, a interval razdelino na 10 delov.
n = 10
gaussLegendre(f, a, b, n)

#' Sestavljeno pravilo lahko sedaj uporabimo za izračun različnih funkcij. Za še en
#' primer lahko vzamemo izračun integrala:
#' $$
#' \int_0^5 \frac{sin(x)}{x}dx
#' $$
#' Ocenimo še koliko izračunov funkcijske vrednosti je potrebnih, da je izračun približka
#' integrala natančen na 10 decimalk.
rez, stKorakov = primerGL()

#' Dobimo rezultat:
rez

#' Dobimo število potrebnih korakov da dosežemo željeno natančnost:
stKorakov