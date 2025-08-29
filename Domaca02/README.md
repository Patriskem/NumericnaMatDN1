# Domača Naloga 2 - Porazdelitvena funkcija za N(0,1) in Gauss-Legendrove kvadrature

##### Avtor: Patrik Schembri

## Opis naloge
### V prvem delu naloge se ukvarjamo s pisanjem učinkovite funkcije za izračun vrednosti porazdelitvene funkcije standardne normalne porazdelitve N(0, 1).
Funkcija normalPor(x) uporablja različne pristope glede na območje.

### Naloge s števili: Gauss-Legendrove kvadrature
V drugem delu naloge se ukvarjamo z izpeljavo Gauss-Legendreovega pravila za numerično integracijo. 

V prvem koraku sem zapisal funkcijo za izračun Gauss-Legendrove kvadrature na dveh točkah. Glavni del te naloge je izpeljati sestavljeno pravilo za integracijo. To je naloga funkcije ```gaussLegendre(f, a, b, n)```. Na koncu pa izračunam še približek za dejanski primer integrala in število izračunov, ki jih moramo opraviti da dosežemo željeno natančnost.

## Uporaba kode
1. V terminalu zaženemo Julio
2. Pritisnemo "]" in izvedemo naslednje ukaze:
    1. ```activate Domaca02```
    2. ```add Distributions```
    3. ```add Plots```
    4. ```add Test```
    5. ```add QuadGK```
3. Nato lahko v Julia terminalu poganjamo posamezne vrstice kode (v Visual Studio Code s stiskom na Shift+Enter na vrstici ki jo želimo izvesti).

## Poganjanje testov
1. V terminalu zaženemo Julio
2. Pritisnemo "]" in izvedemo naslednje ukaze:
    1. ```activate Domaca02```
    2. ```test```
3. Testi se izvedejo

## Generiranje poročila
(0. Porocilo demo.pdf je ze ustvarjeno v mapi scripts)
1. V terminalu zaženemo Julio
2. Pritisnemo "]" in izvedemo: 
    1. ```activate Domaca02```
    2. ```add Weave```
3. Pritisnemo "Backspace", da gremo iz pkg
4. Izvedemo naslednje ukaze:
    1. ```cd("Domaca02\\scripts")```
    2. ```include("makedocs.jl")```
5. Ustvari se poročilo z imenom "demo.pdf"