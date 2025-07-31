# Domača Naloga 1 - SOR iteracija za razpršene matrike

##### Avtor: Patrik Schembri

## Opis naloge
Naj bo A n × n diagonalno dominantna razpršena matrika(velika večina elementov je ničelnih aij=0).

Definiral sem nov podatkovni tip RazprsenaMatrika, ki matriko zaradi prostorskih zahtev hrani v dveh matrikah V in I, kjer sta V in I matriki n×m, tako da velja: ```V(i,j)=A(i,I(i,j))```
V matriki V se torej nahajajo neničelni elementi matrike A. Vsaka vrstica matrike V vsebuje ničelne elemente iz iste vrstice v A. V matriki I pa so shranjeni indeksi stolpcev teh neničelnih elementov.

Za ta podatkovni tip sem defirniral funkcije ```getIndex, setIndex, firstIndex, lastIndex``` in množenje z desne z vektorjem ```*```.

Nato sem definiral funkcijo SOR, ki reši SOR iteracijo. Ta za parametre sprejme razpršeno matriko A, vektor b, začetni približek x0 in relaksacijski vektor omega.

Na koncu sem SOR iteracijo uporabil še za vlaganje grafov v ravnino / prostor s fizikalno metodo. Kjer sem poskušal najti tudi najoptimalnejši parameter omega.

## Uporaba kode
1. V terminalu zaženemo Julio
2. Pritisnemo "]" in izvedemo naslednje ukaze:
    1. ```activate Domaca01```
    2. ```add LinearAlgebra```
    3. ```add Test```
    4. ```add Plots```
3. Nato lahko v Julia terminalu poganjamo posamezne vrstice kode (v Visual Studio Code s stiskom na Shift+Enter na vrstici ki jo želimo izvesti).

## Poganjanje testov
1. V terminalu zaženemo Julio
2. Pritisnemo "]" in izvedemo naslednje ukaze:
    1. ```activate Domaca01```
    2. ```test```
3. Testi se izvedejo

## Generiranje poročila
1. V terminalu zaženemo Julio
2. Pritisnemo "]" in izvedemo: ```activate Domaca01```
3. Pritisnemo "Backspace", da gremo iz pkg
4. Izvedemo naslednje ukaze:
    1. ```cd("Domaca01\\scripts")```
    2. ```include("makedocs.jl")```
5. Ustvari se poročilo z imenom "demo.pdf"