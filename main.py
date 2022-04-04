from Chromosome import Chromosome
import math
import random


def get_length(domeniu, precizie):
    subint = math.log2((domeniu[1] - domeniu[0]) * 10 ** precizie)
    return int(abs(round(subint)))


def f(x):
    return -x ** 2 + x + 2


def selectie(u, lst):
    l = 0
    r = len(lst) - 1
    m = (l + r) // 2
    while l < r:
        m = (l + r) // 2
        if lst[m] == u:
            return m
        elif u < lst[m]:
            r = m - 1
        else:
            l = m + 1
    return l


input = open("./input_files/testInput.txt", 'r')
output = open("output.txt", 'w')
n = 0
for line in input:
    line = line.split(" ")
    if line[0] == 'Populatie':
        n = int(line[1])
    elif line[0] == 'Precizie':
        precizie = int(line[1])
    elif line[0] == 'Recombinare':
        recomb = float(line[1])
    elif line[0] == 'Mutatie':
        mutatie = float(line[1])
    elif line[0] == 'Etape':
        etape = int(line[1])
    elif line[0] == 'Interval':
        domeniu = (int(line[1]), int(line[2]))

length = get_length(domeniu, precizie)
output.write("Populatia initiala \n")
maxList = []
populatie = []
suma_fx = 0
for _ in range(n):
    c = Chromosome(length, precizie, domeniu, [])
    c.set_random_genes()
    x = c.get_value()
    c.set_Fx(f(x))
    suma_fx += c.fx
    populatie.append(c)
    output.write(str(c))
    output.write('\n')

while (etape != 0):
    output.write("\n\nProbabilitati selectie \n")

    for cr in range(0, len(populatie)):
        p = populatie[cr].compute_selection_p(suma_fx)
        output.write(f"cromozom {cr + 1}\t probabilitate  {p}\n")

    s = 0
    output.write("\nIntervale probabilitati selectie \n")
    interv_selectie = [0]
    output.write('0 ')
    for cr in populatie:
        s += cr.get_p()
        interv_selectie.append(s)
        output.write(str(s) + ' ')
    output.write('\n')

    selectati = []
    for i in range(0, len(populatie)):
        u = random.random()
        nr = selectie(u, interv_selectie)
        print(nr)
        selectati.append(populatie[nr-1])
        output.write(f"u = {u}\t selectam cromozomul {nr}\n")

    output.write("\nDupa selectie: \n")
    for cr in range(len(selectati)):
        output.write(f"{cr + 1}: {str(selectati[cr])} \n")

    output.write(f"\nProbabilitatea de incrucisare {recomb}\n")
    recombinati = []
    index = []

    for i, cr in enumerate(selectati):
        u = random.random()
        if u >= recomb:
            output.write(f"{i}: {cr.genes_to_string()}\t u = {u}\n")
        else:
            recombinati.append(cr)
            index.append(i)
            output.write(f"{i+1}: {cr.genes_to_string()}\t u = {u} < {recomb} participa\n")

    if len(recombinati) % 2 != 0:
        recombinati.pop()

    for i in range(0, len(recombinati) - 1, 2):
        cr1 = recombinati[i]
        temp = cr1.genes_to_string()
        cr2 = recombinati[i + 1]
        u = random.randint(0, length + 1)
        output.write(f"Recombinare dintre cromozomul {index[i] + 1} cu cromozomul {index[i + 1] + 1}:\n")
        output.write(f"{cr1.genes_to_string()}\t {cr2.genes_to_string()}\t punctul {u}\n")
        output.write("Rezultat\n")
        recombinati[i].set_genes(cr1.incrucisare(u, cr2.genes_to_string()))
        recombinati[i + 1].set_genes(cr2.incrucisare(u, temp))
        output.write(f"{recombinati[i].genes_to_string()}\t{recombinati[i + 1].genes_to_string()}\n")

    output.write("\nDupa recombinare:\n")
    i = 0
    sumFx = 0
    for i, cr in enumerate(populatie):
        if i < len(index) and cr == index[i]:
            cr = recombinati[i]
            v = cr.get_value()
            cr.set_Fx(f(v))
            i += 1
        sumFx += cr.getFx()
        output.write(f"{i + 1}: {str(cr)} \n")

    for cr in populatie:
        cr.compute_selection_p(sumFx)

    output.write(f"\nProbabilitate de mutatie pentru fiecare gena {mutatie}\n")
    output.write("Au fost modificati cromozomii:\n")
    for i in index:
        output.write(str(i+1) + ' \n')


    output.write("\nDupa mutatie:\n")
    maxFx = populatie[0].getFx()
    sumFx = 0
    for i, cr in enumerate(populatie):
        u = random.random()
        if u < mutatie:
            pos = random.randint(0, length-1)
            cr.set_genes(cr.mutate(pos))
            cr.set_Fx(f(cr.get_value()))
        maxFx = max(maxFx, cr.getFx())
        sumFx += cr.getFx()
        output.write(f"{i + 1}: {str(cr)}\n")

    for cr in populatie:
        cr.compute_selection_p(sumFx)

    maxList.append(maxFx)

    etape -= 1

for i, cr in enumerate(maxList):
    output.write(str(cr) + '\n')
