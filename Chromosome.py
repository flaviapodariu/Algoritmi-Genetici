import math
import random


class Chromosome:
    def __init__(self, length, precision, domain, genes):
        self.length = length
        self.p = 0  # to be computed
        self.precision = precision
        self.domain = domain
        self.genes = genes
        self.value = 0
        self.fx = 0

    def __repr__(self):
        g = self.genes_to_string() + "\t"
        x = "x = " + str(self.value) + "  "
        f = "f = " + str(self.fx) + "  "
        return g + x + f

    def set_length(self):
        subintervals = math.log2((self.domain[1] - self.domain[0]) * 10 ** self.precision)
        self.length = int(abs(round(subintervals)))

    def getFx(self):
        return self.fx

    def get_length(self):
        return self.length

    def get_precision(self):
        return self.precision

    def get_p(self):
        return self.p

    def set_genes(self, genes):
        self.genes = [int(gene) for gene in genes]

    def set_Fx(self, fx):
        self.fx = fx

    def get_value(self):
        coef = (self.domain[1] - self.domain[0]) / (2 ** self.length - 1)
        base10 = ''
        trailing = True
        for gene in self.genes:
            if trailing and gene == 1:
                trailing = False

            if not trailing:
                base10 += str(gene)
        base10 = int(base10, 2)
        self.value = coef * base10 + self.domain[0]
        self.value = round(self.value, self.precision)
        return self.value

    def set_random_genes(self):
        self.genes = [random.randint(0, 1) for _ in range(self.length)]

    def genes_to_string(self):
        toStr = ''
        for gene in self.genes:
            toStr += str(gene)
        return toStr

    def compute_selection_p(self, sumFx):
        self.p = self.fx / sumFx
        return self.p

    def incrucisare(self, punct, schimb):
        gene = self.genes_to_string()
        gene = schimb[0:punct] + gene[punct:]
        return gene

    def mutate(self, pos):
        gene = self.genes_to_string()
        if gene[pos] == '0':
            gene = gene[:pos] + '1' + gene[pos:]
        else:
            gene = gene[:pos] + '0' + gene[pos:]

        return gene
