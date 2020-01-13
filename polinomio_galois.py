from galois import *


class PolinomioGalois(Galois):
    # coeficientes: do maior para o menor grau
    def __init__(self, corpo: Galois, coeficientes: dict):
        if len(coeficientes) == 0:
            raise ArithmeticError()
        self.corpo = corpo
        self.coeficientes = coeficientes

    
    def multiplicacao_polinomios_gf(self, polinomio) -> dict:
        # comparar e lançar exceção caso o corpo de um polinomio seja diferente do outro
        return {Galois.multiplicacao(self.corpo,a,b) for a in self.coeficientes.values() for b in polinomio.values()}


gf_16 = Galois(19,16)        
polinomio = {10: 1, 9: 2, 8: 3, 7: 4, 6: 5, 5: 6, 4: 7, 3: 8, 2: 9, 1: 10, 0: 11}
polinomio2 = {4: 1}
pgf1 = PolinomioGalois(gf_16, polinomio)
pgf1 = PolinomioGalois(gf_16, polinomio2)
multiplicaco = pgf1.multiplicacao_polinomios_gf(polinomio2)
print(multiplicaco)
