from galois import *

class polinomio_galois():
    # coeficientes: do maior para o menor grau
    def __init__(self, corpo: galois, coeficientes: list):
        if len(coeficientes) == 0:
            raise ArithmeticError()
        self.corpo = corpo
        self.coeficientes = coeficientes

    
    def multiplicacao_polinomios(self, polinomio) -> list:
        return [galois.multiplicacao(self.corpo,a,b) for a in self.coeficientes for b in polinomio]


gf_16 = galois(19,16)        
coeficientes = [1,2,3,4,5,6,7,8,9,10,11]
polinomio = polinomio_galois(gf_16, coeficientes)
coeficientes2 = [3]
multiplicaco = polinomio.multiplicacao_polinomios(coeficientes2)
print(gf_16.decimais.index(3))
