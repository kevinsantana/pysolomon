from collections import OrderedDict

class Galois():
    def __init__(self, primitivo: int, tamanho: int) -> list:
        self.primitivo = primitivo
        self.tamanho = tamanho
        self.corpo_galois = OrderedDict()        
        self.inicializar_galois()


    def inicializar_galois(self):
        # chave: alfa; valores: forma decimal        
        x = 1        
        self.corpo_galois["zero"] = 0
        for indice in range(self.tamanho - 1):            
            self.corpo_galois[indice] = x
            x = x * 2
            if x >= self.tamanho:
               x ^= self.primitivo        

                
    def soma_subtracao(self, a: int, b: int) -> int:
        '''
        Implementa a soma e subtração entre dois elementos de um corpo de Galois, que estejam no mesmo GF(tamanho)        
        '''
        return a ^ b


    def multiplicacao(self, a: int, b: int) -> int:
        # a ou b sao decimais em um gf        
        if a == 0 or b == 0:
            return 0 
        return self.corpo_galois[(self.grau(a) + self.grau(b)) % (self.tamanho - 1)]

    
    def elemento_inverso(self, a: int) -> int:
        if a == 0:
            raise ArithmeticError()        
        return self.corpo_galois[(- self.grau(a) % (self.tamanho - 1))]

    
    def divisao(self, a: int, b: int) -> int:
        # a ou b sao decimais em um gf # a / b      
        return self.multiplicacao(a, self.elemento_inverso(b))


    def grau(self, elemento: int):
        return [chave for chave, valor in self.corpo_galois.items() if elemento == valor].pop()


if __name__ == "__main__":
    gf_16 = Galois(19,16)        
    print(gf_16.divisao(11,10))
    print(gf_16.corpo_galois)
