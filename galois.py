class galois():
    def __init__(self, primitivo: int, tamanho: int) -> list:
        self.primitivo = primitivo
        self.tamanho = tamanho
        self.decimais = list()        
        self.inicializar_galois()


    def inicializar_galois(self):
        # indices: alfa; valores: decimais        
        x = 1        
        for _ in range(self.tamanho):
            self.decimais.append(x)
            x = x * 2
            if x >= self.tamanho:
                x ^= self.primitivo
        self.inicializado = True    

                
    def soma_subtracao(self, a: int, b: int) -> int:
        # a ou b sao decimais
        return a ^ b


    def multiplicacao(self, a: int, b: int) -> int:
        # a ou b sao decimais
        # corrigir para suportar multiplicacao entre coeficientes de um polinomio
        if a == 0 or b == 0:
            return 0 
        return self.decimais[(self.decimais.index(a) + self.decimais.index(b)) % (self.tamanho - 1)]

    
    def elemento_inverso(self, a: int) -> int:
        if a == 0:
            raise ArithmeticError()        
        return self.decimais[(- self.decimais.index(a) % (self.tamanho - 1))]


    # def divisao(self, a: int, b: int) -> int:
    #     # a ou b sao decimais
    #     return self.decimais[(self.decimais.index(a) - self.decimais.index(b)) % (self.tamanho - 1)]

    
    def divisao(self, a: int, b: int) -> int:
        # a ou b sao decimais
        # a / b
        return self.multiplicacao(a, self.elemento_inverso(b))

