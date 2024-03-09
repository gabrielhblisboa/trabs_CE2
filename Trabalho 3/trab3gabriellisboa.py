"""
Trabalho 3 - Análise no Tempo e Análise de Circuitos Não-Lineares

Aluno: Gabriel Henrique Braga Lisboa
Disciplina: Circuitos Elétricos II
Professora: Fernanda Duarte Vilela Reis de Oliveira
Turma: EL1
Período: 2023/2
"""

import numpy as np
import matplotlib.pyplot as plt


class AnaliseNodal():

    # Inicializa a Matriz de Condutancias e a matriz de correntes
    # assim que o objeto da classe e criado.
    def __init__(self, numNos):
        self.Gn = np.zeros([numNos + 1, numNos + 1], dtype=np.complex128)
        self.In = np.zeros((numNos + 1), dtype=np.complex128)

    # Estampa do resistor.
    def resistor(self, noA, noB, R):
        self.Gn[noA, noA] += 1/R
        self.Gn[noA, noB] -= 1/R
        self.Gn[noB, noA] -= 1/R
        self.Gn[noB, noB] += 1/R

    # Estampa da fonte de corrente independente.
    def fonteCorrente(self, noA, noB, I):
        self.In[noA] -= I
        self.In[noB] += I

    # Estampa da fonte de corrente controlada por tensao.
    def fonteCorrenteControladaPorTensao(self, noA, noB, noC, noD, Gm):
        self.Gn[noA][noC] += Gm
        self.Gn[noA][noD] -= Gm
        self.Gn[noB][noC] -= Gm
        self.Gn[noB][noD] += Gm

    # Estampa do capacitor.
    def capacitor(self, noA, noB, C, w):
        self.Gn[noA][noA] += 1j*w*C
        self.Gn[noA][noB] -= 1j*w*C
        self.Gn[noB][noA] -= 1j*w*C
        self.Gn[noB][noB] += 1j*w*C

    # Calcula as tensoes nodais do circuito.
    def tensoesNodais(self, numNos):
        self.e = np.zeros((numNos), dtype=np.complex128)
        e = np.linalg.solve(self.Gn[1:, 1:], self.In[1:])
        return e

class AnaliseModificada(AnaliseNodal):

    def __init__(self, tamanhoMatriz):
        super().__init__(tamanhoMatriz)

    def capTrapMod(self, noA, noB, iX, C, passo, vt0, it0):
        self.Gn[noA][iX] += 1
        self.Gn[noB][iX] -= 1
        self.Gn[iX][noA] -= 1
        self.Gn[iX][noB] += 1
        self.Gn[iX][iX] += passo/(2*C)

        self.In[iX] = -(vt0+(passo/(2*C))*it0)

    def indutor(self, noA, noB, iX, L, w):
        self.Gn[noA][iX] += 1
        self.Gn[noB][iX] -= 1
        self.Gn[iX][noA] -= 1
        self.Gn[iX][noB] += 1
        self.Gn[iX][iX] += 1j * w * L

    def indTrapMod(self, noA, noB, iX, L, passo, vt0, it0):
        self.Gn[noA][iX] += 1
        self.Gn[noB][iX] -= 1
        self.Gn[iX][noA] -= 1
        self.Gn[iX][noB] += 1
        self.Gn[iX][iX] += (2*L)/passo

        self.In[iX] = ((2*L)/passo)*it0 + vt0

    def transformador(self, noA, noB, noC, noD, iX, iY, L1, L2, M, w):
        self.Gn[noA][iX] += 1
        self.Gn[noB][iX] -= 1
        self.Gn[noC][iY] += 1
        self.Gn[noD][iY] -= 1
        self.Gn[iX][noA] -= 1
        self.Gn[iX][noB] += 1
        self.Gn[iY][noC] -= 1
        self.Gn[iY][noD] += 1
        self.Gn[iX][iX] += 1j * w * L1
        self.Gn[iX][iY] += 1j * w * M
        self.Gn[iY][iX] += 1j * w * M
        self.Gn[iY][iY] += 1j * w * L2

    def transfTrapMod(self, noA, noB, noC, noD, iX, iY, L1, L2, M, passo, vabt0, vcdt0, iabt0, icdt0):
        self.Gn[noA][iX] += 1
        self.Gn[noB][iX] -= 1
        self.Gn[noC][iY] += 1
        self.Gn[noD][iY] -= 1
        self.Gn[iX][noA] -= 1
        self.Gn[iX][noB] += 1
        self.Gn[iY][noC] -= 1
        self.Gn[iY][noD] += 1
        self.Gn[iX][iX] += (2*L1)/passo
        self.Gn[iX][iY] += (2*M)/passo
        self.Gn[iY][iX] += (2*M)/passo
        self.Gn[iY][iY] += (2*L2)/passo

        self.In[iX] += (2/passo)*(L1*iabt0 + M*icdt0) + vabt0
        self.In[iY] += (2/passo)*(M*iabt0 + L2*icdt0) + vcdt0


    def fonteCorrenteControladaPorCorrente(self, noA, noB, noC, noD, iX, B):
        self.Gn[noA][iX] += B
        self.Gn[noB][iX] -= B
        self.Gn[noC][iX] += 1
        self.Gn[noD][iX] -= 1
        self.Gn[iX][noC] -= 1
        self.Gn[iX][noD] += 1

    def fonteTensaoControladaPorTensao(self, noA, noB, noC, noD, iX, A):
        self.Gn[noA][iX] += 1
        self.Gn[noB][iX] -= 1
        self.Gn[iX][noA] -= 1
        self.Gn[iX][noB] += 1
        self.Gn[iX][noC] += A
        self.Gn[iX][noD] -= A

    def fonteTensaoControladaPorCorrente(self, noA, noB, noC, noD, iX, iY, Rm):
        self.Gn[noA][iY] += 1
        self.Gn[noB][iY] -= 1
        self.Gn[noC][iX] += 1
        self.Gn[noD][iX] -= 1
        self.Gn[iY][noA] -= 1
        self.Gn[iY][noB] += 1
        self.Gn[iX][noC] -= 1
        self.Gn[iX][noD] += 1
        self.Gn[iY][iX]  += Rm

    def fonteTensao(self, noA, noB, iX, V):
        self.Gn[noA][iX] += 1
        self.Gn[noB][iX] -= 1
        self.Gn[iX][noA] -= 1
        self.Gn[iX][noB] += 1

        self.In[iX] -= V


def calcularCircuito(componentes, numNos, tamanhoMatriz, tipo, valIniciais, w, t, passo, primeiroCalculo):

    (resistores, fontesCorrente, fontesCorrenteCtrlTensao, indutores,
     transformadores, capacitores, fontesTensao, fontesTensaoCtrlTensao,
     fontesTensaoCtrlCorrente, fontesCorrenteCtrlCorrente, diodos) = (
        componentes["resistores"], componentes["fontesCorrente"],
        componentes["fontesCorrenteCtrlTensao"], componentes["indutores"],
        componentes["transformadores"], componentes["capacitores"],
        componentes["fontesTensao"], componentes["fontesTensaoCtrlTensao"],
        componentes["fontesTensaoCtrlCorrente"], componentes["fontesCorrenteCtrlCorrente"],
        componentes["diodos"]
    )

    analiseNodalMod = AnaliseModificada(tamanhoMatriz)

    # contadorCorrentes funciona como indice das correntes da matriz Gn
    contadorCorrentes = 1

    # Adicao das estampas dos componentes nas matrizes de condutancias
    # e de correntes.
    for r in resistores:
        noA = int(r.split(' ')[1])
        noB = int(r.split(' ')[2])
        resistencia = float(r.split(' ')[-1])
        analiseNodalMod.resistor(noA, noB, resistencia)

    for i in fontesCorrente:
        noA = int(i.split(' ')[1])
        noB = int(i.split(' ')[2])
        tipoFonte = i.split(' ')[3]
        if tipoFonte == 'SIN':
            if tipo == 'TRAN':
                valorDC = float(i.split(' ')[4])
                amplitude = float(i.split(' ')[5])
                freq = float(i.split(' ')[6])
                fase = float(i.split(' ')[-1]) * (np.pi / 180)
                corrente = (amplitude * np.cos(2 * np.pi * freq * t + fase)) + valorDC
            else:
                corrente = 0
        else:
            if tipo == 'DC':
                corrente = float(i.split(' ')[-1])
            else:
                corrente = 0
        analiseNodalMod.fonteCorrente(noA, noB, corrente)

    for g in fontesCorrenteCtrlTensao:
        noA = int(g.split(' ')[1])
        noB = int(g.split(' ')[2])
        noC = int(g.split(' ')[3])
        noD = int(g.split(' ')[4])
        transcondutancia = float(g.split(' ')[-1])
        analiseNodalMod.fonteCorrenteControladaPorTensao(noA, noB, noC, noD, transcondutancia)

    for l in indutores:
        noA = int(l.split(' ')[1])
        noB = int(l.split(' ')[2])
        iX = numNos + contadorCorrentes
        contadorCorrentes += 1
        indutancia = float(l.split(' ')[3])
        if tipo == 'TRAN':
            vt0 = valIniciais[noA] - valIniciais[noB]
            it0 = valIniciais[iX]
            analiseNodalMod.indTrapMod(noA, noB, iX, indutancia, passo, vt0, it0)
        else:
            analiseNodalMod.indutor(noA, noB, iX, indutancia, w)

    for k in transformadores:
        noA = int(k.split(' ')[1])
        noB = int(k.split(' ')[2])
        noC = int(k.split(' ')[3])
        noD = int(k.split(' ')[4])
        iX = numNos + contadorCorrentes
        iY = numNos + contadorCorrentes + 1
        contadorCorrentes += 2
        L1 = float(k.split(' ')[5])
        L2 = float(k.split(' ')[6])
        M = float(k.split(' ')[-1])
        if tipo == 'TRAN':
            vabt0 = valIniciais[noA] - valIniciais[noB]
            vcdt0 = valIniciais[noC] - valIniciais[noD]
            iabt0 = valIniciais[iX]
            icdt0 = valIniciais[iY]
            analiseNodalMod.transfTrapMod(noA, noB, noC, noD, iX, iY, L1, L2, M, passo, vabt0, vcdt0, iabt0, icdt0)
        else:
            analiseNodalMod.transformador(noA, noB, noC, noD, iX, iY, L1, L2, M, w)

    for c in capacitores:
        noA = int(c.split(' ')[1])
        noB = int(c.split(' ')[2])
        capacitancia = float(c.split(' ')[3])
        if tipo == 'TRAN':
            iX = numNos + contadorCorrentes
            contadorCorrentes += 1
            vt0 = valIniciais[noA] - valIniciais[noB]
            if primeiroCalculo:
                it0 = 0
            else:
                it0 = valIniciais[iX]
            analiseNodalMod.capTrapMod(noA, noB, iX, capacitancia, passo, vt0, it0)
        else:
            analiseNodalMod.capacitor(noA, noB, capacitancia, w)

    for f in fontesCorrenteCtrlCorrente:
        noA = int(f.split(' ')[1])
        noB = int(f.split(' ')[2])
        noC = int(f.split(' ')[3])
        noD = int(f.split(' ')[4])
        iX = numNos + contadorCorrentes
        contadorCorrentes += 1
        ganhoCorrente = float(f.split(' ')[-1])
        analiseNodalMod.fonteCorrenteControladaPorCorrente(noA, noB, noC, noD, iX, ganhoCorrente)

    for e in fontesTensaoCtrlTensao:
        noA = int(e.split(' ')[1])
        noB = int(e.split(' ')[2])
        noC = int(e.split(' ')[3])
        noD = int(e.split(' ')[4])
        iX = numNos + contadorCorrentes
        contadorCorrentes += 1
        ganhoTensao = float(e.split(' ')[-1])
        analiseNodalMod.fonteTensaoControladaPorTensao(noA, noB, noC, noD, iX, ganhoTensao)

    for h in fontesTensaoCtrlCorrente:
        noA = int(h.split(' ')[1])
        noB = int(h.split(' ')[2])
        noC = int(h.split(' ')[3])
        noD = int(h.split(' ')[4])
        iX = numNos + contadorCorrentes
        iY = numNos + contadorCorrentes + 1
        contadorCorrentes += 2
        transresistencia = float(h.split(' ')[-1])
        analiseNodalMod.fonteTensaoControladaPorCorrente(noA, noB, noC, noD, iX, iY, transresistencia)

    for v in fontesTensao:
        noA = int(v.split(' ')[1])
        noB = int(v.split(' ')[2])
        iX = numNos + contadorCorrentes
        contadorCorrentes += 1
        tipoFonte = v.split(' ')[3]
        if tipoFonte == 'SIN':
            if tipo == 'TRAN':
                valorDC = float(v.split(' ')[4])
                amplitude = float(v.split(' ')[5])
                freq = float(v.split(' ')[6])
                fase = float(v.split(' ')[-1]) * (np.pi / 180)
                tensao = (amplitude * np.cos(2*np.pi*freq*t + fase)) + valorDC
            else:
                tensao = 0
        else:
            if tipo == 'DC':
                tensao = float(v.split(' ')[-1])
            else:
                tensao = 0
        analiseNodalMod.fonteTensao(noA, noB, iX, tensao)

    for d in diodos:
        noPositivo = int(d.split(' ')[1])
        noNegativo = int(d.split(' ')[2])
        Is = float(d.split(' ')[3])
        nvT = float(d.split(' ')[-1])
        vn = valIniciais[noPositivo] - valIniciais[noNegativo]
        if vn > 1:
            vn = 1
        elif vn < -2:
            vn = -2
        Go = (Is * np.exp(vn / nvT)) / nvT
        Io = Is * (np.exp(vn / nvT) - 1) - Go * vn
        analiseNodalMod.resistor(noPositivo, noNegativo, 1/Go)
        analiseNodalMod.fonteCorrente(noPositivo, noNegativo, Io)

    # Calculo das tensoes nodais
    vetorMod = analiseNodalMod.tensoesNodais(tamanhoMatriz)

    return vetorMod


def main(arqNetlist, tipo, nos, parametros):

    # Bloco para leitura do arquivo da netlist. 
    # Caso o arquivo não exista, é retornado um código de erro.
    try:
        arquivo = open(arqNetlist, 'r')
    except FileNotFoundError:
        print('Erro: Arquivo não encontrado.')
        return -1
    else:
        netlist = arquivo.readlines()
        arquivo.close()

    # Listas para separar as linhas da netlist de acordo
    # com o tipo de componente associado aquela linha.
    resistores, fontesCorrente, fontesCorrenteCtrlTensao = list(), list(), list()
    indutores, transformadores, capacitores = list(), list(), list()
    fontesTensao, fontesTensaoCtrlTensao, fontesTensaoCtrlCorrente = list(), list(), list()
    fontesCorrenteCtrlCorrente = list()
    diodos = list()

    # contagemNos e uma lista que vai servir para contar a
    # quantidade de nos do circuito, possibilitando determinar
    # os tamanhos das matrizes Gn e In.
    contagemNos = list()

    # variaveisCorrente funciona como um contador de quantas variaveis
    # de corrente terao de ser adicionadas na analise nodal modificada.
    variaveisCorrente = 0

    numCapacitores = 0

    compNaoLin = False

    for linha in netlist:
        linha = linha.strip()
        if linha.startswith('*') or linha == '':
            # Desconsidera os comentários ou as linhas vazias da netlist.
            continue
        else:
            # Determina qual o componente que a respectiva linha da netlist
            # esta descrevendo e conta os nos desse componente (e da variavel de controle, se for
            # uma fonte controlada). Depois, insere a linha na respectiva lista
            # de componentes.
            if linha[0] == 'G':
                contagemNos.append(int(linha.split(' ')[1]))
                contagemNos.append(int(linha.split(' ')[2]))
                contagemNos.append(int(linha.split(' ')[3]))
                contagemNos.append(int(linha.split(' ')[4]))
                fontesCorrenteCtrlTensao.append(linha)
            elif linha[0] == 'I':
                contagemNos.append(int(linha.split(' ')[1]))
                contagemNos.append(int(linha.split(' ')[2]))
                fontesCorrente.append(linha)
            elif linha[0] == 'R':
                contagemNos.append(int(linha.split(' ')[1]))
                contagemNos.append(int(linha.split(' ')[2]))
                resistores.append(linha)
            elif linha[0] == 'L':
                contagemNos.append(int(linha.split(' ')[1]))
                contagemNos.append(int(linha.split(' ')[2]))
                variaveisCorrente += 1
                indutores.append(linha)
            elif linha[0] == 'K':
                contagemNos.append(int(linha.split(' ')[1]))
                contagemNos.append(int(linha.split(' ')[2]))
                contagemNos.append(int(linha.split(' ')[3]))
                contagemNos.append(int(linha.split(' ')[4]))
                variaveisCorrente += 2
                transformadores.append(linha)
            elif linha[0] == 'C':
                contagemNos.append(int(linha.split(' ')[1]))
                contagemNos.append(int(linha.split(' ')[2]))
                numCapacitores += 1
                capacitores.append(linha)
            elif linha[0] == 'F':
                contagemNos.append(int(linha.split(' ')[1]))
                contagemNos.append(int(linha.split(' ')[2]))
                contagemNos.append(int(linha.split(' ')[3]))
                contagemNos.append(int(linha.split(' ')[4]))
                variaveisCorrente += 1
                fontesCorrenteCtrlCorrente.append(linha)
            elif linha[0] == 'E':
                contagemNos.append(int(linha.split(' ')[1]))
                contagemNos.append(int(linha.split(' ')[2]))
                contagemNos.append(int(linha.split(' ')[3]))
                contagemNos.append(int(linha.split(' ')[4]))
                variaveisCorrente += 1
                fontesTensaoCtrlTensao.append(linha)
            elif linha[0] == 'H':
                contagemNos.append(int(linha.split(' ')[1]))
                contagemNos.append(int(linha.split(' ')[2]))
                contagemNos.append(int(linha.split(' ')[3]))
                contagemNos.append(int(linha.split(' ')[4]))
                variaveisCorrente += 2
                fontesTensaoCtrlCorrente.append(linha)
            elif linha[0] == 'V':
                contagemNos.append(int(linha.split(' ')[1]))
                contagemNos.append(int(linha.split(' ')[2]))
                variaveisCorrente += 1
                fontesTensao.append(linha)
            elif linha[0] == 'D':
                compNaoLin = True
                contagemNos.append(int(linha.split(' ')[1]))
                contagemNos.append(int(linha.split(' ')[2]))
                diodos.append(linha)

    componentes = {
        "resistores": resistores,
        "fontesCorrente": fontesCorrente,
        "fontesCorrenteCtrlTensao": fontesCorrenteCtrlTensao,
        "indutores": indutores,
        "transformadores": transformadores,
        "capacitores": capacitores,
        "fontesTensao": fontesTensao,
        "fontesTensaoCtrlTensao": fontesTensaoCtrlTensao,
        "fontesTensaoCtrlCorrente": fontesTensaoCtrlCorrente,
        "fontesCorrenteCtrlCorrente": fontesCorrenteCtrlCorrente,
        "diodos": diodos
    }

    numNos = max(contagemNos)
    tamanhoMatriz = numNos + variaveisCorrente

    primeiroCalculo = True

    if tipo == 'DC':

        tol = parametros[0]
        valIniciais = np.array(parametros[-1])
        w = 0
        t = 0
        passo = 0

        tensoesNodais = np.zeros(len(nos), dtype=np.complex128)

        for x in range(0, 100):

            vetorMod = calcularCircuito(componentes, numNos,
                                        tamanhoMatriz, tipo, valIniciais, w, t, passo, primeiroCalculo)

            if variaveisCorrente != 0:
                e = vetorMod[:numNos]
            else:
                e = vetorMod

            if not compNaoLin:
                break

            if np.max(np.abs(e - valIniciais[1:])) < tol:
                break

            valIniciais = np.hstack((np.array([0]), e))

        indice = 0
        for n in nos:
            tensoesNodais[indice] = e[n - 1]
            indice += 1

        return np.real(tensoesNodais)

    elif tipo == 'TRAN':
        tSimulacao = parametros[0]
        passo = parametros[1]
        tol = parametros[2]
        valIniciais = np.array(parametros[-1])
        w = 0

        tempo = np.linspace(0, tSimulacao, int(tSimulacao/passo))

        tensoesNodais = np.zeros((len(nos), len(tempo)), dtype=np.complex128)

        for x in range(0, 100):
            vetorMod = calcularCircuito(componentes, numNos,
                                        tamanhoMatriz, 'DC', valIniciais, w, 0, passo, primeiroCalculo)

            if variaveisCorrente != 0:
                e = vetorMod[:numNos]
            else:
                e = vetorMod

            auxiliar = np.hstack((np.array([0]), vetorMod))

            valIniciaisTeste = valIniciais[:numNos + 1]

            valIniciais = auxiliar

            if not compNaoLin:
                break

            if np.max(np.abs(e - valIniciaisTeste[1:])) < tol:
                break

        indice = 0
        for n in nos:
            tensoesNodais[indice][0] = e[n - 1]
            indice += 1

        tamanhoMatriz += numCapacitores

        primeiroCalculo = False

        for q, t in enumerate(tempo):
            if t == 0:
                continue

            vetorMod = (
                calcularCircuito(componentes, numNos, tamanhoMatriz, tipo, valIniciais, w, t, passo, primeiroCalculo))

            if variaveisCorrente != 0:
                e = vetorMod[:numNos]
            else:
                e = vetorMod

            valIniciais = np.hstack((np.array([0]), vetorMod))

            indice = 0
            for n in nos:
                tensoesNodais[indice][q] = e[n - 1]
                indice += 1

        return tempo, np.real(tensoesNodais)


if __name__ == '__main__':

    print('Testes DC:')
    print('')

    print('Teste 1:')
    print(main('testesDC/teste1.txt', 'DC', [1, 2], [1e-10, [0, 0.1, 0.1]]))
    print('')
    # Resultado esperado: [2.0, 1.30491003]

    print('Teste 2:')
    print(main('testesDC/teste2.txt', 'DC', [2], [1e-14, [0,3,3]]))
    # Resultado esperado: [0.73068013]
    print('')

    print('Teste 3:')
    print(main('testesDC/teste3.txt', 'DC', [2], [1e-12, [0,0,2]]))
    # Resultado esperado: [2]
    print('')

    print('Teste 4:')
    print(main('testesDC/teste4.txt', 'DC', [2], [1e-12, [0,0,2]]))
    # Resultado esperado: [1.9]
    print('')

    print('Teste 5:')
    print(main('testesDC/teste5.txt', 'DC', [1,2], [1e-12, [0,1,2]]))
    # Resultado esperado: [5, 4.2781936]
    print('')

    print('Teste 6:')
    print(main('testesDC/teste6.txt', 'DC', [1,2,3], [1e-14, [0,5,1,-5]]))
    # Resultado esperado: [ 5,  4.2781936, -5]
    print('')

    print('Teste 7:')
    print(main('testesDC/teste7.txt', 'DC', [1,2,3], [1e-14, [0,-5,4,5]]))
    # Resultado esperado: [ -5,  4.2781936,  5]
    print('')

    tempo, tensoesNodais = main('./testesTRAN/testeTran1.txt', 'TRAN',
                                [1, 2], [2, 0.2e-3, 1e-10, [0, 1, 0.5]])

    # tempo, tensoesNodais = main('./testesTRAN/testeTran2.txt', 'TRAN',
    #                             [1, 2], [2, 0.2e-3, 1e-4, [0, 1, 0]])

    # tempo, tensoesNodais = main('./testesTRAN/testeTran3.txt', 'TRAN',
    #                             [1, 2], [3, 0.2e-3, 1e-4, [0, 0, 0]])
    #
    # tempo, tensoesNodais = main('./testesTRAN/testeTran4.txt', 'TRAN',
    #                             [1, 2], [3, 0.2e-3, 1e-4, [0, 0, 0]])
    #
    # tempo, tensoesNodais = main('./testesTRAN/testeTran5.txt', 'TRAN',
    #                             [1, 2], [3, 0.2e-3, 1e-4, [0, 9.5, 0]])
    #
    # tempo, tensoesNodais = main('./testesTRAN/testeTran6.txt', 'TRAN',
    #                             [1, 3], [2, 0.5e-4, 1e-4, [0, 0, 0, 0]])
    #
    # tempo, tensoesNodais = main('./testesTRAN/testeTran7.txt', 'TRAN',
    #                             [1, 2, 3], [2, 0.5e-4, 1e-4, [0, 0, 0, 0]])
    #
    # tempo, tensoesNodais = main('./testesTRAN/testeTran8.txt', 'TRAN',
    #                             [1, 2, 3], [4, 0.2e-3, 1e-4, [0, 0, 0, 0]])
    #
    # tempo, tensoesNodais = main('./testesTRAN/testeTran9.txt', 'TRAN',
    #                             [1, 2], [4, 0.2e-3, 1e-4, [0, 1, 1]])
    #
    # tempo, tensoesNodais = main('./testesTRAN/testeTran10.txt', 'TRAN',
    #                             [1, 2, 3], [4, 0.2e-3, 1e-4, [0, 5, 5, 0]])

    plt.title('Análise no Tempo')

    linhas, colunas = tensoesNodais.shape

    for i in range(0, linhas):
        plt.plot(tempo, tensoesNodais[i])

    plt.xlabel('Tempo [s]')
    plt.ylabel('Tensão [V]')
    plt.grid(True, which='both', ls='-')

    plt.show()
