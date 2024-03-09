"""
Trabalho 2 - Analise Nodal Modificada e Resposta AC

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

    def indutor(self, noA, noB, iX, L, w):
        self.Gn[noA][iX] += 1
        self.Gn[noB][iX] -= 1
        self.Gn[iX][noA] -= 1
        self.Gn[iX][noB] += 1
        self.Gn[iX][iX] += 1j * w * L

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

def main(arqNetlist, tipo, nos, parametros):

    if tipo == 'AC':
        freqInicial = parametros[0]
        freqFinal = parametros[1]
        ptsPorDec = parametros[2]

        # freqs = np.logspace(int(np.log10(freqInicial)), int(np.log10(freqFinal)),
        #                     num=int(ptsPorDec*(np.log10(freqFinal)-np.log10(freqInicial))))

        freqs = np.logspace(int(np.log10(freqInicial)), int(np.log10(freqFinal)), num=ptsPorDec)

        omegas = 2*np.pi*freqs
    else:
        freqs = [0]
        omegas = [0]

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

    # contagemNos e uma lista que vai servir para contar a
    # quantidade de nos do circuito, possibilitando determinar
    # os tamanhos das matrizes Gn e In.
    contagemNos = list()

    # variaveisCorrente funciona como um contador de quantas variaveis
    # de corrente terao de ser adicionadas na analise nodal modificada.
    variaveisCorrente = 0

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

    numNos = max(contagemNos)
    tamanhoMatriz = numNos + variaveisCorrente

    q = 0

    modulos = np.zeros((len(nos), len(freqs)), dtype=np.complex128)
    fases = np.zeros((len(nos), len(freqs)), dtype=np.complex128)

    # A analise nodal modificada sera feita para cada valor diferente de frequencia.
    for w in omegas:
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
            if tipoFonte == 'AC':
                if tipo == 'AC':
                    amplitude = float(i.split(' ')[4])
                    fase = float(i.split(' ')[-1])*(np.pi/180)
                    corrente = amplitude*np.exp(1j*fase)
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
            indutancia = float(l.split(' ')[-1])
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
            analiseNodalMod.transformador(noA, noB, noC, noD, iX, iY, L1, L2, M, w)

        for c in capacitores:
            noA = int(c.split(' ')[1])
            noB = int(c.split(' ')[2])
            capacitancia = float(c.split(' ')[-1])
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
            if tipoFonte == 'AC':
                if tipo == 'AC':
                    amplitude = float(v.split(' ')[4])
                    fase = float(v.split(' ')[-1])*(np.pi/180)
                    tensao = amplitude*np.exp(1j*fase)
                else:
                    tensao = 0
            else:
                if tipo == 'DC':
                    tensao = float(v.split(' ')[-1])
                else:
                    tensao = 0
            analiseNodalMod.fonteTensao(noA, noB, iX, tensao)

        # Calculo das tensoes nodais
        e = analiseNodalMod.tensoesNodais(tamanhoMatriz)

        j = 0

        # Explicitando as tensoes nodais em modulo e fase.
        for n in nos:
            if tipo == 'AC':
                modulos[j][q] = 20*np.log10(np.abs(e[n - 1]))
                fases[j][q] = np.degrees(np.angle(e[n - 1]))
            else:
                modulos[j][q] = e[n-1]

            j += 1

        q += 1

    if tipo == 'AC':
        # fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(10, 10), sharex=False)
        #
        # ax[0].set_title('Magnitude e Fase ' + arqNetlist)
        #
        # linhas, colunas = modulos.shape
        #
        # for i in range(0, linhas):
        #     # print(np.real(modulos[i]))
        #     ax[0].semilogx(freqs, np.real(modulos[i]))
        #     ax[1].semilogx(freqs, np.real(fases[i]))
        #
        # ax[0].set_xlabel('Freq [Hz]')
        # ax[0].set_ylabel('Magnitude [dB]')
        # ax[0].grid(True, which='both', ls='-')
        #
        # ax[1].set_xlabel('Freq [Hz]')
        # ax[1].set_ylabel('Fase (graus)')
        # ax[1].grid(True, which='both', ls='-')
        #
        # plt.show()

        return freqs, np.real(modulos), np.real(fases)
    else:
        return np.real(modulos)


if __name__ == '__main__':

    # print('Netlist DC 1:')
    # print(main('./Trabalho2NodalModificadaAC/netlistDC1.txt', 'DC', [2], []))
    # # Resultado esperado: [6.]
    #
    # print('')
    #
    # print('Netlist DC 2:')
    # print(main('./Trabalho2NodalModificadaAC/netlistDC2.txt', 'DC', [2, 3, 5, 7, 9, 10], []))
    # # Resultado esperado: [8, 1, 8.4, 0.93333333, -5.6, -3.73333333]
    #
    # print('')
    #
    # print('Netlist DC 3:')
    # print(main('./Trabalho2NodalModificadaAC/netlistDC3.txt', 'DC', [1, 2, 3, 4, 5, 6, 7], []))
    # # Resultado esperado: [10., 5.7554841, 2.49323451, 4.37608413, 2.28100871, 5.7554841, 6.37608413]
    #
    # print('')
    #
    # print('Netlist DC 4:')
    # print(main('./Trabalho2NodalModificadaAC/netlistDC4.txt', 'DC', [2], []))
    # # Resultado esperado: [6.]
    #
    # print('')
    #
    # print('Netlist DC 5:')
    # print(main('./Trabalho2NodalModificadaAC/netlistDC5.txt', 'DC', [2], []))
    # # Resultado esperado: [10.]
    #
    # print('')
    #
    # print('Netlist DC 6:')
    # print(main('./Trabalho2NodalModificadaAC/netlistDC6.txt', 'DC', [3, 4, 5], []))
    # # Resultado esperado: [0.5, 0, 0]
    #
    # print('')

    #print('Netlist AC 1')
    freqs, modulos, fases = main('./Trabalho2NodalModificadaAC/netlistAC8.txt', 'AC', [4], [100, 100e3, 100])
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(10, 10), sharex=False)

    ax[0].set_title('Magnitude e Fase')

    linhas, colunas = modulos.shape

    for i in range(0, linhas):
        # print(np.real(modulos[i]))
        ax[0].semilogx(freqs, np.real(modulos[i]))
        ax[1].semilogx(freqs, np.real(fases[i]))

    ax[0].set_xlabel('Freq [Hz]')
    ax[0].set_ylabel('Magnitude [dB]')
    ax[0].grid(True, which='both', ls='-')

    ax[1].set_xlabel('Freq [Hz]')
    ax[1].set_ylabel('Fase (graus)')
    ax[1].grid(True, which='both', ls='-')

    plt.show()
    #
    # print('Netlist AC 2')
    # main('./Trabalho2NodalModificadaAC/netlistAC2.txt','AC',[1], [0.01, 200, 100])
    #
    # print('Netlist AC 3')
    # main('./Trabalho2NodalModificadaAC/netlistAC3.txt', 'AC', [2], [0.01, 100, 100])
    #
    # print('Netlist AC 4')
    # main('./Trabalho2NodalModificadaAC/netlistAC4.txt', 'AC', [2, 3], [0.01, 1000, 1000])
    #
    # print('Netlist AC 5')
    # main('./Trabalho2NodalModificadaAC/netlistAC5.txt', 'AC', [3], [0.01, 1000, 1000])
    #
    # print('Netlist AC 6')
    # main('./Trabalho2NodalModificadaAC/netlistAC6.txt', 'AC', [2, 5], [0.01, 1e4, 1000])
    #
    # print('Netlist AC 7')
    # main('./Trabalho2NodalModificadaAC/netlistAC7.txt', 'AC', [2, 7], [0.01, 100, 1000])
    #
    # print('Netlist AC 8')
    # main('./Trabalho2NodalModificadaAC/netlistAC8.txt', 'AC', [4], [100, 100e3, 100])
    #
    # print('Netlist AC 9')
    # main('./Trabalho2NodalModificadaAC/netlistAC9.txt', 'AC', [2, 3, 4, 5, 6], [0.01, 100, 1000])
    #
    # print('Netlist AC 10')
    # main('./Trabalho2NodalModificadaAC/netlistAC10.txt', 'AC', [4, 5], [0.01, 1000, 1000])
