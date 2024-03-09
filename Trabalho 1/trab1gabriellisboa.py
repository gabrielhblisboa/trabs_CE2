'''
Trabalho 1

Aluno: Gabriel Henrique Braga Lisboa
Disciplina: Circuitos Elétricos II
Professora: Fernanda Duarte Vilela Reis de Oliveira
Turma: EL1
Período: 2023/2
'''

import numpy as np


class AnaliseNodal():

    # Inicializa a Matriz de Condutancias e a matriz de correntes
    # assim que o objeto da classe e criado.
    def __init__(self, numNos):
        self.Gn = np.zeros([numNos + 1, numNos + 1])
        self.In = np.zeros(numNos + 1)

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
    def fonteCorrenteControlada(self, noA, noB, noC, noD, Gm):
        self.Gn[noA][noC] += Gm
        self.Gn[noA][noD] -= Gm
        self.Gn[noB][noC] -= Gm
        self.Gn[noB][noD] += Gm

    # Calcula as tensoes nodais do circuito.
    def tensoesNodais(self, numNos):
        self.e = np.zeros(numNos - 1)
        e = np.linalg.inv(self.Gn[1:, 1:]) @ self.In [1:]
        # print('')
        # print('Matriz de Tensoes Nodais:')
        # print(e)
        # print('')
        # i=1
        # for tensaoNodal in e:
        #     print('e{} = {} V'.format(i, tensaoNodal))
        #     i += 1
        # print('')
        return e
         

def main(arqNetlist):

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
    resistores, fontes, fontesControladas = list(), list(), list()

    # contagemNos e uma lista que vai servir para contar a
    # quantidade de nos do circuito, possibilitando determinar
    # os tamanhos das matrizes Gn e In.
    contagemNos = list()

    for linha in netlist:
        linha = linha.strip()
        if linha.startswith('*') or linha == '':
            # Desconsidera os comentários ou as linhas vazias da netlist.
            continue
        else:
            # Detecta se a linha da netlist equivale a um resistor, fonte de
            # corrente independente ou fonte de corrente controlada por tensao e
            # conta os nos desse componente (e da variavel de controle, se for
            # uma fonte controlada). Depois, insere a linha na respectiva lista
            # de componentes.
            if linha[0] == 'G':
                contagemNos.append(int(linha.split(' ')[1]))
                contagemNos.append(int(linha.split(' ')[2]))
                contagemNos.append(int(linha.split(' ')[3]))
                contagemNos.append(int(linha.split(' ')[4]))
                fontesControladas.append(linha)
            elif linha[0] == 'I':
                contagemNos.append(int(linha.split(' ')[1]))
                contagemNos.append(int(linha.split(' ')[2]))
                fontes.append(linha)
            elif linha[0] == 'R':
                contagemNos.append(int(linha.split(' ')[1]))
                contagemNos.append(int(linha.split(' ')[2]))
                resistores.append(linha)

    numNos = max(contagemNos)

    analiseNodal = AnaliseNodal(numNos)

    # Adicao da estampa dos resistores na matriz de condutancias.
    for r in resistores:
        noA = int(r.split(' ')[1])
        noB = int(r.split(' ')[2])
        resistencia = int(r.split(' ')[-1])
        analiseNodal.resistor(noA, noB, resistencia)

    # Adicao da estampa das fontes de corrente
    # independentes na matriz de correntes.
    for i in fontes:
        noA = int(i.split(' ')[1])
        noB = int(i.split(' ')[2])
        corrente = int(i.split(' ')[-1])
        analiseNodal.fonteCorrente(noA, noB, corrente)

    # Adicao da estampa das fontes de corrente
    # controladas por tensao na matriz de condutancias.
    for g in fontesControladas:
        noA = int(g.split(' ')[1])
        noB = int(g.split(' ')[2])
        noC = int(g.split(' ')[3])
        noD = int(g.split(' ')[4])
        transcondutancia = int(g.split(' ')[-1])
        analiseNodal.fonteCorrenteControlada(noA, noB, noC, noD, transcondutancia)

    return analiseNodal.tensoesNodais(numNos)


if __name__ == '__main__':
    print('')
    print('Teste de Erro:')
    print(main('teste_de_erro.txt'))

    print('')
    print('--------------------------------------------------------------------------------')

    print('')
    print('Netlist 1:')
    print(main('netlist1.txt'))

    print('--------------------------------------------------------------------------------')
    print('')

    print('Netlist 2:')
    print(main('netlist2.txt'))

    print('--------------------------------------------------------------------------------')
    print('')

    print('Netlist 3:')
    print(main('netlist3.txt'))
