'''
rg.py — Módulo de Histogramas para Definição de Parâmetros de Qualidade
========================================================================
Fornece funções de visualização de histogramas para auxiliar na escolha
dos limites a e b usados nos filtros de qualidade do módulo qt.py.

O workflow típico é:
    1. rg.area(data, ini, fim)   → analisa a distribuição de áreas
    2. Observa o histograma e identifica o corte visual entre sinal e ruído
    3. Usa os limites identificados em qt.area(data, ini, fim, a, b)

O mesmo ciclo se aplica para intt() e emsc().

Todas as funções geram histogramas com 300 bins e grade ligada por padrão.

Funções públicas:
    area(data, ini, fim)   — histograma das áreas integradas de uma banda
    intt(data, b)          — histograma das intensidades num número de onda
    emsc(data, a)          — histograma de um coeficiente do modelo EMSC
'''


def area(data, ini, fim):
    '''
    Plota o histograma das áreas integradas de todos os pixels na faixa
    espectral [ini, fim] cm⁻¹.

    A área integrada (regra do trapézio) é calculada para cada espectro na
    região selecionada. O histograma permite visualizar a distribuição das
    áreas e identificar populações de pixels (ex.: sinal vs. fundo).

    Use este gráfico antes de chamar qt.area() para escolher os limites
    a (mínimo) e b (máximo) de forma informada.

    Parâmetros
    ----------
    data : dict hsp
        Deve conter 'r' (espectros) e 'wn' (números de onda).
    ini : float — limite inferior da faixa de integração em cm⁻¹
    fim : float — limite superior da faixa de integração em cm⁻¹
    '''
    import numpy as np
    import matplotlib.pyplot as plt

    # Seleciona os pontos na faixa de interesse e calcula área por trapézio
    sel  = np.logical_and(data['wn'] > int(ini), data['wn'] < int(fim))
    r    = data['r'][:, sel]
    area = np.trapz(r)   # shape: (n_espectros,)

    plt.figure()
    plt.hist(area, 300)   # 300 bins para boa resolução da distribuição

    linha = 'histograma de área entre ' + str(ini) + ' até ' + str(fim)
    plt.title(linha)
    plt.grid()


def intt(data, b):
    '''
    Plota o histograma das intensidades espectrais de todos os pixels
    no número de onda b (cm⁻¹).

    A intensidade num ponto espectral específico é uma medida pontual da
    absorção naquele número de onda. Use este histograma para visualizar a
    distribuição de intensidades antes de chamar qt.intt().

    Parâmetros
    ----------
    data : dict hsp
        Deve conter 'r' (espectros) e 'wn' (números de onda).
    b : float — número de onda de interesse em cm⁻¹
    '''
    import numpy as np
    import matplotlib.pyplot as plt

    # Seleciona todos os pontos acima de b e pega a primeira coluna
    # (ponto mais próximo do valor b no vetor wn)
    sel = data['wn'] > b
    ver = data['r'][:, sel]
    ver = ver[:, 0]   # intensidade no ponto mais próximo de b

    plt.figure()
    plt.hist(ver, 300)

    linha = 'histograma da intensidade em ' + str(b) + ' cm-1'
    plt.title(linha)
    plt.grid()


def emsc(data, a):
    '''
    Plota o histograma do coeficiente de índice 'a' do modelo EMSC
    para todos os pixels.

    Os coeficientes EMSC têm significado físico:
        - Índice 0 (a₀): escala da referência → proporcional à espessura/concentração
        - Índices 1..n: termos polinomiais de baseline
        - Demais: coeficientes de parafina e água

    O histograma do a₀ é especialmente útil para separar pixels de amostra
    (valores positivos altos) de pixels de fundo (valores próximos de zero).

    Parâmetros
    ----------
    data : dict hsp
        Deve conter 'EMSC_coeff' (gerado por emsc_fit() do módulo emsc.py).
    a : int — índice do coeficiente EMSC a visualizar (0-based)
    '''
    import numpy as np
    import matplotlib.pyplot as plt

    ver = data['EMSC_coeff'][:, a]   # coeficiente 'a' de todos os espectros

    plt.figure()
    plt.hist(ver, 300)

    linha = 'histograma do coeficiente ' + str(a) + ' do modelo de EMSC'
    plt.title(linha)