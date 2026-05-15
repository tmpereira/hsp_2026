'''
qt.py — Módulo de Controle de Qualidade (Quality Test) para Micro-FTIR
=======================================================================
Contém filtros para remoção automática de pixels de baixa qualidade
espectral em imagens hiperespectrais de FTIR.

A qualidade dos pixels pode ser comprometida por:
    - Bolhas de ar ou falhas na preparação da amostra
    - Bordas da amostra com espessura insuficiente
    - Pixels de substrato (sem amostra)
    - Saturação do detector em regiões de alta absorção

Todos os filtros trabalham com critérios de inclusão/exclusão baseados em
limiares mínimo (a) e máximo (b) de alguma métrica espectral. Pixels cuja
métrica cai fora do intervalo [a, b] são removidos.

A máscara 'sel' do dicionário hsp é atualizada a cada operação, permitindo
rastrear quais pixels originais sobreviveram a todos os filtros aplicados.

Funções públicas:
    area(data, ini, fim, a, b)     — filtra por área integrada de uma banda
    intt(data, ini, a, b)          — filtra por intensidade num número de onda
    emsc(data, ini, a, b)          — filtra por coeficiente do modelo EMSC
    mean(data, ini, fim, a, b)     — filtra por escala do espectro médio
    otsu_area(data, ini1, fim1, k) — limiar automático Otsu sobre a área
    otsu_emsc(data)                — limiar automático Otsu sobre coef. EMSC a₀
'''


def area(data, ini, fim, a, b):
    '''
    Remove pixels cujo espectro tem área integrada fora do intervalo [a, b]
    na faixa espectral [ini, fim] cm⁻¹.

    A área integrada é calculada pela regra do trapézio (np.trapz) sobre os
    pontos espectrais na faixa selecionada. É uma medida robusta da quantidade
    de absorbância total na região e pode ser interpretada como proporcional
    à concentração e espessura local da amostra.

    Critério de remoção:
        - área < a → pixel muito fraco (substrato, bolha, borda fina)
        - área > b → pixel muito forte (saturação, precipitado, artefato)

    Parâmetros
    ----------
    data : dict hsp
    ini : float — limite inferior da faixa de integração em cm⁻¹
    fim : float — limite superior da faixa de integração em cm⁻¹
    a   : float — valor mínimo aceitável de área
    b   : float — valor máximo aceitável de área

    Retorna
    -------
    data : dict hsp
        'r'   — espectros filtrados
        'sel' — máscara atualizada (False nos pixels removidos)
    '''
    import numpy as np

    oldnspc = data['r'].shape[0]   # número de espectros antes do filtro

    # Seleciona os pontos na faixa de integração e calcula a área pelo trapézio
    sel  = np.logical_and(data['wn'] > int(ini), data['wn'] < int(fim))
    r    = data['r'][:, sel]
    area = np.trapz(r)             # integração numérica por regra do trapézio

    # Máscara: True para pixels dentro do intervalo aceitável
    sel = np.logical_and(area > float(a), area < float(b))

    data['r']            = data['r'][sel, :]
    data['sel'][data['sel']] = sel    # propaga a remoção para a máscara global

    newnspc = data['r'].shape[0]

    linha  = '\n teste de qualidade usando a área'
    linha += '\n região ' + str(ini) + ' até ' + str(fim)
    linha += '\n min_value: ' + str(a) + '\n max_value: ' + str(b)
    linha += '\n espectros removidos: ' + str(oldnspc - newnspc)
    print(linha, end='')
    data['log'] = np.char.add(data['log'], linha)
    return data


def intt(data, ini, a, b):
    '''
    Remove pixels cuja intensidade em um número de onda específico está
    fora do intervalo [a, b].

    Útil para filtrar pixels com base em uma banda característica específica,
    como a banda de amida I (≈1650 cm⁻¹) que indica presença de proteína,
    ou a banda de C=O (≈1740 cm⁻¹) que indica ésteres lipídicos.

    Parâmetros
    ----------
    data : dict hsp
    ini : float — número de onda de interesse em cm⁻¹ (ponto mais próximo)
    a   : float — intensidade mínima aceitável
    b   : float — intensidade máxima aceitável

    Retorna
    -------
    data : dict hsp
        'r'   — espectros filtrados
        'sel' — máscara atualizada
    '''
    import numpy as np

    oldnspc = data['r'].shape[0]

    # Encontra o índice do número de onda mais próximo de 'ini'
    sel  = data['wn'] == ini
    area = data['r'][:, sel]   # intensidade de cada espectro nesse ponto

    # Máscara de pixels dentro do intervalo aceitável
    sel = np.logical_and(area > float(a), area < float(b))
    sel = np.reshape(sel, (-1,))   # garante vetor 1D

    data['r']            = data['r'][sel, :]
    data['sel'][data['sel']] = sel

    newnspc = data['r'].shape[0]

    linha  = '\n teste de qualidade usando a intensidade'
    linha += '\n pico em ' + str(ini)
    linha += '\n min_value: ' + str(a) + '\n max_value: ' + str(b)
    linha += '\n espectros removidos: ' + str(oldnspc - newnspc)
    print(linha, end='')
    data['log'] = np.char.add(data['log'], linha)
    return data


def emsc(data, ini, a, b):
    '''
    Remove pixels cujo coeficiente EMSC de índice 'ini' está fora de [a, b].

    O modelo EMSC decompõe cada espectro em coeficientes interpretáveis:
        - Índice 0: coeficiente do espectro de referência (a₀) → relacionado
          à espessura/concentração da amostra
        - Índices 1..n: coeficientes de baseline polinomial
        - Demais: coeficientes da parafina e da água

    O coeficiente a₀ (índice 0) é o mais usado para controle de qualidade,
    pois é proporcional à quantidade total de material biológico no pixel.
    Pixels com a₀ muito baixo tipicamente correspondem a bordas ou bolhas.

    Parâmetros
    ----------
    data : dict hsp
        Deve conter 'EMSC_coeff' (gerado por emsc_fit() do módulo emsc.py).
    ini : int — índice do coeficiente EMSC a usar como critério (0-based)
    a   : float — valor mínimo aceitável
    b   : float — valor máximo aceitável

    Retorna
    -------
    data : dict hsp
        'r'   — espectros filtrados
        'sel' — máscara atualizada
    '''
    import numpy as np

    oldnspc = data['r'].shape[0]

    # Extrai o coeficiente de índice 'ini' para todos os espectros
    area = data['EMSC_coeff'][:, ini]

    sel = np.logical_and(area > float(a), area < float(b))
    sel = np.reshape(sel, (-1,))

    data['r']            = data['r'][sel, :]
    data['sel'][data['sel']] = sel

    newnspc = data['r'].shape[0]

    linha  = '\n teste de qualidade usando os coeficientes do EMSC'
    linha += '\n coef número ' + str(ini)
    linha += '\n min_value: ' + str(a) + '\n max_value: ' + str(b)
    linha += '\n espectros removidos: ' + str(oldnspc - newnspc)
    print(linha, end='')
    return data


def mean(data, ini, fim, a, b):
    '''
    Remove pixels cuja escala em relação ao espectro médio está fora de [a, b].

    Projeta cada espectro sobre o espectro médio da região via regressão linear
    mínimos quadrados (y = α · média + c), obtendo o coeficiente de escala α.
    Pixels com α muito diferente de 1 têm forma espectral ou amplitude muito
    diferente do padrão, podendo ser artefatos ou regiões heterogêneas.

    Algoritmo:
        Para cada espectro y na região:
            [α, c] = lstsq([média, 1], y)
        O coeficiente α mede o quanto y precisa ser escalado para se aproximar
        do espectro médio.

    Parâmetros
    ----------
    data : dict hsp
    ini : float — limite inferior da região de análise em cm⁻¹
    fim : float — limite superior da região de análise em cm⁻¹
    a   : float — coeficiente de escala mínimo aceitável
    b   : float — coeficiente de escala máximo aceitável

    Retorna
    -------
    data : dict hsp
        'r'   — espectros filtrados
        'sel' — máscara atualizada
    '''
    import numpy as np

    oldnspc = data['r'].shape[0]

    # Seleciona a região espectral de interesse
    sel = np.logical_and(data['wn'] > ini, data['wn'] < fim)
    r1  = data['r'][:, sel]

    # Monta o sistema de regressão: X = [média, vetor de uns]
    # alpha[0] = coeficiente de escala de cada espectro em relação à média
    y      = r1[:, :].T
    xx     = np.vstack((r1.mean(axis=0), np.ones_like(r1[1, :]))).T
    alpha  = np.linalg.lstsq(xx, y, rcond=-1)[0][0].T
    area   = alpha   # coeficiente de escala

    sel = np.logical_and(area > float(a), area < float(b))

    data['r']            = data['r'][sel, :]
    data['sel'][data['sel']] = sel

    newnspc = data['r'].shape[0]

    linha  = '\n teste de qualidade usando a meanspc'
    linha += '\n região ' + str(ini) + ' até ' + str(fim)
    linha += '\n min_value: ' + str(a) + '\n max_value: ' + str(b)
    linha += '\n espectros removidos: ' + str(oldnspc - newnspc)
    print(linha, end='')
    data['log'] = np.char.add(data['log'], linha)
    return data


def otsu_area(data, ini1, fim1, k=1):
    '''
    Remove pixels de baixa qualidade usando limiar automático de Otsu
    calculado sobre a área espectral na faixa [ini1, fim1] cm⁻¹.

    O método de Otsu determina automaticamente um limiar ótimo que maximiza
    a separação entre dois grupos (pixels de amostra vs. pixels de fundo/ruído),
    sem necessidade de definir manualmente os valores a e b.

    O parâmetro k permite ajustar o limiar: k=1 usa o limiar padrão do Otsu,
    k<1 torna o critério mais permissivo (retém mais pixels), k>1 mais restritivo.

    Parâmetros
    ----------
    data : dict hsp
    ini1 : float — limite inferior da faixa de integração em cm⁻¹
    fim1 : float — limite superior da faixa de integração em cm⁻¹
    k    : float, opcional (padrão=1) — fator multiplicativo do limiar Otsu

    Retorna
    -------
    data : dict hsp
        'r'   — espectros que passaram o limiar Otsu
        'sel' — máscara atualizada
    '''
    import numpy as np
    from skimage import filters

    oldnspc = data['r'].shape[0]

    # Calcula o coeficiente de escala pela mesma projeção de mean()
    sel  = np.logical_and(data['wn'] > ini1, data['wn'] < fim1)
    r1   = data['r'][:, sel]
    y    = r1[:, :].T
    xx   = np.vstack((r1.mean(axis=0), np.ones_like(r1[1, :]))).T
    alpha = np.linalg.lstsq(xx, y, rcond=-1)[0][0].T
    area  = alpha

    # Calcula o limiar de Otsu e aplica o fator k
    val = k * filters.threshold_otsu(area)
    print('\nlimiar de corte do Otsu foi ', k * val, '\n')

    # Mantém apenas pixels acima do limiar (pertencentes à amostra)
    sel = area > float(val)

    data['r']            = data['r'][sel, :]
    data['sel'][data['sel']] = sel

    newnspc = data['r'].shape[0]
    return data


def otsu_emsc(data):
    '''
    Remove pixels de baixa qualidade usando limiar automático de Otsu
    calculado sobre o coeficiente EMSC a₀ (índice 0).

    O coeficiente a₀ do modelo EMSC é proporcional à quantidade de material
    biológico no pixel. O método de Otsu separa automaticamente pixels de
    amostra (a₀ alto) de pixels de fundo/bolha/borda (a₀ baixo).

    Não requer parâmetros adicionais — o limiar é determinado completamente
    de forma automática pela distribuição dos coeficientes.

    Parâmetros
    ----------
    data : dict hsp
        Deve conter 'EMSC_coeff' (gerado por emsc_fit() do módulo emsc.py).

    Retorna
    -------
    data : dict hsp
        'r'   — espectros que passaram o limiar Otsu
        'sel' — máscara atualizada
    '''
    from skimage import filters

    # Extrai o coeficiente a₀ (índice 0) de todos os espectros
    d   = data['EMSC_coeff'][:, 0]
    ver = filters.threshold_otsu(d)   # limiar automático de Otsu

    print('limiar definido pelo Otsu para a0 foi ', ver)

    # Mantém apenas pixels com a₀ acima do limiar (amostra presente)
    sel = data['EMSC_coeff'][:, 0] > ver

    data['r']            = data['r'][sel, :]
    data['sel'][data['sel']] = sel
    return data
