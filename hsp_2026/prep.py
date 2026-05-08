'''
prep.py — Módulo de Pré-Processamento de Imagens de Micro-FTIR
==============================================================
Contém todas as funções de pré-processamento espectral necessárias antes
da análise multivariada (K-Means, PCA, EMSC, etc.).

Todas as funções recebem e retornam o dicionário hsp padrão, modificando-o
in-place e registrando cada operação no campo 'log'. Isso permite rastrear
exatamente quais etapas foram aplicadas em cada conjunto de dados.

Pipeline típico de pré-processamento:
    1. cut()    — restringe a faixa espectral de interesse
    2. golay()  — suavização e/ou derivadas (Savitzky-Golay)
    3. norm()   — normalização vetorial (SNV)
    4. offset() — remoção de offset em região plana

Funções públicas:
    cut(data, a, b)                    — recorte espectral entre a e b cm⁻¹
    norm(data)  / snv(data)            — normalização vetorial (SNV)
    golay(data, diff, order, win)      — filtro / derivada Savitzky-Golay
    norm2r(data, ini1, fim1, ini2, fim2)— SNV em duas regiões concatenadas
    pcares(data, n)                    — remoção de ruído por PCA
    napc(dados, noise, npcs)           — remoção de ruído por NAPC
    offset(data, ini, fim)             — remoção de offset pelo mínimo
    binned(data)                       — binagem 2×2 pixels
    rand(data, k)                      — seleção aleatória de k pixels
    dsample(data)                      — subamostragem espacial
'''


def cut(data, a, b):
    '''
    Recorta o espectro mantendo apenas os números de onda entre a e b cm⁻¹.

    Útil para focar a análise em regiões espectrais de interesse biológico
    (ex.: região de impressão digital 900–1800 cm⁻¹) e remover regiões
    com alta absorção de água ou com baixo sinal-ruído.

    Parâmetros
    ----------
    data : dict hsp
        Deve conter 'r' (espectros) e 'wn' (números de onda).
    a : float
        Limite inferior do intervalo espectral em cm⁻¹ (inclusive).
    b : float
        Limite superior do intervalo espectral em cm⁻¹ (inclusive).

    Retorna
    -------
    data : dict hsp
        'r'  — reduzido para n_pontos_selecionados colunas
        'wn' — vetor reduzido ao intervalo [a, b]
    '''
    import numpy as np

    # Cria seletores booleanos para cada limite
    sel1 = (data['wn'] > a)   # True para wn acima do limite inferior
    sel2 = (data['wn'] < b)   # True para wn abaixo do limite superior
    # Interseção: True apenas onde ambas as condições são satisfeitas
    # A operação aritmética equivale a um AND lógico de forma vetorizada
    ver = (sel1.astype(int) + sel2.astype(int)) - 1
    sel = ver.astype(bool)

    data['r']  = data['r'][:, sel]    # seleciona apenas as colunas desejadas
    data['wn'] = data['wn'][sel]

    linha = '\n restrição espectral de ' + str(a) + ' cm-1 até ' + str(b) + ' cm-1'
    print(linha, end='')
    data['log'] = np.char.add(data['log'], linha)
    return data


def norm(data):
    '''
    Aplica normalização vetorial SNV (Standard Normal Variate) a todos os espectros.

    A SNV é a normalização mais comum em espectroscopia vibracional.
    Para cada espectro individualmente, subtrai a média e divide pelo
    desvio padrão:
        z_norm = (z - μ) / σ

    Isso remove variações multiplicativas (diferenças de espessura,
    espalhamento físico) e variações de offset (diferenças de linha base),
    deixando apenas a forma espectral relativa.

    Parâmetros
    ----------
    data : dict hsp
        Deve conter 'r' (espectros).

    Retorna
    -------
    data : dict hsp
        'r' — espectros normalizados (média=0, desvio=1 por espectro)
    '''
    import numpy as np

    spc   = data['r']
    media = np.mean(spc, axis=1)   # média de cada espectro (shape: n_espectros,)
    std   = np.std(spc, axis=1)    # desvio padrão de cada espectro

    # [:,None] expande para (n_espectros, 1) para broadcasting correto com (n_espectros, n_pontos)
    data['r'] = np.divide((spc - media[:, None]), std[:, None])

    linha = '\n normalização vetorial(SNV) em única região'
    print(linha, end='')
    data['log'] = np.char.add(data['log'], linha)
    return data


def snv(data):
    '''
    Alias para norm(). Aplica SNV (Standard Normal Variate).

    Mantido para compatibilidade com código que usa a nomenclatura 'snv'
    diretamente.

    Parâmetros
    ----------
    data : dict hsp

    Retorna
    -------
    data : dict hsp (ver norm())
    '''
    return norm(data)


def golay(data, diff, order, win):
    '''
    Aplica o filtro de Savitzky-Golay para suavização e/ou cálculo de derivadas.

    O filtro Savitzky-Golay ajusta um polinômio local de grau 'order' a uma
    janela deslizante de 'win' pontos e calcula a derivada de ordem 'diff'.
    É a técnica padrão em espectroscopia para:
        - diff=0: suavização sem distorção das bandas espectrais
        - diff=1: primeira derivada (remove offset linear, realça bandas)
        - diff=2: segunda derivada (remove baseline curva, resolve sobreposições)

    Implementação matricial eficiente:
        1. Calcula os coeficientes SG usando savgol_coeffs
        2. Monta uma matriz esparsa D de convolução com diagonais deslocadas
        3. Aplica via produto matricial: r_filtered = r · D

    Parâmetros
    ----------
    data : dict hsp
        Deve conter 'r' (espectros).
    diff : int
        Ordem da derivada (0 = apenas suavização, 1 = 1ª derivada, etc.)
    order : int
        Grau do polinômio de ajuste local. Valores típicos: 2 ou 3.
    win : int
        Tamanho da janela (número de pontos). Deve ser ímpar e > order.
        Janelas maiores → mais suavização, menos resolução espectral.

    Retorna
    -------
    data : dict hsp
        'r' — espectros filtrados/derivados
    '''
    import numpy as np
    from scipy.signal import savgol_coeffs
    from scipy.sparse import spdiags
    import numpy.matlib

    n = int((win - 1) / 2)   # metade da janela (raio)

    # Calcula os coeficientes do filtro SG para a janela e derivada especificadas
    sgcoeff = savgol_coeffs(win, order, deriv=diff)[:, None]
    # Replica os coeficientes para cada ponto espectral (n_pontos cópias)
    sgcoeff = np.matlib.repmat(sgcoeff, 1, data['r'].shape[1])

    # Monta a matriz de convolução usando diagonais:
    # cada diagonal i contém o coeficiente para o deslocamento correspondente
    diags = np.arange(-n, n + 1)
    D = spdiags(sgcoeff, diags, data['r'].shape[1], data['r'].shape[1]).toarray()

    # Zera as bordas da matriz para evitar artefatos de borda
    D[:, 0:n] = 0
    D[:, data['r'].shape[1] - 5:data['r'].shape[1]] = 0

    # Aplica o filtro via multiplicação matricial (equivalente a convolução)
    data['r'] = np.dot(data['r'], D)

    linha  = '\n filtro savitz golay usando \n'
    linha += '  >> derivada ordem: ' + str(diff) + '\n'
    linha += '  >> janela: ' + str(win) + '\n'
    linha += '  >> polinômio: ' + str(order) + ' ordem'
    print(linha, end='')
    data['log'] = np.char.add(data['log'], linha)
    return data


def norm2r(data, ini1, fim1, ini2, fim2):
    '''
    Aplica SNV separadamente em duas regiões espectrais e as concatena.

    Útil quando se quer normalizar e comparar duas faixas espectrais distintas
    que têm intensidades muito diferentes (ex.: região de proteínas 1500–1700 cm⁻¹
    + região de lipídios 2800–3000 cm⁻¹), evitando que uma domine a normalização.

    Cada região é normalizada independentemente pela sua própria média e desvio
    padrão, depois os espectros são concatenados.

    Parâmetros
    ----------
    data : dict hsp
    ini1, fim1 : float — limites em cm⁻¹ da primeira região
    ini2, fim2 : float — limites em cm⁻¹ da segunda região

    Retorna
    -------
    data : dict hsp
        'r'  — espectros com as duas regiões normalizadas e concatenadas
        'wn' — vetor de números de onda concatenado (região1 + região2)
    '''
    import numpy as np

    # --- Normaliza a primeira região ---
    sel  = np.logical_and(data['wn'] > int(ini1), data['wn'] < int(fim1))
    r1   = data['r'][:, sel]
    wn1  = data['wn'][sel][:, None]
    media = np.mean(r1, axis=1)
    std   = np.std(r1, axis=1)
    r1 = np.divide((r1 - media[:, None]), std[:, None])

    # --- Normaliza a segunda região ---
    sel  = np.logical_and(data['wn'] > int(ini2), data['wn'] < int(fim2))
    r2   = data['r'][:, sel]
    wn2  = data['wn'][sel][:, None]
    media = np.mean(r2, axis=1)
    std   = np.std(r2, axis=1)
    r2 = np.divide((r2 - media[:, None]), std[:, None])

    # Concatena as duas regiões normalizadas
    data['r']  = np.column_stack((r1, r2))
    data['wn'] = np.vstack((wn1, wn2)).reshape(-1)

    linha  = '\n normalização vetorial(SNV) em 2 regiões\n'
    linha += '  >> r1: ' + str(ini1) + ' cm-1 até ' + str(fim1) + ' cm-1\n'
    linha += '  >> r2: ' + str(ini2) + ' cm-1 até ' + str(fim2) + ' cm-1'
    print(linha, end='')
    data['log'] = np.char.add(data['log'], linha)
    return data


def pcares(data, n):
    '''
    Remove ruído espectral por reconstrução PCA com n componentes principais.

    A ideia é que os primeiros n PCs capturam a variância "real" (sinal
    biológico), enquanto os PCs seguintes capturam ruído aleatório do detector.
    Ao zerar os scores dos PCs além do n-ésimo e reconstruir os espectros,
    obtém-se uma versão filtrada com ruído reduzido.

    Algoritmo:
        1. Centraliza os espectros pela média global
        2. Projeta no espaço de PCA (calcula scores)
        3. Zera scores das componentes > n
        4. Reconstrói: r_denoised = média + scores_filtrados · loadings

    Parâmetros
    ----------
    data : dict hsp
        Deve conter 'r' (espectros).
    n : int
        Número de componentes principais a reter.
        Valores maiores → mais sinal preservado, menos ruído removido.

    Retorna
    -------
    data : dict hsp
        'r' — espectros reconstruídos com n PCs
    '''
    import numpy as np
    from sklearn.decomposition import PCA

    print('inicializando o pcares')
    pca   = PCA()
    media = np.mean(data['r'], axis=0)   # espectro médio global

    # Ajusta o modelo PCA nos espectros centrados
    pca.fit(data['r'] - media)
    scoress = pca.transform(data['r'])   # projeta no espaço de PCA

    # Zera scores das componentes além da n-ésima (filtro de baixa dimensionalidade)
    scoress[:, n - 1:-1] = 0

    coeff     = pca.components_   # loadings (vetores espectrais do PCA)
    data['r'] = media + np.dot(scoress, coeff)   # reconstrói espectros filtrados

    linha = '\n remoção de ruído usando somente redução de PCA com ' + str(n) + ' pcs'
    print(linha, end='')
    data['log'] = np.char.add(data['log'], linha)
    return data


def napc(dados, noise, npcs):
    '''
    Remove ruído espectral pelo método NAPC (Noise-Adjusted Principal Components).

    O NAPC é uma extensão do PCA que leva em conta explicitamente a estrutura
    de covariância do ruído instrumental (medido separadamente), produzindo
    componentes que maximizam a razão sinal-ruído em vez da variância total.

    Algoritmo (baseado em Green et al., 1988):
        1. Calcula as matrizes de covariância dos dados (Σ_d) e do ruído (Σ_n)
        2. Diagonaliza Σ_n via SVD: Σ_n = E₁ · diag(s₁) · E₁ᵀ
        3. Calcula a matriz de branqueamento do ruído: F = E₁ / √s₁
        4. Diagonaliza a covariância dos dados branqueada: F^T · Σ_d · F = G · diag(b) · G^T
        5. A transformação NAPC é: H = F · G
        6. Projeta os dados: scores = H^T · (dados - média)
        7. Zera os scores além de npcs e reconstrói

    Parâmetros
    ----------
    dados : dict hsp
        Dados a serem filtrados. Deve conter 'r'.
    noise : dict hsp
        Espectros de ruído puro (ex.: espectros de substrato sem amostra,
        ou diferença entre espectros de uma mesma região adquiridos duas vezes).
    npcs : int
        Número de componentes NAPC a reter (correspondentes às direções
        espectrais com maior razão sinal-ruído).

    Retorna
    -------
    dados : dict hsp
        'r' — espectros filtrados pelo NAPC
    '''
    import numpy as np

    d = dados['r']
    n = noise['r']

    print('inicializando NAPC')

    # Matrizes de covariância dos dados e do ruído
    sigmad = np.cov(d.T)
    sigman = np.cov(n.T)

    # SVD da covariância do ruído para obter sua decomposição espectral
    [a, s1, e1] = np.linalg.svd(sigman)
    e1 = e1.T   # autovetores como colunas

    # Matriz de branqueamento do ruído: F transforma Σ_n → identidade
    # Cada coluna de e1 é dividida pela raiz quadrada do autovalor correspondente
    F = e1 / np.sqrt(s1)

    # Covariância dos dados no espaço branqueado pelo ruído
    sigma_adj = F.T @ sigmad @ F

    # SVD da covariância ajustada para obter as direções NAPC
    [a, b, G] = np.linalg.svd(sigma_adj)
    G = G.T

    # Transformação NAPC completa: H mapeia do espaço original para o espaço NAPC
    H = F.dot(G)

    # Centraliza os dados pela média
    meanspc  = np.tile(d.mean(0).reshape(1, -1), (d.shape[0], 1))
    meandata = d - meanspc

    # Projeta e filtra: zera as componentes além de npcs
    scoresNAPC          = H.T @ meandata.T
    scoresNAPC[npcs:, :] = 0

    # Reconstrói usando a transformação inversa: H^T · z = scores → z = (H^T)⁻¹ · scores
    zcorr = np.linalg.solve(H.T, scoresNAPC).T

    dados['r'] = zcorr + meanspc   # reinsere a média

    linha = '\n remoção de ruído usando NAPC com ' + str(npcs) + ' pcs'
    print(linha, end='')
    dados['log'] = np.char.add(dados['log'], linha)
    return dados


def offset(data, ini, fim):
    '''
    Remove offset vertical dos espectros usando o valor mínimo em uma região plana.

    Para cada espectro, encontra o valor mínimo na região [ini, fim] cm⁻¹
    (onde não deve haver bandas de absorção) e subtrai esse valor de todo
    o espectro. Isso garante que a linha de base na região de referência
    fique próxima de zero.

    Útil para corrigir deslocamentos verticais causados por diferenças de
    espessura, emissão de fundo do substrato ou deriva do instrumento.

    Parâmetros
    ----------
    data : dict hsp
    ini : float — limite inferior da região de referência em cm⁻¹
    fim : float — limite superior da região de referência em cm⁻¹

    Retorna
    -------
    data : dict hsp
        'r' — espectros com offset removido (mínimo ≈ 0 na região [ini, fim])
    '''
    import numpy as np

    # Seleciona a região de referência para cálculo do offset
    sel    = np.logical_and(data['wn'] > int(ini), data['wn'] < int(fim))
    r      = data['r'][:, sel]

    # O offset de cada espectro é o seu valor mínimo nessa região
    minimo = np.min(r, axis=1)
    minimo = np.reshape(minimo, (-1, 1))
    # Replica o offset para subtrair de todos os pontos do espectro
    minimo = np.tile(minimo, data['r'].shape[1])

    data['r'] = data['r'] - minimo

    linha = '\n remoção de offset usando o mínimo valor entre ' + str(ini) + ' e ' + str(fim)
    print(linha, end='')
    data['log'] = np.char.add(data['log'], linha)
    return data


def binned(data):
    '''
    Aplica binagem 2×2: agrupa blocos de 4 pixels em um único pixel médio.

    A binagem reduz a resolução espacial pela metade em cada dimensão,
    mas aumenta significativamente a razão sinal-ruído (SNR ≈ 2× em amplitude)
    porque a média de 4 pixels independentes reduz o ruído aleatório por √4=2.

    Útil quando:
        - O mapa original tem resolução espacial maior que a difração permite
        - É necessário reduzir o número de pontos para análises mais rápidas
        - O SNR individual dos pixels é muito baixo

    Algoritmo:
        - Redimensiona r de (n_espectros, n_z) para (dx, dy, n_z)
        - Para cada bloco 2×2, calcula a média dos 4 espectros
        - Resultado: imagem de (⌊dx/2⌋-1, ⌊dy/2⌋-1) pixels

    Parâmetros
    ----------
    data : dict hsp
        Deve conter 'r', 'dx', 'dy'.

    Retorna
    -------
    data : dict hsp
        'r'   — espectros binados: shape ((dxbin*dybin), n_z)
        'sel' — nova máscara (todos True)
        'dx'  — nova dimensão vertical (≈ dx/2)
        'dy'  — nova dimensão horizontal (≈ dy/2)
    '''
    import numpy as np

    r  = data['r']
    # Redimensiona para grade 3D (linhas, colunas, pontos espectrais)
    r  = r.reshape(data['dx'], data['dy'], -1)
    dx = r.shape[0]
    dy = r.shape[1]
    dz = r.shape[2]

    # Calcula dimensões da imagem binada (-1 para evitar problema de borda)
    dxbin = int(np.floor(dx / 2)) - 1
    dybin = int(np.floor(dy / 2)) - 1
    rbin  = np.ones((dxbin, dybin, dz))

    jj = 0
    ii = 0
    # Percorre a imagem em passos de 2 pixels (blocos 2×2)
    for i in range(0, dy - 2, 2):
        for j in range(0, dx - 2, 2):
            sel = r[j:j+2, i:i+2, :]          # seleciona bloco 2×2×dz
            sel = np.mean(sel.reshape(4, dz), axis=0)  # média dos 4 pixels
            rbin[jj, ii, :] = sel
            jj += 1
        jj = 0
        ii += 1

    data['r']   = rbin.reshape((dxbin * dybin, dz))
    data['sel'] = np.ones((dxbin * dybin,)).astype('bool')
    data['dx']  = dxbin
    data['dy']  = dybin

    linha = '\n dados binados 2x2'
    print(linha, end='')
    data['log'] = np.char.add(data['log'], linha)
    return data


def rand(data, k):
    '''
    Seleciona aleatoriamente k pixels da imagem para análise.

    Útil para testar parâmetros de processamento em um subconjunto pequeno
    antes de aplicar em toda a imagem, ou para reduzir o custo computacional
    de algoritmos lentos (ex.: K-Means com muitos clusters em imagens grandes).

    A seleção é feita com reposição (randint), portanto pode haver repetições
    se k > n_espectros. Para seleção sem reposição, considere np.random.choice.

    Parâmetros
    ----------
    data : dict hsp
    k : int
        Número de pixels a selecionar aleatoriamente.

    Retorna
    -------
    data : dict hsp
        'r'   — subconjunto de k espectros selecionados
        'sel' — máscara atualizada (True apenas nos pixels selecionados)
    '''
    import numpy as np

    # Gera índices aleatórios no intervalo [0, n_espectros)
    points = np.random.randint(data['r'].shape[0], size=k)
    sel    = np.zeros_like(data['sel'])
    sel[points] = True          # marca os pixels selecionados

    data['r']           = data['r'][sel][:]
    data['sel'][data['sel']] = sel   # atualiza a máscara global

    linha = '\n seleção de ' + str(k) + ' pixels de maneira aleatória'
    print(linha, end='')
    data['log'] = np.char.add(data['log'], linha)
    return data


def dsample(data):
    '''
    Realiza subamostragem espacial (downsampling) com passo de 2 pixels.

    Remove pixels alternados da grade espacial, reduzindo o número de espectros
    aproximadamente pela metade. Diferentemente de binned(), não faz média —
    simplesmente descarta pixels, o que preserva a resolução espectral mas
    reduz a densidade espacial.

    A seleção é feita zerando linhas e colunas alternadas na grade 2D
    e usando a máscara resultante para filtrar os espectros.

    Parâmetros
    ----------
    data : dict hsp
        Deve conter 'r', 'dx', 'dy', 'sel'.

    Retorna
    -------
    data : dict hsp
        'r'   — espectros subasmostrados
        'sel' — máscara atualizada
    '''
    import numpy as np

    n   = 2
    sel = np.ones((data['dx'], data['dy']))

    # Índices das linhas e colunas a manter (passo n=2)
    XX = list(range(0, sel.shape[0] - 1, n))
    YY = list(range(0, sel.shape[1] - 2, n))

    # Zera as linhas e colunas dos índices selecionados (mantém as demais)
    sel[XX, :] = 0
    sel[:, YY] = 0

    sel = sel.reshape(-1,).astype('bool')

    data['r']           = data['r'][sel, :]
    data['sel'][data['sel']] = sel
    return data