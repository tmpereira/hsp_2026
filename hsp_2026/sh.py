'''
sh.py — Módulo de Visualização de Imagens de Micro-FTIR
========================================================
Fornece funções para gerar imagens hiperespectrais a partir de diferentes
métricas espectrais, permitindo explorar visualmente a distribuição espacial
de componentes químicos na amostra.

Cada função reconstrói a grade 2D (dx × dy) a partir do vetor de espectros
selecionados (data['sel']), posicionando os valores calculados nos pixels
corretos e deixando os pixels removidos como zero.

Tipos de imagens disponíveis:
    - Intensidade num ponto espectral (intt, int_plt)
    - Área integrada de uma banda (area)
    - Escala relativa ao espectro médio (mean)
    - Coeficiente do modelo EMSC (emsc)
    - Score de componente principal (pc)
    - Espectros aleatórios (pplot)

Funções públicas:
    intt(data, b)              — imagem da intensidade no número de onda b
    area(data, a, b)           — imagem da área integrada entre a e b cm⁻¹
    mean(data, ini1, fim1)     — imagem da escala relativa ao espectro médio
    pplot(data, nspc)          — plota nspc espectros aleatórios
    emsc(data, b)              — imagem do coeficiente EMSC de índice b
    int_plt(data, b)           — imagem interativa com seleção de pixel por clique
    pc(data, n, k=10)          — imagem do score da PCn + loading espectral
'''


def intt(datta, b):
    '''
    Gera e exibe uma imagem baseada na intensidade espectral no número de onda b.

    A intensidade no ponto b é usada como o valor de cada pixel na imagem,
    resultando em um mapa de distribuição espacial da absorção naquela
    frequência. Útil para visualizar a distribuição de um componente específico:
        ex.: 1650 cm⁻¹ → proteínas (amida I)
             1740 cm⁻¹ → lipídios (C=O de ésteres)
             1085 cm⁻¹ → ácidos nucleicos (P=O)

    Parâmetros
    ----------
    datta : dict hsp
        Deve conter 'r', 'wn', 'dx', 'dy', 'sel', 'filename'.
    b : float
        Número de onda em cm⁻¹ para geração da imagem.
    '''
    import matplotlib.pyplot as plt
    import numpy as np

    # Seleciona todos os pontos com wn > b e pega o primeiro (mais próximo de b)
    sel = datta['wn'] > b
    ver = datta['r'][:, sel]
    ver = ver[:, 0]   # intensidade no ponto mais próximo de b

    # Reconstrói o mapa 2D: zeros no fundo, valores dos pixels válidos
    dplot = np.zeros(datta['dx'] * datta['dy'])
    dplot[datta['sel']] = ver
    dplot = dplot.reshape(datta['dx'], datta['dy'])

    plt.figure()
    # vmin/vmax definem a escala de cor com base nos dados (sem outliers por padrão)
    plt.pcolor(dplot, vmin=np.min(ver), vmax=np.max(ver))
    plt.clim(np.min(ver), np.max(ver))
    plt.colorbar()

    l = 'imagem da intensidade ' + str(b) + ' cm-1 \n ' + str(datta['filename'])[:-4]
    plt.title(l)
    plt.show()


def area(data, a, b):
    '''
    Gera e exibe uma imagem baseada na área integrada da banda entre a e b cm⁻¹.

    A área integrada (regra do trapézio) é proporcional à concentração do
    composto que absorve naquela região espectral, considerando a espessura
    da amostra. Produz imagens mais robustas que a intensidade pontual por
    ser menos sensível a deslocamentos de banda.

    Parâmetros
    ----------
    data : dict hsp
        Deve conter 'r', 'wn', 'dx', 'dy', 'sel', 'filename'.
    a : float — limite inferior da banda em cm⁻¹
    b : float — limite superior da banda em cm⁻¹
    '''
    import numpy as np
    import matplotlib.pyplot as plt

    # Seleciona pontos na banda de interesse e calcula área por trapézio
    sel  = (data['wn'] > a) & (data['wn'] < b)
    r    = data['r'][:, sel]
    area = np.trapz(r)   # área de cada pixel

    print(area.min)   # informação de debug: área mínima

    # Reconstrói o mapa 2D
    dplot = np.zeros(data['dx'] * data['dy'])
    dplot[data['sel']] = area
    dplot = dplot.reshape(data['dx'], data['dy'])

    plt.figure()
    plt.pcolor(dplot, vmin=np.min(area), vmax=np.max(area))
    plt.clim(np.min(area), np.max(area))
    plt.colorbar()

    l = ('imagem da área da banda entre ' + str(a) + ' cm-1' +
         str(b) + ' cm-1 \n ' + str(data['filename'])[:-4])
    plt.title(l)
    plt.show()


def mean(data, ini1, fim1):
    '''
    Gera e exibe uma imagem baseada na escala relativa ao espectro médio
    na região [ini1, fim1] cm⁻¹.

    Para cada pixel, ajusta a equação:
        z ≈ α · z̄ + c
    via mínimos quadrados, onde z̄ é o espectro médio da região. O coeficiente
    α (escala) mede o quanto aquele pixel se assemelha ao padrão médio em termos
    de intensidade. Pixels com α > 1 têm mais material do que a média; α < 1,
    menos material.

    Parâmetros
    ----------
    data : dict hsp
        Deve conter 'r', 'wn', 'dx', 'dy', 'sel', 'filename'.
    ini1 : float — limite inferior da região de análise em cm⁻¹
    fim1 : float — limite superior da região de análise em cm⁻¹

    Retorna
    -------
    dplot : ndarray (dx, dy)
        A imagem gerada, para uso posterior em subplots se necessário.
    '''
    import numpy as np
    import matplotlib.pyplot as plt

    # Seleciona a região espectral de interesse
    sel = np.logical_and(data['wn'] > ini1, data['wn'] < fim1)
    r1  = data['r'][:, sel]

    # Regressão linear de cada espectro contra o espectro médio da região:
    # sistema: y = α · média + c (mínimos quadrados)
    y     = r1[:, :].T
    xx    = np.vstack((r1.mean(axis=0), np.ones_like(r1[1, :]))).T
    alpha = np.linalg.lstsq(xx, y, rcond=-1)[0][0].T   # coeficiente α de cada pixel
    meansvalue = alpha

    # Reconstrói o mapa 2D
    dplot = np.zeros(data['dx'] * data['dy'])
    dplot[data['sel']] = meansvalue
    dplot = dplot.reshape(data['dx'], data['dy'])

    plt.pcolor(dplot, vmin=np.min(meansvalue), vmax=np.max(meansvalue))
    plt.clim(np.min(meansvalue), np.max(meansvalue))
    plt.colorbar()

    l = ('imagem da meanspc entre ' + str(ini1) + ' cm-1' +
         str(fim1) + ' cm-1 \n ' + str(data['filename'])[:-4])
    plt.title(l)
    return dplot


def pplot(data, nspc):
    '''
    Plota nspc espectros selecionados aleatoriamente da imagem.

    Útil para inspeção visual rápida da qualidade espectral após
    cada etapa de pré-processamento, permitindo verificar se bandas
    características estão presentes, se há artefatos, etc.

    Parâmetros
    ----------
    data : dict hsp
        Deve conter 'r' (espectros) e 'wn' (números de onda).
    nspc : int
        Número de espectros a plotar aleatoriamente.
    '''
    import numpy as np
    import matplotlib.pyplot as plt

    r = data['r']
    # Gera nspc índices aleatórios no intervalo [0, n_espectros)
    k = np.random.randint(0, r.shape[0], (nspc,), dtype='uint32')

    plt.figure()
    for i in k:
        plt.plot(data['wn'], r[i][:])   # plota o espectro de índice i
    plt.xlabel('Número de onda (cm⁻¹)')
    plt.show()


def emsc(datta, b):
    '''
    Gera e exibe uma imagem do coeficiente EMSC de índice b.

    Os coeficientes do modelo EMSC têm significado físico:
        - b=0 (a₀): escala da referência, proporcional à espessura/concentração
        - b=1,2,... (termos de baseline): variações lentas de fundo
        - demais: contribuição da parafina e da água

    O mapa do coeficiente a₀ é especialmente informativo para avaliar
    a homogeneidade da espessura da amostra.

    Parâmetros
    ----------
    datta : dict hsp
        Deve conter 'EMSC_coeff', 'sel', 'dx', 'dy', 'filename'.
    b : int
        Índice do coeficiente EMSC a mapear (0-based).
    '''
    import matplotlib.pyplot as plt
    import numpy as np

    # Extrai os coeficientes do índice b apenas para os pixels válidos (sel=True)
    ver = datta['EMSC_coeff'][datta['sel'], b]

    # Reconstrói o mapa 2D
    dplot = np.zeros(datta['dx'] * datta['dy'])
    dplot[datta['sel']] = ver
    dplot = dplot.reshape(datta['dx'], datta['dy'])

    plt.figure()
    plt.pcolor(dplot)
    plt.clim(np.min(ver), np.max(ver))
    plt.colorbar()

    l = 'imagem do coeficiente EMSC [' + str(b) + '] — ' + str(datta['filename'])[:-4]
    plt.title(l)
    plt.show()


def int_plt(datta, b):
    '''
    Versão interativa de intt(): exibe a imagem e permite clicar num pixel
    para visualizar seu espectro completo.

    O usuário clica na imagem para selecionar um pixel; o espectro daquele
    pixel é então exibido numa segunda figura. O processo repete até que
    o usuário clique num pixel próximo da borda (xx ≤ 5 ou yy ≤ 5).

    Útil para inspecionar espectros individuais de regiões de interesse
    identificadas visualmente no mapa de imagem.

    Parâmetros
    ----------
    datta : dict hsp
        Deve conter 'r', 'wn', 'dx', 'dy', 'sel'.
    b : float
        Número de onda em cm⁻¹ usado para gerar a imagem base.
    '''
    import matplotlib.pyplot as plt
    import numpy as np

    # Gera a imagem de intensidade (igual a intt, mas sem título)
    sel = datta['wn'] > b
    ver = datta['r'][:, sel]
    # Reconstrói a grade 3D original para acesso por coordenada (linha, coluna)
    dc  = datta['r'].reshape(datta['dx'], datta['dy'], -1)
    ver = ver[:, 0]
    data = ver.reshape(datta['dx'], datta['dy'])

    # Exibe a imagem e aguarda o primeiro clique do usuário
    plt.figure(1)
    plt.pcolor(data)
    x  = plt.ginput(1)    # captura 1 ponto de clique: retorna [(x, y)]
    xx = int(x[0][0])     # coordenada horizontal (coluna)
    yy = int(x[0][1])     # coordenada vertical (linha)
    print(xx, yy)

    # Loop interativo: continua até clicar próximo da borda (sinal de saída)
    while yy > 5 and xx > 5:
        plt.figure(1)
        plt.pcolor(data)
        x  = plt.ginput(1)
        # Nota: xx recebe y[1] e yy recebe x[0] para corrigir a orientação
        # da grade matplotlib vs. convenção matricial (linha, coluna)
        xx = int(x[0][1])
        yy = int(x[0][0])

        plt.close(2)
        plt.figure(2)
        plt.title(str(xx) + '  ' + str(yy))
        # dc[xx, yy, :] acessa o espectro do pixel na posição (linha=xx, coluna=yy)
        plt.plot(datta['wn'], dc[xx, yy, :].reshape(-1))
        print(xx, yy)


def pca(data, n, k=10):
    '''
    Gera e exibe uma imagem do score da componente principal n (PCn).

    A Análise de Componentes Principais (PCA) decompõe a matriz de espectros
    X (n_pixels × n_pontos) em:

        X = T · P^T + E

    onde:
        T — matriz de scores (n_pixels × k): coordenada de cada pixel no
            espaço das PCs; T[:, 0] = scores da PC1, T[:, 1] = PC2, etc.
        P — matriz de loadings (k × n_pontos): direção de cada PC no
            espaço espectral; indica quais bandas mais contribuem para a PC.
        E — resíduo (variância não explicada pelas k PCs).

    O mapa de scores da PCn revela a distribuição espacial do padrão
    espectral capturado por aquela componente. Componentes iniciais (PC1, PC2)
    costumam representar as principais variações físico-químicas da amostra
    (espessura, composição); componentes tardias tendem a capturar ruído.

    Algoritmo:
        1. Centraliza os dados: X_c = X - média_espectral
        2. Calcula a decomposição SVD da matriz de covariância: X_c = U·S·V^T
        3. Os scores são T = X_c · V (projeção dos espectros nas direções PCs)
        4. O score da PCn de cada pixel é mapeado de volta à grade 2D (dx × dy)

    Parâmetros
    ----------
    data : dict hsp
        Deve conter 'r', 'wn', 'dx', 'dy', 'sel', 'filename'.
    n : int
        Índice da componente principal a visualizar (1-based: n=1 → PC1).
    k : int, opcional (padrão=10)
        Número de PCs calculadas pelo SVD. Deve ser ≥ n. Aumentar k não
        altera os scores das primeiras PCs, mas aumenta o custo computacional.

    Retorna
    -------
    scores : ndarray (n_pixels_válidos,)
        Vetor dos scores da PCn para cada pixel válido.
    loadings : ndarray (n_pontos,)
        Loading da PCn: contribuição de cada número de onda para esta PC.
    var_pct : float
        Percentual da variância total explicado pela PCn.

    Exemplo
    -------
    >>> scores, loadings, var = sh.pc(data, n=1)
    # Exibe imagem da PC1 e retorna os scores para uso posterior
    '''
    import numpy as np
    import matplotlib.pyplot as plt

    r = data['r']                      # matriz (n_pixels × n_pontos)
    n_pixels, n_pts = r.shape

    # ── 1. Centralização ────────────────────────────────────────────────────
    # Subtrai a média espectral de cada variável (número de onda), tornando
    # o centroide da nuvem de pontos na origem.
    media = r.mean(axis=0)             # espectro médio global (n_pontos,)
    Xc    = r - media                  # matriz centrada

    # ── 2. SVD da matriz de covariância ─────────────────────────────────────
    # A decomposição SVD de Xc fornece U (scores não normalizados), S (valores
    # singulares) e Vt (loadings transposto). Usamos apenas as k primeiras PCs
    # para economizar memória e tempo de cálculo (SVD econômico).
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    # Restringe às k primeiras componentes
    k_eff = min(k, n_pixels, n_pts)
    U     = U[:, :k_eff]               # (n_pixels × k)
    S     = S[:k_eff]                  # (k,)
    Vt    = Vt[:k_eff, :]             # (k × n_pontos)

    # ── 3. Scores: projeção dos espectros nas direções PCs ───────────────────
    # T = Xc · V = U · diag(S)  →  cada coluna é os scores de uma PC
    scores_all = U * S                 # (n_pixels × k), broadcasting de S

    # ── 4. Variância explicada ───────────────────────────────────────────────
    var_total = (S ** 2).sum()
    var_pct   = (S[n - 1] ** 2 / var_total) * 100  # % da PCn (n é 1-based)

    print(f'PC{n}: {var_pct:.2f}% da variância total')

    # ── 5. Extrai os scores da PCn solicitada (n é 1-based → índice n-1) ────
    scores   = scores_all[:, n - 1]   # (n_pixels,)
    loadings = Vt[n - 1, :]           # (n_pontos,)

    # ── 6. Reconstrói o mapa 2D ──────────────────────────────────────────────
    # Posiciona cada score no pixel correto; pixels removidos (sel=False) ficam 0
    dplot = np.zeros(data['dx'] * data['dy'])
    dplot[data['sel']] = scores
    dplot = dplot.reshape(data['dx'], data['dy'])

    # ── 7. Visualização ──────────────────────────────────────────────────────
    # Usa colormap divergente (RdBu_r) pois scores são centrados em zero:
    # vermelho = score positivo (semelhante ao padrão da PC)
    # azul     = score negativo (oposto ao padrão da PC)
    vmax = np.abs(scores).max()        # escala simétrica em torno de zero

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    # Painel esquerdo: imagem dos scores
    im = axes[0].pcolor(dplot, cmap='RdBu_r', vmin=-vmax, vmax=vmax)
    plt.colorbar(im, ax=axes[0])
    axes[0].set_title(f'PC{n} scores ({var_pct:.2f}% var)\n{str(data["filename"])[:-4]}')
    axes[0].set_aspect('equal')

    # Painel direito: loading da PCn (quais bandas definem esta PC)
    axes[1].plot(data['wn'], loadings)
    axes[1].axhline(0, color='k', linewidth=0.5, linestyle='--')
    axes[1].set_xlabel('Número de onda (cm⁻¹)')
    axes[1].set_ylabel(f'Loading PC{n}')
    axes[1].set_title(f'Loading PC{n}')
    axes[1].invert_xaxis()             # convenção FTIR: eixo x decrescente

    plt.tight_layout()
    plt.show()

    return scores, loadings, var_pct