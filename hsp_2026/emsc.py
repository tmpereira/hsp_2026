"""
emsc.py — Módulo de Correção Multiplicativa de Espalhamento Estendida (EMSC)
=============================================================================
A EMSC (Extended Multiplicative Scatter Correction) é uma técnica de
pré-processamento espectral que remove variações físicas (espalhamento de luz,
diferenças de espessura, efeito de Mie) de espectros de infravermelho,
preservando apenas a informação química de interesse.

O modelo EMSC descreve cada espectro observado como:

    z(ν) = a₀ · r(ν) + a₁ + a₂·ν + a₃·ν² + ... + c₁·p₁(ν) + ... + c_n·hₙ(ν)

onde:
    r(ν)   — espectro de referência (média dos alvos)
    aᵢ     — coeficientes polinomiais de baseline
    pᵢ(ν)  — componentes principais da parafina (contaminante)
    hᵢ(ν)  — componentes principais da água (contaminante)

Após a correção, o espectro corrigido é:
    z_corr(ν) = (z(ν) - contribuições de interferentes) / a₀

Funções públicas:
    pca(p, k)                              — PCA auxiliar via SVD
    create_model_h2o(...)                  — monta modelo EMSC com parafina + água
    create_model(...)                      — monta modelo EMSC somente com parafina
    emsc_model_view(emsc)                  — visualiza os vetores do modelo
    emsc_fit(data, model)                  — aplica a correção EMSC nos dados
"""

import numpy as np
import matplotlib.pyplot as plt


# ─────────────────────────────────────────────────────────────────────────────
# FUNÇÃO AUXILIAR DE PCA
# ─────────────────────────────────────────────────────────────────────────────

def pca(p, k=10):
    """
    Realiza Análise de Componentes Principais (PCA) via Decomposição em
    Valores Singulares (SVD) da matriz de covariância.

    A PCA é usada internamente para extrair os padrões espectrais mais
    representativos da parafina e da água, que serão usados como regressores
    no modelo EMSC.

    Parâmetros
    ----------
    p : ndarray, shape (n_espectros, n_pontos)
        Matriz de espectros (linhas = amostras, colunas = variáveis espectrais).
    k : int, opcional (padrão=10)
        Número de componentes principais a reter.

    Retorna
    -------
    coeff : ndarray, shape (n_pontos, k)
        Loadings das k primeiras componentes principais (vetores de projeção).
        Cada coluna é um autovetor da matriz de covariância.
    perc : float
        Porcentagem da variância total explicada pelas k componentes retidas.
    """
    # Calcula a matriz de covariância entre variáveis espectrais (n_pontos x n_pontos)
    # Transpõe p para que cada coluna seja uma amostra, conforme exigido por np.cov
    pcov = np.cov(p.T)

    # SVD da matriz de covariância:
    #   coeff  — autovetores (loadings), colunas ortogonais
    #   latent — autovalores (variância explicada por cada PC), em ordem decrescente
    coeff, latent = np.linalg.svd(pcov, full_matrices=True)[:2]

    # Retém apenas as k primeiras componentes
    coeff = coeff[:, :k]

    # Calcula a porcentagem de variância explicada pelas k componentes
    perc = 100 * latent[:k].sum() / latent.sum()

    return (coeff, perc)


# ─────────────────────────────────────────────────────────────────────────────
# CONSTRUÇÃO DO MODELO EMSC — PARAFINA + ÁGUA
# ─────────────────────────────────────────────────────────────────────────────

def create_model_h2o(target, parafin, para_pcs, h2o, h2o_pcs, polyorder):
    """
    Constrói a matriz do modelo EMSC incluindo parafina e vapor de água
    como interferentes espectrais.

    A matriz do modelo (memsc) é montada coluna a coluna, na ordem:
        1. Espectro de referência (média dos alvos, normalizado)
        2. Polinômios de baseline (graus 0 até polyorder)
        3. Componentes principais da parafina (região 1300–1500 cm⁻¹)
        4. Componentes principais da água (região > 1350 cm⁻¹)

    Parâmetros
    ----------
    target : dict hsp
        Dados do tecido-alvo. Usado para calcular o espectro de referência.
        Deve conter 'r' (espectros) e 'wn' (números de onda).
    parafin : dict hsp
        Espectros puros de parafina para modelar a contaminação por este
        reagente de inclusão histológica.
    para_pcs : int
        Número de componentes principais da parafina a incluir no modelo.
    h2o : dict hsp
        Espectros puros de água para modelar a absorção do vapor atmosférico.
    h2o_pcs : int
        Número de componentes principais da água a incluir no modelo.
    polyorder : int
        Grau máximo do polinômio de correção de baseline.
        Ex: polyorder=2 inclui termos constante, linear e quadrático.

    Retorna
    -------
    emsc : dict
        Dicionário com o modelo EMSC contendo:
        'emsc_matrix'   — ndarray (n_pontos, n_regressores): matriz de regressão
        'legend'        — array de strings com nome de cada regressor
        'wn'            — vetor de números de onda cm⁻¹
        'perc_val_para' — % variância da parafina explicada pelos para_pcs
        'perc_val_h2o'  — % variância da água explicada pelos h2o_pcs
    """
    wn = target['wn']     # vetor de números de onda
    emsc = {}
    r = target['r']       # matriz de espectros do alvo

    # --- 1. Espectro de referência ---
    # Calcula a média espectral de todos os pixels do alvo
    memsc = r.mean(0)
    # Normaliza para norma unitária (vetor unitário) para estabilidade numérica
    memsc = memsc / (memsc @ memsc) ** 0.5
    memsc = memsc.reshape(-1, 1)   # shape: (n_pontos, 1)
    legend = ['mean_spc']

    # --- 2. Polinômios de baseline ---
    # Cria uma base polinomial no intervalo [-1, 1] (normalizado para estabilidade)
    # para modelar variações lentas de baseline causadas por espalhamento físico
    base = np.linspace(-1, 1, r.shape[1]).reshape(-1, 1)
    for i in range(polyorder + 1):
        legend.append('poly ' + str(i))
    polyorder = np.arange(0, polyorder + 1).reshape(1, -1)
    # Expande a base para ter uma coluna por grau: base^0, base^1, ..., base^n
    base = np.tile(base, (1, polyorder.shape[1]))
    base = base ** polyorder                     # shape: (n_pontos, polyorder+1)
    memsc = np.hstack((memsc, base))

    # --- 3. Componentes principais da parafina ---
    # Seleciona apenas a região espectral da parafina (1300–1500 cm⁻¹),
    # onde suas bandas C-H são dominantes; fora dela, zera para não interferir
    sel = (wn > 1300) & (wn < 1500)
    r = parafin['r'].copy()
    r[:, ~sel] = 0
    for i in range(para_pcs):
        legend.append('para ' + str(i + 1))
    ver = pca(r, para_pcs)
    emsc['perc_val_para'] = ver[1]    # salva % variância explicada
    memsc = np.hstack((memsc, ver[0]))

    # --- 4. Componentes principais da água ---
    # Seleciona apenas a região acima de 1350 cm⁻¹ onde as bandas de água
    # (vapor atmosférico) são mais pronunciadas no infravermelho médio
    sel = (wn > 1350)
    r = h2o['r'].copy()
    r[:, ~sel] = 0
    ver = pca(r, h2o_pcs)
    emsc['perc_val_h2o'] = ver[1]     # salva % variância explicada
    memsc = np.hstack((memsc, ver[0]))
    for i in range(h2o_pcs):
        legend.append('h2o ' + str(i + 1))

    # --- Monta o dicionário de saída ---
    emsc['emsc_matrix'] = memsc
    emsc['legend'] = np.array(legend)
    emsc['wn'] = wn
    return emsc


# ─────────────────────────────────────────────────────────────────────────────
# VISUALIZAÇÃO DO MODELO EMSC
# ─────────────────────────────────────────────────────────────────────────────

def emsc_model_view(emsc):
    """
    Visualiza cada coluna (regressor) da matriz do modelo EMSC em gráficos
    individuais, permitindo inspecionar a forma espectral de cada componente.

    Útil para validar se o modelo capturou corretamente os padrões de
    referência, baseline e interferentes (parafina, água).

    Parâmetros
    ----------
    emsc : dict
        Dicionário retornado por create_model() ou create_model_h2o(),
        contendo 'emsc_matrix', 'legend' e 'wn'.
    """
    memsc = emsc['emsc_matrix']   # shape: (n_pontos, n_regressores)
    legend = emsc['legend']        # rótulos de cada regressor
    wn = emsc['wn']               # números de onda cm⁻¹

    # Itera sobre cada regressor e plota seu perfil espectral
    for i in range(memsc.shape[1]):
        plt.figure()
        plt.plot(wn, memsc[:, i])
        plt.title(legend[i])      # título = nome do regressor


# ─────────────────────────────────────────────────────────────────────────────
# APLICAÇÃO DA CORREÇÃO EMSC
# ─────────────────────────────────────────────────────────────────────────────

def emsc_fit(data, model):
    """
    Aplica o modelo EMSC aos espectros de dados, corrigindo espalhamento
    e removendo interferentes (parafina, água, baseline).

    Algoritmo:
        1. Resolve o sistema linear:  XX · β = z(ν)ᵀ  →  β = lstsq(XX, zᵀ)
           Onde XX é a matriz do modelo e z são os espectros observados.
        2. O espectro corrigido é calculado subtraindo todas as contribuições
           dos interferentes (colunas 1..n, exceto a referência, coluna 0):
               z_corr = z - XX[:,1:] · β[1:,:]ᵀ
        3. Divide pelo coeficiente da referência (β[0,:]) para normalizar
           a intensidade absoluta:
               z_final = z_corr / β[0,:]

    Parâmetros
    ----------
    data : dict hsp
        Dicionário com os dados a corrigir. Deve conter 'r' (espectros).
        Modificado in-place com as chaves adicionadas.
    model : dict
        Modelo EMSC gerado por create_model() ou create_model_h2o().
        Deve conter 'emsc_matrix'.

    Retorna
    -------
    data : dict hsp
        O mesmo dicionário de entrada, atualizado com:
        'r'           — espectros corrigidos (n_espectros x n_pontos)
        'EMSC_model'  — matriz XX do modelo usada
        'EMSC_coeff'  — coeficientes β para cada espectro (n_espectros x n_regressores)
    """
    datar = data['r']             # shape: (n_espectros, n_pontos)
    XX = model['emsc_matrix']     # shape: (n_pontos, n_regressores)

    # Resolve mínimos quadrados: XX · β = datar.T
    # beta shape: (n_regressores, n_espectros)
    beta = np.linalg.lstsq(XX, datar.T, rcond=-1)[0]

    # Subtrai contribuição de todos os regressores EXCETO o espectro de referência
    # (coluna 0 = referência; colunas 1..n = baseline + interferentes)
    spccorr = datar - (XX[:, 1:].dot(beta[1:, :]).T)

    # Divide pelo coeficiente da referência (β[0,:]) para normalizar amplitude
    # tile repete o vetor para ter o mesmo shape dos espectros
    div = np.tile(beta[0, :].T.reshape(-1, 1), (1, spccorr.shape[1]))
    spccorr = spccorr / div

    # Armazena resultados no dicionário de dados
    data['EMSC_model'] = XX
    data['EMSC_coeff'] = beta.T   # transpõe para (n_espectros, n_regressores)
    data['r'] = spccorr
    return data


# ─────────────────────────────────────────────────────────────────────────────
# CONSTRUÇÃO DO MODELO EMSC — SOMENTE PARAFINA
# ─────────────────────────────────────────────────────────────────────────────

def create_model(target, parafin, para_pcs, polyorder):
    """
    Constrói a matriz do modelo EMSC incluindo apenas parafina como
    interferente espectral (sem correção de água).

    Versão simplificada de create_model_h2o(), usada quando o vapor de água
    não é uma fonte significativa de variação nos dados ou quando não há
    espectros de referência de água disponíveis.

    A matriz do modelo (memsc) é montada coluna a coluna, na ordem:
        1. Espectro de referência (média dos alvos, normalizado)
        2. Polinômios de baseline (graus 0 até polyorder)
        3. Componentes principais da parafina (região 1300–1500 cm⁻¹)

    Parâmetros
    ----------
    target : dict hsp
        Dados do tecido-alvo. Usado para calcular o espectro de referência.
        Deve conter 'r' (espectros) e 'wn' (números de onda).
    parafin : dict hsp
        Espectros puros de parafina.
    para_pcs : int
        Número de componentes principais da parafina a incluir no modelo.
    polyorder : int
        Grau máximo do polinômio de correção de baseline.

    Retorna
    -------
    emsc : dict
        Dicionário com o modelo EMSC contendo:
        'emsc_matrix'   — ndarray (n_pontos, n_regressores): matriz de regressão
        'legend'        — array de strings com nome de cada regressor
        'wn'            — vetor de números de onda cm⁻¹
        'perc_val_para' — % variância da parafina explicada pelos para_pcs
    """
    wn = target['wn']     # vetor de números de onda
    emsc = {}
    r = target['r']       # matriz de espectros do alvo

    # --- 1. Espectro de referência ---
    memsc = r.mean(0)
    memsc = memsc / (memsc @ memsc) ** 0.5
    memsc = memsc.reshape(-1, 1)
    legend = ['mean_spc']

    # --- 2. Polinômios de baseline ---
    base = np.linspace(-1, 1, r.shape[1]).reshape(-1, 1)
    for i in range(polyorder + 1):
        legend.append('poly ' + str(i))
    polyorder = np.arange(0, polyorder + 1).reshape(1, -1)
    base = np.tile(base, (1, polyorder.shape[1]))
    base = base ** polyorder
    memsc = np.hstack((memsc, base))

    # --- 3. Componentes principais da parafina ---
    # Restringe à região 1300–1500 cm⁻¹ (bandas C-H da parafina)
    sel = (wn > 1300) & (wn < 1500)
    r = parafin['r'].copy()
    r[:, ~sel] = 0
    for i in range(para_pcs):
        legend.append('para ' + str(i + 1))
    ver = pca(r, para_pcs)
    emsc['perc_val_para'] = ver[1]
    memsc = np.hstack((memsc, ver[0]))

    # --- Monta o dicionário de saída ---
    emsc['emsc_matrix'] = memsc
    emsc['legend'] = np.array(legend)
    emsc['wn'] = wn
    return emsc

