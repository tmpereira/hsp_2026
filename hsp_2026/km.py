'''
km.py — Módulo de Clusterização por K-Means para Imagens de Micro-FTIR
=======================================================================
Implementa análise de clusterização não supervisionada usando o algoritmo
K-Means (sklearn) aplicado a imagens hiperespectrais de FTIR.

O K-Means agrupa pixels espectralmente semelhantes em k clusters, permitindo
identificar regiões de composição química homogênea na amostra sem a necessidade
de um modelo supervisionado ou classes pré-definidas.

Estratégia de uso típica:
    1. Executar fit(data, k_max) para testar k de 2 até k_max clusters
    2. Visualizar o mapa com sh(data, k) para o k de interesse
    3. Visualizar os espectros médios por cluster com spc(data, k)

Funções públicas:
    fit(data, k, fold)           — K-Means para uma imagem (k = 2..k clusters)
    sh(data, k)                  — Mapa de clusters com barra de cores
    spc(data, k)                 — Espectros médios por cluster
    fit_common(data, k, fold)    — K-Means comum a múltiplas imagens
    sh2(data, k)                 — Mapa de clusters sem eixos (ideal para subplot)

Paleta de cores padrão (15 cores, compatível com até 14 clusters + fundo):
    0 — preto (pixels não selecionados / fundo)
    1 — vermelho,  2 — verde,  3 — azul,  4 — cinza,  5 — ciano,
    6 — roxo,  7 — verde-escuro,  8 — salmão,  9 — marfim,
    10 — azul-aço,  11 — oliva,  12 — pêssego,  13 — bege,  14 — ciano
'''


def fit(data, k, fold=15):
    '''
    Executa K-Means para todos os valores de k no intervalo [2, k].

    Para cada valor de k, o algoritmo é repetido 'fold' vezes com inicializações
    aleatórias diferentes (n_init=fold) e a melhor solução (menor inércia) é
    retida. Isso reduz a sensibilidade ao ponto de inicialização do K-Means.

    Os rótulos e centroides de TODOS os valores de k são armazenados
    simultaneamente no dicionário data, permitindo comparar soluções
    sem precisar re-executar o algoritmo.

    Parâmetros
    ----------
    data : dict hsp
        Dicionário com os dados espectrais. Deve conter 'r' (espectros
        pré-processados) e 'log'. Modificado in-place.
    k : int
        Número máximo de clusters. O ajuste é feito para k = 2, 3, ..., k.
    fold : int, opcional (padrão=15)
        Número de inicializações aleatórias do K-Means para cada valor de k.
        Valores maiores aumentam a robustez mas também o tempo de execução.

    Retorna
    -------
    data : dict hsp
        O mesmo dicionário de entrada, atualizado com:
        'km_label'      — ndarray (n_espectros, k-1): rótulos de cluster para
                          cada k testado. Coluna i contém os rótulos para k=i+2.
                          Os rótulos variam de 1 a k (não de 0 a k-1).
        'km_centroid'   — ndarray: centroides de todos os k testados empilhados.
        'km_k_centroid' — ndarray: identificador de k para cada centroide.
    '''
    from sklearn.cluster import KMeans
    import timeit
    import numpy as np

    # Pré-aloca arrays para guardar resultados de todos os k testados
    label     = np.ones((data['r'].shape[0], k - 1))  # shape: (n_espectros, k-1)
    centroid  = np.ones((1, data['r'].shape[1]))       # acumulador de centroides
    k_centroid = np.ones((1, 1))                       # identificador de k por centroide

    X = data['r']   # matriz de espectros: cada linha é um pixel
    print(str(X.shape[0]) + " spectra")
    tt = timeit.default_timer()   # marca o tempo total

    for i in range(2, k + 1):
        t = timeit.default_timer()
        # Executa K-Means com 'fold' inicializações; random_state=0 garante reprodutibilidade
        kmeans = KMeans(n_clusters=i, n_init=fold, random_state=0).fit(X)
        # Armazena rótulos na coluna i-2 (k=2 → coluna 0, k=3 → coluna 1, ...)
        # +1 para que os rótulos comecem em 1 (0 é reservado para o fundo)
        label[:, i - 2] = (kmeans.labels_ + 1)
        # Empilha os centroides deste k no acumulador
        centroid = np.vstack((centroid, kmeans.cluster_centers_))
        # Registra qual k gerou cada bloco de centroides
        k_centroid = np.vstack((k_centroid, np.tile(i, (i, 1))))
        print(str(i) + " cluster")
        print(str(np.round(timeit.default_timer() - t, decimals=2)) + " seconds")

    print("total time: " + str((np.round(timeit.default_timer() - tt, decimals=2))) + " seconds")

    # Salva resultados no dicionário de dados
    data['km_label']      = label
    data['km_centroid']   = centroid
    data['km_k_centroid'] = k_centroid

    # Registra operação no log
    linha = '\n kmeans de 2 até ' + str(k) + ' fazendo ' + str(fold) + ' repetições'
    print(linha, end='')
    data['log'] = np.char.add(data['log'], linha)
    return data


def sh(data, k):
    '''
    Exibe o mapa de clusterização K-Means para k clusters com barra de cores.

    O mapa é gerado redimensionando os rótulos de cluster do vetor de pixels
    selecionados (data['sel']) para a grade 2D (dx × dy) da imagem.
    Pixels não selecionados (data['sel'] == False) recebem rótulo 0 (preto).

    O gráfico usa proporção de aspecto correta com base nas dimensões
    da imagem, e inclui uma barra de cores lateral com as 14 cores disponíveis.

    Parâmetros
    ----------
    data : dict hsp
        Deve conter: 'km_label', 'sel', 'dx', 'dy'.
    k : int
        Número de clusters a visualizar (deve ter sido calculado em fit()).
    '''
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors

    # Reconstrói o mapa 2D: começa com zeros (fundo preto) e preenche os pixels válidos
    dplot = np.zeros((data['dx'] * data['dy']))
    dplot[data['sel']] = data['km_label'][:, k - 2]   # k-2 = índice da coluna para este k
    dplot = dplot.reshape(data['dx'], data['dy'])       # reshape para imagem 2D

    # Calcula a dimensão maior para normalizar o tamanho do gráfico
    n = np.max([data['dx'], data['dy']])

    # Paleta de cores discreta com 15 cores (0 = preto/fundo, 1-14 = clusters)
    colmap = [(0,0,0), (1,0,0), (0,1,0), (0,0,1), (0.41,0.41,0.41), (0,1,1),
              (0.58,0,0.82), (0,0.50,0), (0.98,0.50,0.44), (1, 1,0.87),
              (0.39,0.58,0.92), (0.50,0.50,0), (1,0.89,0.76), (0.96,0.96,0.86),
              (0,1,1)]
    cmap = colors.ListedColormap(colmap)
    boundaries = list(range(15))
    norm = colors.BoundaryNorm(boundaries, cmap.N, clip=True)

    # Posiciona o mapa com proporção de aspecto correta
    plt.axes([0.1, 0.1, (data['dy'] / n) * 0.7, (data['dx'] / n) * 0.7])
    plt.pcolor(dplot, cmap=cmap, norm=norm)

    # Adiciona barra de cores lateral mostrando as 14 cores disponíveis
    plt.axes([0.85, 0.1, 0.05, 0.75])
    plt.pcolor(np.arange(14)[:, None], cmap=cmap, norm=norm)


def spc(data, k):
    '''
    Plota o espectro médio de cada cluster para a solução de k clusters.

    Cada cluster é representado por sua cor correspondente na paleta padrão,
    facilitando a associação visual entre o mapa (sh) e os espectros.
    O espectro médio é calculado sobre todos os pixels pertencentes ao cluster.

    Parâmetros
    ----------
    data : dict hsp
        Deve conter: 'km_label', 'r' (espectros), 'wn' (números de onda).
    k : int
        Número de clusters cujos espectros serão plotados.
    '''
    import numpy as np
    import matplotlib.pyplot as plt

    # Paleta de cores (mesma de sh() para consistência visual)
    cmap = [(0,0,0), (1,0,0), (0,1,0), (0,0,1), (0.41,0.41,0.41), (0,1,1),
            (0.58,0,0.82), (0,0.50,0), (0.98,0.50,0.44), (1, 1,0.87),
            (0.39,0.58,0.92), (0.50,0.50,0), (1,0.89,0.76), (0.96,0.96,0.86),
            (0,1,1)]

    for i in list(range(k)):
        # Máscara booleana para selecionar pixels do cluster i+1
        sel = data['km_label'][:, k - 2] == i + 1
        # Plota o espectro médio do cluster com a cor correspondente
        plt.plot(data['wn'], np.mean(data['r'][sel, :], axis=0),
                 color=cmap[i], linewidth=2)


def fit_common(data, k, fold=10):
    '''
    Executa K-Means em múltiplas imagens concatenadas usando um espaço
    espectral comum.

    Esta função é útil quando se quer comparar diferentes amostras ou
    condições com os mesmos clusters, garantindo que a separação em grupos
    leve em conta a variabilidade de TODAS as imagens simultaneamente.

    Algoritmo:
        1. Concatena todos os espectros de todas as imagens em uma única
           matriz X, registrando a origem de cada espectro em 'imglabel'.
        2. Executa K-Means em X para k = 2..k.
        3. Redistribui os rótulos de cada espectro de volta para o
           dicionário da imagem de origem.

    Parâmetros
    ----------
    data : list of dict hsp
        Lista de dicionários hsp, um por imagem. Cada dicionário deve
        conter 'r' (espectros pré-processados) e 'log'.
    k : int
        Número máximo de clusters.
    fold : int, opcional (padrão=10)
        Número de inicializações aleatórias para cada valor de k.

    Retorna
    -------
    data : list of dict hsp
        A mesma lista de entrada, com cada dicionário atualizado com:
        'km_label'      — rótulos para os pixels daquela imagem
        'km_centroid'   — centroides comuns a todas as imagens
        'km_k_centroid' — identificador de k por centroide
    '''
    from sklearn.cluster import KMeans
    import timeit
    import numpy as np

    # Concatena os espectros de todas as imagens e registra a origem de cada um
    X = data[0]['r']
    imglabel = np.zeros((data[0]['r'].shape[0], 1))   # imagem 0 → rótulo 0
    for i in range(1, len(data)):
        print(i)
        X = np.vstack((X, data[i]['r']))
        # Marca cada espectro com o índice da imagem de origem
        imglabel = np.vstack((imglabel, int(i) * np.ones((data[i]['r'].shape[0], 1))))

    # Pré-aloca arrays para todos os k testados
    label     = np.ones((X.shape[0], k - 1))
    centroid  = np.ones((1, X.shape[1]))
    k_centroid = np.ones((1, 1))

    print(str(X.shape[0]) + " spectra")
    for i in range(2, k + 1):
        t = timeit.default_timer()
        kmeans = KMeans(n_clusters=i, n_init=fold, random_state=0).fit(X)
        label[:, i - 2] = (kmeans.labels_ + 1)
        centroid = np.vstack((centroid, kmeans.cluster_centers_))
        k_centroid = np.vstack((k_centroid, np.tile(i, (i, 1))))
        print(str(i) + " cluster")
        print(str(np.round(timeit.default_timer() - t, decimals=2)) + " seconds")

    # Distribui os rótulos de volta para cada imagem original
    j = 0
    linha = '\n common kmeans de 2 até ' + str(k) + ' fazendo ' + str(fold) + ' repetições'
    for i in data:
        # Filtra as linhas do array global que pertencem à imagem j
        sel = imglabel == j
        i['km_label']      = label[sel.reshape(-1), :]
        i['km_centroid']   = centroid
        i['km_k_centroid'] = k_centroid
        i['log'] = np.char.add(i['log'], linha)
        j += 1

    print(linha, end='')
    return data


def sh2(data, k):
    '''
    Exibe o mapa de clusterização K-Means sem eixos e sem barra de cores.

    Versão simplificada de sh(), ideal para uso em subplots onde múltiplos
    mapas são exibidos lado a lado. Remove todos os elementos decorativos
    (eixos, ticks, colorbar) para maximizar o espaço da imagem.

    Parâmetros
    ----------
    data : dict hsp
        Deve conter: 'km_label', 'sel', 'dx', 'dy'.
    k : int
        Número de clusters a visualizar.
    '''
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors

    # Reconstrói o mapa 2D com zeros no fundo e rótulos nos pixels válidos
    dplot = np.zeros((data['dx'] * data['dy']))
    dplot[data['sel']] = data['km_label'][:, k - 2]
    dplot = dplot.reshape(data['dx'], data['dy'])

    # Paleta e normalização idênticas a sh()
    colmap = [(0,0,0), (1,0,0), (0,1,0), (0,0,1), (0.41,0.41,0.41), (0,1,1),
              (0.58,0,0.82), (0,0.50,0), (0.98,0.50,0.44), (1, 1,0.87),
              (0.39,0.58,0.92), (0.50,0.50,0), (1,0.89,0.76), (0.96,0.96,0.86),
              (0,1,1)]
    cmap = colors.ListedColormap(colmap)
    boundaries = list(range(15))
    norm = colors.BoundaryNorm(boundaries, cmap.N, clip=True)

    plt.pcolor(dplot, cmap=cmap, norm=norm)
    plt.axis('off')   # remove eixos para visualização limpa em subplots