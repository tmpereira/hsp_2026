> ⚠️ **Documentação em construção.** Algumas informações podem estar imprecisas ou incompletas.

> ⚠️ **Documentação em construção.** Algumas informações podem estar imprecisas ou incompletas.

# Módulo `km.py` — Clusterização K-Means

O módulo `km` implementa **análise de agrupamento não supervisionada**
usando o algoritmo K-Means aplicado aos espectros da imagem hiperespectral.
O objetivo é identificar automaticamente regiões de composição química
homogênea sem a necessidade de classes pré-definidas ou espectros de referência.

---

## O que é K-Means?

K-Means é um algoritmo de aprendizado de máquina não supervisionado que
agrupa \(N\) amostras em \(k\) clusters, minimizando a **soma das
distâncias quadradas** de cada ponto ao centroide do seu cluster:

$$
\min_{\{C_j\}} \sum_{j=1}^{k} \sum_{\mathbf{x}_i \in C_j}
\|\mathbf{x}_i - \boldsymbol{\mu}_j\|^2
$$

onde \(\boldsymbol{\mu}_j\) é o centroide (espectro médio) do cluster \(j\).

No contexto de Micro-FTIR:
- Cada **pixel** é um ponto em um espaço de \(n\) dimensões (uma por número de onda)
- Pixels espectralmente semelhantes → mesmo cluster
- Cada cluster corresponde a uma **região de composição química homogênea**

### K-Means em imagens hiperespectrais

```
Imagem: 100 × 100 pixels × 900 números de onda
  ↓
Matriz X: (10.000 × 900)   ← entrada do K-Means
  ↓
K-Means k=4
  ↓
Rótulos: vetor (10.000,) com valores 1, 2, 3 ou 4
  ↓
Reshape + reconstrução do mapa 2D → imagem colorida 100 × 100
```

### Limitações do K-Means

| Limitação | O que fazer |
|-----------|-------------|
| Número de clusters \(k\) deve ser definido pelo usuário | Use `fit(data, k_max)` para testar vários `k` de uma vez |
| Sensível à inicialização aleatória | Use `fold ≥ 10` repetições |
| Pressupõe clusters esféricos no espaço espectral | Pode ser ineficaz para grupos muito alongados ou irregulares |
| Não escala bem sem pré-processamento | Sempre aplique `prep.snv()` ou EMSC antes |

---

## A paleta de cores

Todos os mapas do módulo `km` usam uma **paleta discreta de 15 cores**:

| Rótulo | Cor | RGB aproximado |
|--------|-----|----------------|
| 0 | Preto (fundo) | (0, 0, 0) |
| 1 | Vermelho | (1, 0, 0) |
| 2 | Verde | (0, 1, 0) |
| 3 | Azul | (0, 0, 1) |
| 4 | Cinza | (0.41, 0.41, 0.41) |
| 5 | Ciano | (0, 1, 1) |
| 6 | Roxo | (0.58, 0, 0.82) |
| 7 | Verde-escuro | (0, 0.50, 0) |
| 8 | Salmão | (0.98, 0.50, 0.44) |
| 9 | Marfim | (1, 1, 0.87) |
| 10 | Azul-aço | (0.39, 0.58, 0.92) |
| 11 | Oliva | (0.50, 0.50, 0) |
| 12 | Pêssego | (1, 0.89, 0.76) |
| 13 | Bege | (0.96, 0.96, 0.86) |
| 14 | Ciano claro | (0, 1, 1) |

O rótulo 0 (preto) é sempre reservado para os **pixels removidos**
pelos filtros do `qt.py`. Por isso os clusters recebem rótulos 1, 2, 3…
e nunca 0.

---

## Funções do módulo

### `fit(data, k, fold=15)` — Executa K-Means para múltiplos valores de k

```python
data = km.fit(data, k=6)        # testa k = 2, 3, 4, 5, 6
data = km.fit(data, k=8, fold=20)  # com 20 inicializações por k
```

**O que faz:**

Executa o K-Means para **todos os valores de k no intervalo `[2, k]`**,
de uma vez só. Para cada k, repete o ajuste `fold` vezes com inicializações
aleatórias diferentes e retém a solução com a menor inércia.

**Por que testar vários k?**

Não existe um critério objetivo universal para escolher o "melhor" k.
A estratégia recomendada é:
1. Rodar `fit(data, k_max)` uma única vez (custa tempo, mas calcula tudo)
2. Visualizar os mapas com `sh(data, k)` para cada k
3. Interpretar os espectros médios com `spc(data, k)`
4. Escolher o k que produz clusters biologicamente interpretáveis

**Parâmetros:**

| Parâmetro | Tipo | Descrição |
|-----------|------|-----------|
| `data` | dict | Dicionário hsp com `'r'` e `'log'` |
| `k` | int | Número máximo de clusters (calcula de 2 até k) |
| `fold` | int | Número de inicializações aleatórias por k (padrão = 15) |

**Retorna:** `data` atualizado com:

| Chave | Shape | Conteúdo |
|-------|-------|---------|
| `data['km_label']` | `(n_pixels, k-1)` | Coluna `i` contém os rótulos para k = i+2 |
| `data['km_centroid']` | `(total_centroides, n_pontos)` | Centroides de todos os k empilhados |
| `data['km_k_centroid']` | `(total_centroides, 1)` | Índice k de cada centroide |

**Acesso aos rótulos:**

```python
# k=2 → coluna 0, k=3 → coluna 1, k=4 → coluna 2...
rotulos_k4 = data['km_label'][:, 4 - 2]   # coluna índice k-2
```

As funções `sh()` e `spc()` fazem esse cálculo automaticamente quando
você passa o parâmetro `k`.

**Saída no terminal durante a execução:**

```
10000 spectra
2 cluster
0.84 seconds
3 cluster
1.12 seconds
4 cluster
1.35 seconds
...
total time: 8.21 seconds
```

**Dica sobre o `fold`:**

Valores maiores de `fold` aumentam as chances de encontrar o mínimo
global do K-Means, mas aumentam proporcionalmente o tempo de execução.
Para exploração inicial use `fold=5`; para resultados finais use `fold=20`.

---

### `sh(data, k)` — Mapa de clusters com barra de cores

```python
km.sh(data, k=4)
```

**O que faz:**

Reconstrói o mapa 2D de clusters para a solução de `k` grupos e exibe
com a paleta discreta de 15 cores e uma barra de cores lateral.

**Parâmetros:**

| Parâmetro | Tipo | Descrição |
|-----------|------|-----------|
| `data` | dict | Dicionário hsp com `'km_label'`, `'sel'`, `'dx'`, `'dy'` |
| `k` | int | Número de clusters a visualizar |

**Como interpretar o mapa:**

- Cada cor representa um grupo de pixels espectralmente semelhantes
- Pixels pretos (rótulo 0) são os removidos pelo `qt.py`
- Regiões de mesma cor têm composição química similar

**Proporção de aspecto:**

A função calcula automaticamente a proporção correta entre largura e
altura com base em `data['dx']` e `data['dy']`, evitando distorção da
imagem.

---

### `spc(data, k)` — Espectros médios por cluster

```python
km.spc(data, k=4)
```

**O que faz:**

Para cada cluster de 1 a `k`, calcula o **espectro médio** de todos os
pixels pertencentes ao cluster e plota sobrepostos na mesma figura.
Cada espectro usa a cor correspondente da paleta, facilitando a
associação com o mapa gerado por `sh()`.

**Parâmetros:**

| Parâmetro | Tipo | Descrição |
|-----------|------|-----------|
| `data` | dict | Dicionário hsp com `'km_label'`, `'r'`, `'wn'` |
| `k` | int | Número de clusters |

**Como interpretar os espectros médios:**

Os espectros médios são a "assinatura química" de cada cluster. Compare
as diferenças entre eles para atribuir identidade biológica:

| Diferença espectral | Interpretação possível |
|--------------------|-----------------------|
| Cluster com Amida I/II mais intensas (1650, 1546 cm⁻¹) | Região rica em proteína |
| Cluster com C=O éster intenso (1740 cm⁻¹) | Região rica em lipídio |
| Cluster com P=O intenso (1085, 1240 cm⁻¹) | Região rica em ácido nucleico (núcleo celular) |
| Cluster com espectro de baixa intensidade | Possivelmente fundo residual |

**Exemplo de workflow visual:**

```python
# Executa K-Means
km.fit(data, k=5)

# Visualiza o mapa para k=4
km.sh(data, k=4)

# Em outra figura, plota os espectros médios dos 4 clusters
km.spc(data, k=4)
# → Identifica: cluster 1 (vermelho) = epitelio (rico em proteína)
#               cluster 2 (verde)    = estroma (colágeno)
#               cluster 3 (azul)    = núcleos (alto P=O)
#               cluster 4 (cinza)   = borda/fundo residual
```

---

### `fit_common(data_list, k, fold=10)` — K-Means comum a múltiplas imagens

```python
# data é uma lista de dicionários hsp
dados = [data1, data2, data3]
dados = km.fit_common(dados, k=5)

# Agora cada dicionário tem km_label com os rótulos baseados
# no espaço espectral comum das 3 imagens
km.sh(dados[0], k=4)
km.sh(dados[1], k=4)
km.sh(dados[2], k=4)
```

**O que faz:**

Concatena os espectros de **todas as imagens** numa única matriz, executa
o K-Means nesse espaço espectral comum, e redistribui os rótulos de volta
para cada imagem original.

**Por que isso é necessário?**

Se você rodar `km.fit()` separadamente em cada imagem, o K-Means de
cada uma vai definir seus próprios clusters independentemente. Isso
significa que o "Cluster 1" da imagem A pode não corresponder
quimicamente ao "Cluster 1" da imagem B — a mesma cor pode representar
tecidos diferentes nas duas imagens.

Com `fit_common()`, todos os dados são analisados no mesmo espaço:
os clusters têm **significado consistente** entre todas as imagens,
permitindo comparação direta.

**Parâmetros:**

| Parâmetro | Tipo | Descrição |
|-----------|------|-----------|
| `data` | list de dict | Lista de dicionários hsp, um por imagem |
| `k` | int | Número máximo de clusters |
| `fold` | int | Número de inicializações por k (padrão = 10) |

**Retorna:** a mesma lista, com cada dicionário atualizado com `'km_label'`,
`'km_centroid'` e `'km_k_centroid'`.

**Quando usar `fit_common()` em vez de `fit()`:**

| Situação | Recomendação |
|----------|-------------|
| Uma única imagem | `fit()` |
| Várias imagens do mesmo experimento (ex.: controle vs. tratado) | `fit_common()` |
| Réplicas biológicas que devem ser comparadas | `fit_common()` |
| Imagens de diferentes experimentos sem relação entre si | `fit()` separado |

**Requisito:** Todas as imagens devem ter o **mesmo vetor de números de onda**
(mesma faixa espectral, mesma resolução). Garanta isso aplicando
`prep.cut(data, ini, fim)` com os mesmos valores em todas antes de
chamar `fit_common()`.

---

### `sh2(data, k)` — Mapa de clusters para subplots (sem eixos)

```python
# Comparação de vários k lado a lado
fig, axes = plt.subplots(1, 4, figsize=(16, 4))
for i, k_val in enumerate([2, 3, 4, 5]):
    plt.sca(axes[i])
    km.sh2(data, k=k_val)
    axes[i].set_title(f'k = {k_val}')
plt.tight_layout()
```

**O que faz:**

Idêntico a `sh()`, mas **sem eixos e sem barra de cores**. Útil para
compor figuras com múltiplos mapas lado a lado.

**Parâmetros:**

| Parâmetro | Tipo | Descrição |
|-----------|------|-----------|
| `data` | dict | Dicionário hsp com `'km_label'`, `'sel'`, `'dx'`, `'dy'` |
| `k` | int | Número de clusters a visualizar |

**Diferença em relação a `sh()`:**

| Aspecto | `sh()` | `sh2()` |
|---------|--------|---------|
| Eixos | Visíveis | Ocultos (`plt.axis('off')`) |
| Barra de cores | Presente | Ausente |
| Proporção de aspecto | Calculada automaticamente | Depende do axes atual |
| Uso ideal | Figura standalone | Subplots comparativos |

---

## Comparativo das funções

| Função | Modifica `data`? | Entrada | Saída principal |
|--------|-----------------|---------|----------------|
| `fit(data, k)` | **Sim** | Um dict | `data['km_label']`, centroides |
| `sh(data, k)` | Não | Um dict | Mapa colorido com barra |
| `spc(data, k)` | Não | Um dict | Espectros médios sobrepostos |
| `fit_common(list, k)` | **Sim** | Lista de dicts | `km_label` em cada dict |
| `sh2(data, k)` | Não | Um dict | Mapa colorido sem eixos |

---

## Workflow completo recomendado

### Uma única imagem

```python
import hsp_2026.file as file
import hsp_2026.prep as prep
import hsp_2026.qt   as qt
import hsp_2026.sh   as sh
import hsp_2026.km   as km

# ── 1. Carrega e pré-processa ────────────────────────────────────────────
data = file.fsm('tecido.fsm')
data = qt.otsu_area(data, 900, 1800)
data = prep.cut(data, 900, 1800)
data = prep.snv(data)

# ── 2. Exploração inicial com PCA ────────────────────────────────────────
sh.pc(data, n=1)   # quantas fontes de variação existem?
sh.pc(data, n=2)

# ── 3. K-Means para vários k ─────────────────────────────────────────────
data = km.fit(data, k=6, fold=15)

# ── 4. Inspeciona os mapas e espectros para cada k ───────────────────────
for k_val in [3, 4, 5, 6]:
    plt.figure(figsize=(10, 4))
    plt.subplot(1, 2, 1)
    km.sh2(data, k=k_val)
    plt.title(f'Mapa k={k_val}')
    plt.subplot(1, 2, 2)
    km.spc(data, k=k_val)
    plt.title(f'Espectros médios k={k_val}')
    plt.tight_layout()
```

### Múltiplas imagens com espaço espectral comum

```python
# Carrega e pré-processa todas as imagens com os mesmos parâmetros
dados = []
for arquivo in ['controle.fsm', 'tratado.fsm', 'tratado2.fsm']:
    d = file.fsm(arquivo)
    d = qt.otsu_area(d, 900, 1800)
    d = prep.cut(d, 900, 1800)
    d = prep.snv(d)
    dados.append(d)

# K-Means no espaço comum
dados = km.fit_common(dados, k=5, fold=15)

# Visualiza comparativamente (mesma paleta, mesmos clusters)
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
titulos = ['Controle', 'Tratado 1', 'Tratado 2']
for i, (d, titulo) in enumerate(zip(dados, titulos)):
    plt.sca(axes[i])
    km.sh2(d, k=4)
    axes[i].set_title(titulo)
plt.tight_layout()
```
