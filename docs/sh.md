> ⚠️ **Documentação em construção.** Algumas informações podem estar imprecisas ou incompletas.

> ⚠️ **Documentação em construção.** Algumas informações podem estar imprecisas ou incompletas.

# Módulo `sh.py` — Visualização de Imagens Hiperespectrais

O módulo `sh` (*show*) é a **janela visual** da biblioteca `hsp_2026`.
Ele transforma o vetor de espectros armazenado em `data['r']` em mapas
de imagem bidimensionais, revelando a distribuição espacial de compostos
químicos, espessuras e padrões multivariados na amostra.

---

## O princípio da reconstrução 2D

Após os filtros do `qt.py`, os pixels válidos ficam compactados em
`data['r']` (linhas = pixels, colunas = números de onda), mas a **posição
espacial original** de cada pixel é preservada na máscara `data['sel']`.

Todas as funções do `sh.py` usam o mesmo procedimento para reconstruir
a grade 2D:

```python
# 1. Calcula um valor escalar para cada pixel válido
valor = ...                                 # vetor de comprimento n_válidos

# 2. Cria uma grade com zeros (inclui pixels removidos)
dplot = np.zeros(data['dx'] * data['dy'])

# 3. Insere os valores nos pixels corretos
dplot[data['sel']] = valor

# 4. Reshapa para a grade 2D original
dplot = dplot.reshape(data['dx'], data['dy'])
```

Os pixels removidos aparecem como zero (preto na maioria dos colormaps),
delimitando visualmente a região de amostra válida.

---

## Funções do módulo

### `intt(data, b)` — Imagem de intensidade num ponto espectral

```python
sh.intt(data, b=1650)
```

**O que faz:**

Cria um mapa 2D onde o valor de cada pixel é a **absorbância no número
de onda mais próximo de `b`** cm⁻¹.

**Parâmetros:**

| Parâmetro | Tipo | Descrição |
|-----------|------|-----------|
| `data` | dict | Dicionário hsp com `'r'`, `'wn'`, `'dx'`, `'dy'`, `'sel'`, `'filename'` |
| `b` | float | Número de onda em cm⁻¹ |

**Bandas de referência típicas:**

| Banda (cm⁻¹) | Molécula | Aplicação |
|--------------|----------|-----------|
| 1650 | Amida I (C=O de proteínas) | Distribuição de proteínas |
| 1546 | Amida II (N–H de proteínas) | Confirmação de proteínas |
| 1740 | C=O de ésteres (lipídios) | Distribuição de lipídios |
| 2924 | CH₂ assimétrico (cadeia carbônica) | Lipídios insaturados |
| 1085 | P=O simétrico (ácidos nucleicos) | DNA/RNA |
| 1240 | P=O assimétrico | Ácidos nucleicos |

**Exemplo:**

```python
import hsp_2026.sh as sh

sh.intt(data, b=1650)   # mapa de distribuição de proteínas
sh.intt(data, b=1740)   # mapa de distribuição de lipídios
```

**Detalhes técnicos:**

A função usa `data['wn'] > b` e pega a **primeira coluna**, ou seja,
o ponto imediatamente acima de `b`. Para espectros com resolução de
4 cm⁻¹, o erro de posicionamento é ≤ 4 cm⁻¹, irrelevante na prática.

---

### `area(data, a, b)` — Imagem de área integrada de uma banda

```python
sh.area(data, a=1600, b=1700)
```

**O que faz:**

Calcula a **área integrada** (regra do trapézio) entre `a` e `b` cm⁻¹
para cada pixel e gera o mapa 2D dessas áreas:

$$
A_i = \int_{a}^{b} z_i(\tilde{\nu})\,d\tilde{\nu}
\;\approx\;
\sum_{k} \frac{z_i(\tilde{\nu}_k) + z_i(\tilde{\nu}_{k+1})}{2}
\,\Delta\tilde{\nu}
$$

**Parâmetros:**

| Parâmetro | Tipo | Descrição |
|-----------|------|-----------|
| `data` | dict | Dicionário hsp com `'r'`, `'wn'`, `'dx'`, `'dy'`, `'sel'`, `'filename'` |
| `a` | float | Limite inferior da banda (cm⁻¹) |
| `b` | float | Limite superior da banda (cm⁻¹) |

**Vantagem sobre `intt()`:**

Ao integrar sobre uma faixa, `area()` é mais resistente a:
- Ruído espectral pontual
- Pequenas variações no centro de banda entre pixels (deslocamento químico)
- Sobreposição de bandas vizinhas estreitas

Prefira `area()` quando a banda de interesse tem largura significativa
(ex.: Amida I 1600–1700 cm⁻¹, fingerprint 900–1800 cm⁻¹).

**Exemplo:**

```python
sh.area(data, a=1600, b=1700)   # mapa de Amida I integrada
sh.area(data, a=900,  b=1800)   # mapa da fingerprint completa (~ espessura total)
```

---

### `mean(data, ini1, fim1)` — Imagem de escala relativa ao espectro médio

```python
dplot = sh.mean(data, ini1=900, fim1=1800)
```

**O que faz:**

Para cada pixel, ajusta o modelo:

$$
\mathbf{z}_i \approx \alpha_i \cdot \bar{\mathbf{z}} + c_i
$$

via mínimos quadrados, onde \(\bar{\mathbf{z}}\) é o **espectro médio
global** da região selecionada. O coeficiente \(\alpha_i\) é mapeado
como imagem.

**Parâmetros:**

| Parâmetro | Tipo | Descrição |
|-----------|------|-----------|
| `data` | dict | Dicionário hsp com `'r'`, `'wn'`, `'dx'`, `'dy'`, `'sel'`, `'filename'` |
| `ini1` | float | Limite inferior da região (cm⁻¹) |
| `fim1` | float | Limite superior da região (cm⁻¹) |

**Retorna:** `dplot` — a imagem 2D como ndarray `(dx, dy)`

**Interpretação do coeficiente α:**

| Valor de α | Significado |
|-----------|-------------|
| α ≈ 1 | Pixel "médio" — forma e intensidade típicas |
| α > 1 | Pixel mais intenso que a média (amostra mais espessa) |
| α < 1 | Pixel menos intenso que a média (amostra mais fina) |
| α ≈ 0 | Pixel de fundo ou borda (quase sem sinal) |
| α < 0 | Artefato (forma espectral invertida — artefato de Mie) |

**Quando usar `mean()` em vez de `area()`?**

`mean()` normaliza pela **forma** do espectro médio, não apenas pela
magnitude. É mais informativo em amostras heterogêneas onde diferentes
regiões têm composições diferentes: enquanto `area()` mede apenas "quanto
absorve", `mean()` mede "quanto aquele pixel se parece com o padrão
típico da amostra, em termos de escala".

---

### `pplot(data, nspc)` — Espectros aleatórios

```python
sh.pplot(data, nspc=20)
```

**O que faz:**

Seleciona `nspc` pixels aleatoriamente e plota seus espectros sobrepostos
numa única figura. É uma ferramenta de **inspeção rápida** da qualidade
espectral.

**Parâmetros:**

| Parâmetro | Tipo | Descrição |
|-----------|------|-----------|
| `data` | dict | Dicionário hsp com `'r'` e `'wn'` |
| `nspc` | int | Número de espectros a plotar |

**Quando usar:**

```python
# Depois de cada etapa do pipeline, use pplot() para verificar
data = file.fsm('tecido.fsm')
sh.pplot(data, 10)     # espectros brutos — verifica a aquisição

data = prep.cut(data, 900, 1800)
sh.pplot(data, 10)     # após corte — confirma a faixa espectral

data = prep.snv(data)
sh.pplot(data, 10)     # após normalização — verifica se baselines estão corrigidas
```

**Dica:** Use `nspc=5` para uma visão rápida e `nspc=50` para uma
inspeção mais estatística da variabilidade da imagem.

---

### `emsc(data, b)` — Imagem de um coeficiente do modelo EMSC

```python
sh.emsc(data, b=0)   # mapa do coeficiente a₀
```

**O que faz:**

Mapeia o coeficiente de índice `b` do modelo EMSC como imagem 2D.
Requer que `emsc.emsc_fit()` tenha sido executado antes.

**Parâmetros:**

| Parâmetro | Tipo | Descrição |
|-----------|------|-----------|
| `data` | dict | Dicionário hsp com `'EMSC_coeff'`, `'sel'`, `'dx'`, `'dy'`, `'filename'` |
| `b` | int | Índice do coeficiente EMSC (0-based) |

**Índices dos coeficientes:**

| Índice `b` | Coeficiente | Mapa representa |
|-----------|-------------|----------------|
| 0 | \(a_0\) | Distribuição de "quantidade de amostra biológica" |
| 1 | \(b_0\) | Variação do offset de baseline pela imagem |
| 2 | \(b_1\) | Variação da inclinação da baseline |
| 3+ | \(c_k\) | Distribuição espacial de cada interferente (parafina, água) |

O mapa de \(a_0\) (índice 0) é o mais frequentemente usado porque
corresponde diretamente à **espessura efetiva de material biológico**
em cada pixel, já corrigida para parafina e variações de baseline.

---

### `int_plt(data, b)` — Inspeção interativa pixel a pixel

```python
sh.int_plt(data, b=1650)
```

**O que faz:**

Exibe o mapa de intensidade no número de onda `b` e entra em **modo
interativo**: você clica num pixel e o espectro completo daquele pixel
é exibido numa segunda figura. O ciclo repete até clicar próximo da
borda (x ≤ 5 ou y ≤ 5), que encerra o loop.

**Parâmetros:**

| Parâmetro | Tipo | Descrição |
|-----------|------|-----------|
| `data` | dict | Dicionário hsp com `'r'`, `'wn'`, `'dx'`, `'dy'`, `'sel'` |
| `b` | float | Número de onda (cm⁻¹) usado para gerar o mapa base |

**Fluxo de uso:**

```
1. sh.int_plt(data, 1650) → exibe o mapa de Amida I
2. Clique numa região de interesse na imagem
3. Uma segunda janela abre com o espectro daquele pixel
4. Clique em outra região para ver outro espectro
5. Clique na borda do mapa (canto inferior esquerdo) para sair
```

**Quando usar:**

`int_plt()` é ideal para **investigar pixels anômalos** identificados
nos mapas de imagem: um cluster com cor diferente na imagem de K-Means,
uma região com \(\alpha\) muito alto no mapa `mean()`, etc. Ao clicar
naquele pixel, você vê imediatamente se o espectro apresenta alguma
anomalia (saturação, artefato de Mie, contaminação).

**Nota técnica:**

Internamente, a função reconstrói o tensor 3D `(dx, dy, n_pontos)` para
acesso por coordenada `[linha, coluna]`. Há uma troca de eixos (`xx` e
`yy`) para corrigir a diferença entre a convenção de `plt.ginput()`
(retorna `x` = coluna, `y` = linha) e a convenção matricial (linha, coluna).

---

### `pc(data, n, k=10)` — Imagem e loading de componente principal (PCA)

```python
scores, loadings, var_pct = sh.pc(data, n=1)
scores, loadings, var_pct = sh.pc(data, n=2, k=15)
```

**O que faz:**

Realiza a **Análise de Componentes Principais** (PCA) via Decomposição
em Valores Singulares (SVD) e exibe:

- **Painel esquerdo**: mapa 2D dos scores da PC`n` (colormap RdBu_r divergente)
- **Painel direito**: loading espectral da PC`n` (contribuição de cada banda)

**Parâmetros:**

| Parâmetro | Tipo | Descrição |
|-----------|------|-----------|
| `data` | dict | Dicionário hsp com `'r'`, `'wn'`, `'dx'`, `'dy'`, `'sel'`, `'filename'` |
| `n` | int | Índice da PC a visualizar (**1-based**: `n=1` → PC1, `n=2` → PC2…) |
| `k` | int | Número de PCs calculadas pelo SVD (padrão = 10) |

**Retorna:**

| Nome | Tipo | Conteúdo |
|------|------|---------|
| `scores` | ndarray `(n_pixels,)` | Score da PC`n` para cada pixel válido |
| `loadings` | ndarray `(n_pontos,)` | Loading espectral da PC`n` |
| `var_pct` | float | % da variância total explicada pela PC`n` |

**O que é PCA?**

PCA encontra as **direções de máxima variância** no espaço espectral.
A matriz de espectros \(\mathbf{X}\) (pixels × números de onda) é decomposta como:

$$
\mathbf{X} = \mathbf{T}\,\mathbf{P}^T + \mathbf{E}
$$

onde:
- \(\mathbf{T}\) — **scores** (pixels × PCs): coordenada de cada pixel no espaço das PCs
- \(\mathbf{P}\) — **loadings** (PCs × números de onda): direção espectral de cada PC
- \(\mathbf{E}\) — resíduo não explicado

O algoritmo de `pc()`:

1. **Centralização**: \(\mathbf{X}_c = \mathbf{X} - \bar{\mathbf{z}}\) — subtrai o espectro médio
2. **SVD**: \(\mathbf{X}_c = \mathbf{U}\,\mathbf{S}\,\mathbf{V}^T\)
3. **Scores**: \(\mathbf{T} = \mathbf{U} \cdot \text{diag}(\mathbf{S})\)
4. **Variância da PC**\(n\): \(\text{var}_n = S_n^2 \,/\, \sum_i S_i^2\)

**Interpretação dos scores:**

O colormap **RdBu_r** (vermelho–branco–azul) é usado porque os scores
são centrados em zero:

- **Vermelho** (score positivo): pixels cujo espectro é **semelhante ao loading**
- **Azul** (score negativo): pixels cujo espectro é **oposto ao loading**
- **Branco** (score ≈ 0): pixels neutros em relação a essa PC

**Interpretação dos loadings:**

O loading espectral indica **quais bandas definem aquela PC**:

- Pico positivo no loading → banda associada às regiões vermelhas do mapa
- Pico negativo no loading → banda associada às regiões azuis do mapa

**Exemplo de interpretação:**

```
PC1 — 45% da variância
  Loading: pico positivo em 1650 (Amida I) e 1546 (Amida II)
           pico negativo em 1740 (C=O lipídico)
  → PC1 separa regiões ricas em proteína (vermelho) de regiões
    ricas em lipídio (azul)
```

**Quantas PCs calcular (`k`)?**

O parâmetro `k` controla o número de PCs no SVD econômico. Você sempre
pode usar `n < k`. Aumentar `k` não muda as PCs anteriores, mas aumenta
o tempo de cálculo. Valores típicos:

| Situação | k recomendado |
|----------|--------------|
| Exploração inicial | 5–10 |
| Análise detalhada | 10–20 |
| Comparação de muitas PCs | 20–30 |

**Exemplo completo:**

```python
import hsp_2026.sh as sh

# Explora as primeiras 3 PCs
for i in range(1, 4):
    scores, loadings, var = sh.pc(data, n=i)
    print(f'PC{i}: {var:.1f}% da variância')

# Salva os scores da PC1 para análise posterior
scores_pc1, _, _ = sh.pc(data, n=1)
```

---

## Comparativo das funções de imagem

| Função | Métrica por pixel | Precisa de EMSC? | Interativa? | Retorna valor? |
|--------|------------------|-----------------|-------------|----------------|
| `intt(data, b)` | Intensidade num ponto | Não | Não | Não |
| `area(data, a, b)` | Área integrada (trapézio) | Não | Não | Não |
| `mean(data, ini, fim)` | Coeficiente α (regressão) | Não | Não | Sim (`dplot`) |
| `emsc(data, b)` | Coeficiente EMSC | **Sim** | Não | Não |
| `int_plt(data, b)` | Intensidade num ponto | Não | **Sim** | Não |
| `pc(data, n, k)` | Score da PCn (SVD) | Não | Não | Sim (`scores, loadings, var`) |
| `pplot(data, nspc)` | — (espectros, não imagem) | Não | Não | Não |

---

## Workflow de exploração visual recomendado

Um pipeline típico de exploração de uma nova amostra usa várias funções
do `sh.py` em sequência:

```python
import hsp_2026.file as file
import hsp_2026.qt   as qt
import hsp_2026.prep as prep
import hsp_2026.sh   as sh

# ── 1. Carrega e verifica a qualidade ────────────────────────────────────
data = file.fsm('tecido.fsm')
sh.pplot(data, nspc=10)           # inspeção de espectros brutos

# ── 2. Controle de qualidade ─────────────────────────────────────────────
data = qt.otsu_area(data, 900, 1800)
sh.intt(data, 1650)               # confirma que apenas amostra está presente

# ── 3. Pré-processamento ─────────────────────────────────────────────────
data = prep.cut(data, 900, 1800)
data = prep.snv(data)
sh.pplot(data, nspc=10)           # verifica espectros após normalização

# ── 4. Exploração química ─────────────────────────────────────────────────
sh.intt(data, 1650)               # proteínas
sh.intt(data, 1740)               # lipídios
sh.area(data, 1600, 1700)         # Amida I integrada

# ── 5. PCA para padrões multivariados ────────────────────────────────────
sh.pc(data, n=1)                  # principal fonte de variação
sh.pc(data, n=2)                  # segunda fonte de variação
sh.pc(data, n=3)                  # terceira...

# ── 6. Inspeção pontual de regiões de interesse ──────────────────────────
sh.int_plt(data, 1650)            # clica em pixels específicos para ver espectros
```
