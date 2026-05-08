> ⚠️ **Documentação em construção.** Algumas informações podem estar imprecisas ou incompletas.

> ⚠️ **Documentação em construção.** Algumas informações podem estar imprecisas ou incompletas.

# Módulo `emsc.py` — Correção Multiplicativa de Espalhamento Estendida

O módulo `emsc` implementa a **EMSC** (*Extended Multiplicative Scatter
Correction*), uma técnica de pré-processamento espectral desenvolvida
especificamente para espectroscopia de infravermelho de tecidos biológicos.

---

## O problema do espalhamento em Micro-FTIR

Quando a luz infravermelha interage com uma amostra biológica heterogênea,
além de ser **absorvida** pelas moléculas (o sinal químico que queremos),
ela também é **espalhada**. Esse espalhamento físico distorce o espectro
de várias formas:

| Efeito físico | Causa | Distorção no espectro |
|---------------|-------|----------------------|
| Espessura variável | Cortes histológicos não perfeitamente uniformes | Escala multiplicativa diferente em cada pixel |
| Efeito de Mie | Células com tamanho comparável ao comprimento de onda | Linha de base ondulada (oscilações de longa escala) |
| Inclusão em parafina | Resíduo de parafina não removido na histologia | Bandas C–H sobrepostas (1300–1500 cm⁻¹) |
| Vapor de água | Purga incompleta do caminho óptico | Bandas de rotação-vibração da água |

Sem correção, esses efeitos **dominam a variância** dos dados e
dificultam ou impossibilitam a análise multivariada (K-Means, PCA).

---

## O modelo matemático da EMSC

A EMSC modela cada espectro observado \(z(\tilde{\nu})\) como soma de
contribuições físicas e interferentes:

$$
z(\tilde{\nu}) = a_0 \cdot r(\tilde{\nu})
  + \underbrace{b_0 + b_1\tilde{\nu} + b_2\tilde{\nu}^2 + \cdots}_{\text{baseline polinomial}}
  + \underbrace{c_1 p_1(\tilde{\nu}) + \cdots + c_m p_m(\tilde{\nu})}_{\text{parafina}}
  + \underbrace{d_1 h_1(\tilde{\nu}) + \cdots + d_n h_n(\tilde{\nu})}_{\text{água}}
  + \varepsilon(\tilde{\nu})
$$

onde:
- \(r(\tilde{\nu})\) — espectro de referência (média espectral do alvo)
- \(a_0\) — escala da referência (proporcional à espessura × concentração)
- \(b_i\) — coeficientes polinomiais de baseline
- \(p_i(\tilde{\nu})\) — componentes principais dos espectros de parafina
- \(h_i(\tilde{\nu})\) — componentes principais dos espectros de água
- \(\varepsilon\) — resíduo (variação química real + ruído)

O modelo é **ajustado por mínimos quadrados** para cada pixel,
estimando os coeficientes \(a_0, b_i, c_i, d_i\).

### A correção

Após o ajuste, o espectro corrigido é calculado **subtraindo os
interferentes e dividindo pela escala**:

$$
z_{\text{corr}}(\tilde{\nu}) = \frac{z(\tilde{\nu}) - \text{(baseline + parafina + água)}}{a_0}
$$

O resultado é um espectro normalizado que representa apenas a
**composição química relativa** do pixel, livre de efeitos físicos.

---

## Fluxo de trabalho completo

```
1. Coletar espectros de referência:
   • Espectros do tecido-alvo (para a referência r)
   • Espectros puros de parafina (para modelar a contaminação)
   • [Opcional] Espectros de vapor de água

2. Construir o modelo:
   emsc.create_model()       — sem água
   emsc.create_model_h2o()   — com água

3. Validar o modelo:
   emsc.emsc_model_view()    — visualiza cada componente

4. Aplicar a correção:
   emsc.emsc_fit()           — corrige todos os pixels

5. Controle de qualidade pós-correção:
   rg.emsc(data, a=0)        — histograma de a₀
   qt.otsu_emsc(data)        — remove pixels com a₀ baixo
```

---

## Funções do módulo

### `pca(p, k=10)` — PCA auxiliar para extração de padrões

```python
loadings, variancia_pct = emsc.pca(espectros_parafina, k=3)
```

**O que faz:**

Realiza PCA via SVD da **matriz de covariância** e retorna os loadings
(autovetores) das `k` primeiras componentes. Esta função é usada
**internamente** por `create_model()` e `create_model_h2o()` para
extrair os padrões espectrais mais representativos da parafina e da água.

**Por que PCA para modelar os interferentes?**

Em vez de usar um único espectro de parafina como regressor, a EMSC usa
as `k` primeiras PCs dos espectros de parafina. Isso é mais robusto porque:
- A parafina pode ter variações entre diferentes lotes ou cortes
- As PCs capturam toda a variabilidade conhecida do interferente
- Com 2–3 PCs de parafina, tipicamente ≥ 99% da variância é coberta

**Parâmetros:**

| Parâmetro | Tipo | Descrição |
|-----------|------|-----------|
| `p` | ndarray `(n, m)` | Matriz de espectros (linhas = amostras, colunas = números de onda) |
| `k` | int | Número de PCs a reter (padrão = 10) |

**Retorna:**

| Nome | Shape | Conteúdo |
|------|-------|---------|
| `loadings` | `(n_pontos, k)` | Autovetores — um loading por coluna |
| `variancia_pct` | float | % da variância total explicada pelas `k` PCs |

**Exemplo:**

```python
# Verifica quantas PCs são suficientes para a parafina
for k in [1, 2, 3, 5]:
    _, var = emsc.pca(parafina['r'], k=k)
    print(f'k={k}: {var:.1f}% da variância')
# k=1: 85.3%
# k=2: 97.6%
# k=3: 99.1%   ← geralmente suficiente
```

---

### `create_model(target, parafin, para_pcs, polyorder)` — Modelo EMSC sem água

```python
model = emsc.create_model(
    target    = data_tecido,
    parafin   = data_parafina,
    para_pcs  = 3,
    polyorder = 2
)
```

**O que faz:**

Constrói a matriz do modelo EMSC montando coluna a coluna os seguintes
regressores:

```
Coluna 0:        Espectro de referência (média dos alvos, normalizado)
Colunas 1..p+1:  Polinômios de baseline (grau 0, 1, …, polyorder)
Colunas seguintes: PCs da parafina (para_pcs colunas)
```

**Parâmetros:**

| Parâmetro | Tipo | Descrição |
|-----------|------|-----------|
| `target` | dict | Dados do tecido-alvo — usados para calcular o espectro de referência |
| `parafin` | dict | Espectros puros de parafina |
| `para_pcs` | int | Número de PCs da parafina a incluir |
| `polyorder` | int | Grau máximo do polinômio de baseline |

**Retorna:** dicionário `model` com:

| Chave | Conteúdo |
|-------|---------|
| `'emsc_matrix'` | ndarray `(n_pontos, n_regressores)` — a matriz de regressão |
| `'legend'` | array de strings com o nome de cada coluna |
| `'wn'` | vetor de números de onda |
| `'perc_val_para'` | % da variância da parafina coberta pelos `para_pcs` |

**Detalhes técnicos importantes:**

*Espectro de referência normalizado:*

$$
\mathbf{r} \leftarrow \frac{\bar{\mathbf{z}}}{\|\bar{\mathbf{z}}\|_2}
$$

A normalização pela norma euclidiana garante que o coeficiente \(a_0\)
tenha escala consistente entre experimentos diferentes.

*Polinômios em base normalizada:*

O polinômio é construído sobre uma grade \(\lambda \in [-1, +1]\)
(e não sobre os valores brutos de cm⁻¹). Isso evita problemas numéricos
de mal-condicionamento que ocorrem quando se usa \(\tilde{\nu}^2\) com
valores na casa de \(10^6\).

*Região da parafina:*

As PCs da parafina são calculadas apenas na janela **1300–1500 cm⁻¹**
(bandas C–H deformação de parafina), com zeros fora dessa janela.
Isso evita que as PCs da parafina capturem variações do tecido em
outras regiões espectrais.

**Como escolher `para_pcs` e `polyorder`?**

| Parâmetro | Valor mínimo | Valor típico | Valor máximo |
|-----------|-------------|-------------|-------------|
| `para_pcs` | 1 | 3 | 5 |
| `polyorder` | 0 | 2 | 4 |

Use `emsc.pca(parafina['r'], k=5)` para verificar quantas PCs são
necessárias para cobrir ≥ 99% da variância da parafina.
Polinômios de grau 2 (constante + linear + quadrático) são suficientes
para a maioria das variações de baseline em Micro-FTIR.

---

### `create_model_h2o(target, parafin, para_pcs, h2o, h2o_pcs, polyorder)` — Modelo EMSC com água

```python
model = emsc.create_model_h2o(
    target    = data_tecido,
    parafin   = data_parafina,
    para_pcs  = 3,
    h2o       = data_agua,
    h2o_pcs   = 2,
    polyorder = 2
)
```

**O que faz:**

Igual a `create_model()`, mas adiciona as PCs do vapor de água como
regressores extras:

```
Coluna 0:        Espectro de referência
Colunas 1..p+1:  Polinômios de baseline
Colunas seguintes: PCs da parafina
Últimas colunas: PCs da água (h2o_pcs colunas)
```

**Quando usar `create_model_h2o()` em vez de `create_model()`?**

| Situação | Recomendação |
|----------|-------------|
| Instrumento purga com N₂/ar seco e aquisição estável | `create_model()` suficiente |
| Bandas de vapor de água visíveis (~1350–2000 cm⁻¹ e ~3500–4000 cm⁻¹) | `create_model_h2o()` |
| Amostra em meio aquoso ou tecido hidratado | `create_model_h2o()` |
| Análise na região 1350–1800 cm⁻¹ onde a água interfere | `create_model_h2o()` |

**Parâmetros adicionais em relação a `create_model()`:**

| Parâmetro | Tipo | Descrição |
|-----------|------|-----------|
| `h2o` | dict | Espectros de referência de vapor de água |
| `h2o_pcs` | int | Número de PCs da água a incluir (tipicamente 1–3) |

**Região da água:**

As PCs da água são calculadas para `wn > 1350 cm⁻¹`, zerando fora dessa
janela. A região mais problemática do vapor de água no infravermelho médio
são as bandas de rotação-vibração entre 1350–2000 cm⁻¹.

---

### `emsc_model_view(model)` — Visualização dos componentes do modelo

```python
emsc.emsc_model_view(model)
```

**O que faz:**

Plota cada coluna da matriz do modelo em uma figura separada, com o
nome do regressor como título. Permite inspecionar visualmente se os
componentes têm a forma espectral esperada.

**Parâmetros:**

| Parâmetro | Tipo | Descrição |
|-----------|------|-----------|
| `model` | dict | Modelo retornado por `create_model()` ou `create_model_h2o()` |

**O que verificar em cada componente:**

| Componente | Forma esperada |
|-----------|----------------|
| `mean_spc` | Espectro típico do tecido com bandas de Amida I (1650), Amida II (1546), etc. |
| `poly 0` | Linha horizontal constante |
| `poly 1` | Linha diagonal (rampa) |
| `poly 2` | Parábola (forma de U ou ∩) |
| `para 1` | Picos em ~1378 e ~1468 cm⁻¹ (C–H deformação da parafina) |
| `h2o 1` | Oscilações características da água em 1350–2000 cm⁻¹ |

Se os componentes não tiverem a forma esperada, verifique se os dados de
parafina e água foram carregados corretamente e se os números de onda
são compatíveis com `target`.

---

### `emsc_fit(data, model)` — Aplicação da correção EMSC

```python
data = emsc.emsc_fit(data, model)
```

**O que faz:**

Aplica o modelo a todos os pixels de `data['r']`, gerando os espectros
corrigidos. É a função central do módulo.

**Parâmetros:**

| Parâmetro | Tipo | Descrição |
|-----------|------|-----------|
| `data` | dict | Dados hsp com `'r'` (modificado in-place) |
| `model` | dict | Modelo retornado por `create_model()` ou `create_model_h2o()` |

**Retorna:** o mesmo `data`, agora com:

| Chave adicionada | Conteúdo |
|-----------------|---------|
| `data['r']` | Espectros corrigidos (substitui os originais) |
| `data['EMSC_model']` | Matriz do modelo utilizada |
| `data['EMSC_coeff']` | Coeficientes β — shape `(n_pixels, n_regressores)` |

**O algoritmo passo a passo:**

```
Para cada pixel i:
    1. Resolve:  model['emsc_matrix'] · β_i = z_i
       via mínimos quadrados (np.linalg.lstsq)

    2. Subtrai os interferentes:
       z_corr_i = z_i - Σ(col_j · β_j)   para j = 1..n (tudo exceto col 0)

    3. Normaliza pela referência:
       z_final_i = z_corr_i / β_0
```

Na prática, o passo 1 é resolvido em lote para todos os pixels
simultaneamente (operação matricial), o que torna o cálculo eficiente
mesmo para imagens grandes (dezenas de milhares de pixels).

**Depois do `emsc_fit()` — verificação obrigatória:**

```python
# 1. Plota espectros aleatórios — baselines devem estar planas
import hsp_2026.sh as sh
sh.pplot(data, nspc=20)

# 2. Visualiza o mapa de a₀ — deve refletir a espessura da amostra
sh.emsc(data, b=0)

# 3. Histograma de a₀ para escolher limiar de qualidade
import hsp_2026.rg as rg
rg.emsc(data, a=0)

# 4. Controle de qualidade automático baseado em a₀
import hsp_2026.qt as qt
data = qt.otsu_emsc(data)
```

---

## Índices dos coeficientes em `data['EMSC_coeff']`

Após `emsc_fit()`, `data['EMSC_coeff']` tem shape `(n_pixels, n_regressores)`.
A coluna de cada coeficiente depende dos parâmetros usados na criação do modelo:

**Para `create_model(target, parafin, para_pcs=3, polyorder=2)`:**

| Coluna | Nome | Interpretação |
|--------|------|--------------|
| 0 | \(a_0\) | Escala da referência (espessura × conc.) |
| 1 | \(b_0\) | Offset de baseline (constante) |
| 2 | \(b_1\) | Inclinação linear de baseline |
| 3 | \(b_2\) | Curvatura quadrática de baseline |
| 4 | \(c_1\) | Contribuição da PC1 da parafina |
| 5 | \(c_2\) | Contribuição da PC2 da parafina |
| 6 | \(c_3\) | Contribuição da PC3 da parafina |

A fórmula geral do índice de início das PCs de parafina é
`1 + polyorder + 1` (1 referência + (polyorder+1) termos polinomiais).

---

## Exemplo completo

```python
import hsp_2026.file as file
import hsp_2026.prep as prep
import hsp_2026.emsc as emsc
import hsp_2026.rg   as rg
import hsp_2026.qt   as qt
import hsp_2026.sh   as sh

# ── 1. Carrega os dados ─────────────────────────────────────────────────
data     = file.fsm('tecido.fsm')
parafina = file.fsm('parafina_pura.fsm')

# ── 2. Corte espectral (antes do EMSC) ─────────────────────────────────
data     = prep.cut(data,     900, 3800)
parafina = prep.cut(parafina, 900, 3800)

# ── 3. Verifica quantas PCs cobrem a variância da parafina ─────────────
for k in [1, 2, 3, 5]:
    _, var = emsc.pca(parafina['r'], k=k)
    print(f'k={k}: {var:.1f}%')

# ── 4. Constrói o modelo ────────────────────────────────────────────────
model = emsc.create_model(
    target    = data,
    parafin   = parafina,
    para_pcs  = 3,
    polyorder = 2
)

# ── 5. Valida o modelo visualmente ─────────────────────────────────────
emsc.emsc_model_view(model)

# ── 6. Aplica a correção ────────────────────────────────────────────────
data = emsc.emsc_fit(data, model)

# ── 7. Verifica os espectros corrigidos ────────────────────────────────
sh.pplot(data, nspc=20)   # baselines devem estar planas e centradas

# ── 8. Controle de qualidade baseado no coeficiente a₀ ─────────────────
rg.emsc(data, a=0)
data = qt.otsu_emsc(data)

# ── 9. Análise multivariada ────────────────────────────────────────────
sh.pc(data, n=1)
```

---

## Comparativo `create_model` vs `create_model_h2o`

| Aspecto | `create_model` | `create_model_h2o` |
|---------|---------------|-------------------|
| Interferentes removidos | Parafina + baseline | Parafina + água + baseline |
| Dados extras necessários | Apenas parafina | Parafina + espectros de água |
| Número de regressores | 1 + (poly+1) + para_pcs | 1 + (poly+1) + para_pcs + h2o_pcs |
| Quando usar | Instrumento bem purgado | Bandas de água visíveis nos espectros |
| Complexidade | Menor | Maior |

Para a maioria dos experimentos de Micro-FTIR de tecidos em parafina
com boa purga do instrumento, `create_model()` com `para_pcs=3` e
`polyorder=2` é suficiente.
