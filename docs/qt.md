> ⚠️ **Documentação em construção.** Algumas informações podem estar imprecisas ou incompletas.

> ⚠️ **Documentação em construção.** Algumas informações podem estar imprecisas ou incompletas.

# Módulo `qt.py` — Controle de Qualidade Espectral

Este módulo é responsável por **remover pixels defeituosos** da imagem
hiperespectral antes das análises estatísticas. O nome "qt" vem de
*Quality Test* (teste de qualidade).

---

## Por que remover pixels?

Quando o microscópio de FTIR varre a amostra, nem todos os pixels
adquiridos contêm informação espectral válida. Vários problemas podem
ocorrer durante a preparação e a medição:

| Problema | Causa | Sintoma espectral |
|----------|-------|-------------------|
| Pixel de substrato | Região sem amostra (ex.: CaF₂ puro) | Absorção muito baixa em toda a faixa |
| Bolha de ar | Falha na preparação da lâmina | Espectro com forma irregular |
| Borda da amostra | Pixel parcialmente fora da amostra | Absorção muito baixa |
| Saturação | Amostra muito espessa | Bandas com absorção "cortada" no topo |
| Artefato de espalhamento | Tamanho de célula compatível com o comprimento de onda | Forma espectral ondulada |

Se esses pixels não forem removidos, eles **distorcem seriamente** os
algoritmos de K-Means, PCA e EMSC, pois representam uma classe espúria
(o "fundo") que vai competir com as classes biológicas reais.

---

## Como os filtros funcionam

Todos os filtros do `qt.py` seguem o mesmo princípio:

1. **Calcula uma métrica** para cada pixel (área, intensidade, coeficiente, etc.)
2. **Compara com um limiar**: pixels fora do intervalo `[a, b]` são marcados como ruins
3. **Atualiza a máscara `sel`**: os pixels removidos recebem `sel = False`
4. **Remove fisicamente de `data['r']`**: o array `r` é filtrado para conter apenas pixels válidos

```
Imagem original (100 × 100 = 10.000 pixels)
    │
    ├─ qt.otsu_area()  →  remove fundo e bordas  →  8.500 pixels válidos
    │
    └─ Análise multivariada com apenas pixels de amostra real
```

### A máscara `sel` — entendendo o mecanismo

A máscara `sel` é um vetor booleano que rastreia quais pixels do mapa
**original** ainda estão presentes. Isso é essencial para reconstruir
a imagem 2D depois da análise:

```python
# Pixels originais: dx × dy = 100 × 100 = 10.000
# Após filtros: data['r'].shape[0] == 8.500

# Para visualizar a imagem corretamente:
dplot = np.zeros(data['dx'] * data['dy'])   # grid com zeros no fundo
dplot[data['sel']] = valores_calculados      # preenche apenas os pixels válidos
dplot = dplot.reshape(data['dx'], data['dy']) # reconstrói a grade 2D
```

Todas as funções de visualização em `sh.py` fazem esse processo
automaticamente.

---

## Workflow recomendado: use `rg.py` para escolher os limiares

Antes de aplicar qualquer filtro com limiar manual (`a` e `b`), use o
módulo `rg.py` para visualizar a distribuição da métrica como histograma:

```python
import hsp_2026.rg as rg
import hsp_2026.qt as qt

# 1. Olha o histograma da área na região fingerprint
rg.area(data, 900, 1800)
# → histograma mostra dois grupos: fundo (~área < 50) e amostra (área > 50)

# 2. Aplica o filtro com os limiares identificados no histograma
data = qt.area(data, 900, 1800, a=50, b=500)
```

Para não ter que ajustar manualmente, use as funções `otsu_area()` ou
`otsu_emsc()`, que determinam o limiar automaticamente.

---

## Funções do módulo

### `area(data, ini, fim, a, b)` — Filtrar por área integrada

```python
data = qt.area(data, ini=900, fim=1800, a=50, b=500)
```

**O que faz:**

Calcula a **área integrada** da banda espectral entre `ini` e `fim` cm⁻¹
para cada pixel (usando a regra do trapézio), e remove os pixels cuja
área está fora do intervalo `[a, b]`.

**Por que a área é uma boa métrica?**

A área integrada é proporcional à concentração × espessura do material
que absorve naquela faixa. Na região fingerprint (900–1800 cm⁻¹):
- Áreas muito baixas → pixel de substrato ou borda fina (pouco material)
- Áreas muito altas → amostra muito espessa ou detector saturado

**Como usar `rg.area()` para escolher `a` e `b`:**

```python
rg.area(data, 900, 1800)
# O histograma mostrará (tipicamente) um pico estreito próximo de zero
# (fundo/bolhas) e uma distribuição mais larga (amostra).
# Escolha 'a' logo após o primeiro pico e 'b' na cauda direita.
```

**Saída no terminal:**

```
 teste de qualidade usando a área
 região 900 até 1800
 min_value: 50
 max_value: 500
 espectros removidos: 1523
```

---

### `intt(data, ini, a, b)` — Filtrar por intensidade num ponto

```python
data = qt.intt(data, ini=1650, a=0.05, b=1.5)
```

**O que faz:**

Extrai a **intensidade de absorbância no número de onda mais próximo de
`ini`** e remove pixels fora de `[a, b]`.

**Quando usar:**

Útil quando você tem uma **banda de referência confiável** que deve
estar presente em todos os pixels válidos. Por exemplo:

- `ini=1650` cm⁻¹ (Amida I): remove pixels sem proteína — adequado para
  filtrar fundo em imagens de tecido
- `ini=1740` cm⁻¹ (C=O lipídico): remove pixels sem lipídio — adequado
  para estudos de lipossomas ou membranas

**Diferença em relação a `area()`:**

`intt()` usa apenas **um ponto** do espectro, o que o torna mais sensível
a ruído pontual. `area()` usa a integral de uma região inteira, sendo mais
robusto. Prefira `area()` quando possível; use `intt()` quando tiver uma
banda muito característica e específica.

**⚠️ Atenção:**

A função busca o número de onda exatamente igual a `ini` em `data['wn']`.
Como os números de onda são pontos flutuantes, é possível que nenhum ponto
corresponda exatamente. Use `rg.intt(data, ini)` antes para confirmar que
o ponto existe, ou prefira `area()` com uma janela estreita em torno do pico.

---

### `emsc(data, ini, a, b)` — Filtrar por coeficiente EMSC

```python
# Filtra pelo coeficiente a₀ (índice 0): escala da referência
data = qt.emsc(data, ini=0, a=0.3, b=2.0)
```

**O que faz:**

Remove pixels cujo **coeficiente EMSC de índice `ini`** está fora de
`[a, b]`. Requer que `emsc.emsc_fit()` tenha sido executado antes,
pois o campo `data['EMSC_coeff']` precisa existir.

**Entendendo os índices dos coeficientes:**

O modelo EMSC decompõe cada espectro como:

$$
z = a_0 \cdot \mathbf{r}_{\text{ref}} + b_0 + b_1\lambda + b_2\lambda^2 + \ldots + c_k\mathbf{p}_k + \varepsilon
$$

| Índice | Coeficiente | Significado |
|--------|-------------|-------------|
| 0 | \(a_0\) | Escala da referência — proporcional à quantidade de material biológico |
| 1 | \(b_0\) | Offset da baseline polinomial |
| 2 | \(b_1\) | Inclinação linear da baseline |
| 3 | \(b_2\) | Curvatura da baseline |
| 4+ | \(c_k\) | Contribuição da parafina, água, etc. |

**O caso mais comum — filtrar pelo \(a_0\):**

```python
# Após aplicar EMSC:
data = emsc.emsc_fit(data, model)

# Visualiza o histograma do coeficiente a0
rg.emsc(data, a=0)

# Aplica o filtro
data = qt.emsc(data, ini=0, a=0.2, b=3.0)
```

O \(a_0\) é a melhor métrica de qualidade depois do EMSC porque
representa diretamente a quantidade de sinal biológico no pixel,
normalizada pelo espectro de referência e corrigida para parafina e água.

---

### `mean(data, ini, fim, a, b)` — Filtrar por escala relativa ao espectro médio

```python
data = qt.mean(data, ini=900, fim=1800, a=0.5, b=2.0)
```

**O que faz:**

Para cada pixel, ajusta a equação:
$$
\mathbf{z} \approx \alpha \cdot \bar{\mathbf{z}} + c
$$

via mínimos quadrados, onde \(\bar{\mathbf{z}}\) é o espectro médio
global da região. O coeficiente \(\alpha\) mede o quanto aquele pixel
precisa ser **escalado** para se aproximar da média da imagem.

Remove pixels cujo \(\alpha\) está fora de `[a, b]`.

**Interpretação do coeficiente α:**

| Valor de α | Interpretação |
|-----------|---------------|
| α ≈ 1 | Pixel "típico" — similar ao padrão médio da imagem |
| α < 0.5 | Pixel muito fraco: borda, bolha, pouco material |
| α > 2.0 | Pixel muito intenso: espessura anormal, artefato |
| α < 0 | Espectro com forma invertida: artefato grave |

**Vantagem em relação a `area()`:**

`mean()` é menos sensível ao tipo de tecido e mais robusto a variações
de composição, porque usa o **espectro inteiro** (forma + amplitude)
como referência. É especialmente útil em amostras heterogêneas onde a
área integrada pode variar muito por razões biológicas reais.

**Exemplo de uso com `rg.py`:**

```python
# Não existe rg.mean() — use sh.mean() para visualizar o mapa
# e depois defina os limiares visualmente

# Inspete o mapa de α
import hsp_2026.sh as sh
sh.mean(data, 900, 1800)

# Depois aplica o filtro
data = qt.mean(data, 900, 1800, a=0.3, b=1.8)
```

---

### `otsu_area(data, ini1, fim1, k=1)` — Limiar automático de Otsu sobre a área

```python
# Sem ajuste manual — limiar determinado automaticamente
data = qt.otsu_area(data, ini1=900, fim1=1800)

# Com fator k para afinar o limiar
data = qt.otsu_area(data, ini1=900, fim1=1800, k=0.8)
```

**O que faz:**

Aplica o **método de Otsu** para determinar automaticamente o limiar
de corte da área espectral, sem precisar inspecionar histogramas ou
definir `a` e `b` manualmente.

**O que é o método de Otsu?**

Otsu é um algoritmo clássico de visão computacional (originalmente para
binarização de imagens em preto e branco) que encontra o limiar \(t\)
que **maximiza a separância entre dois grupos**:

$$
t^* = \operatorname*{arg\,max}_t \left[ w_0(t)\,\sigma_0^2(t) + w_1(t)\,\sigma_1^2(t) \right]
$$

onde \(w_0, w_1\) são as frações de pixels em cada grupo e
\(\sigma_0^2, \sigma_1^2\) são suas variâncias.

No contexto de Micro-FTIR, os dois grupos naturais são:
- **Grupo 0** (fundo): pixels de substrato, bolhas e bordas → áreas baixas
- **Grupo 1** (amostra): pixels com tecido ou células → áreas altas

O Otsu encontra o limiar que melhor separa esses dois grupos de forma
totalmente automática.

**O parâmetro `k`:**

Multiplica o limiar calculado pelo Otsu:
- `k = 1.0` → usa exatamente o limiar de Otsu (padrão)
- `k = 0.8` → limiar 20% mais baixo → retém mais pixels (mais permissivo)
- `k = 1.2` → limiar 20% mais alto → remove mais pixels (mais restritivo)

Use `k < 1` quando o Otsu for muito agressivo (remover pixels de amostra
real); use `k > 1` quando ainda houver pixels de fundo passando.

**Esta é a função de controle de qualidade mais recomendada para iniciantes**
por não exigir nenhum parâmetro manual.

**Exemplo:**

```python
import hsp_2026.file as file
import hsp_2026.qt   as qt
import hsp_2026.sh   as sh

data = file.fsm('tecido.fsm')

# Antes: visualiza a imagem com todos os pixels (incluindo fundo)
sh.intt(data, 1650)

# Aplica controle de qualidade automático
data = qt.otsu_area(data, ini1=900, fim1=1800)

# Depois: apenas pixels de amostra
sh.intt(data, 1650)
```

---

### `otsu_emsc(data)` — Limiar automático de Otsu sobre o coeficiente EMSC \(a_0\)

```python
# Requer que emsc_fit() já tenha sido aplicado
data = emsc.emsc_fit(data, model)
data = qt.otsu_emsc(data)
```

**O que faz:**

Aplica o método de Otsu automaticamente sobre o **coeficiente \(a_0\)** do
modelo EMSC (índice 0 em `data['EMSC_coeff']`), que representa a escala
da referência após a correção de espalhamento.

**Por que usar \(a_0\) em vez da área?**

Depois da aplicação do EMSC, o coeficiente \(a_0\) é uma métrica mais
precisa de qualidade do que a área integrada bruta porque:

1. Já está corrigido para variações de parafina e água (nas amostras embebidas)
2. Corresponde diretamente à "quantidade de amostra biológica" no pixel,
   sem interferências físicas
3. A separação entre fundo (a₀ ≈ 0) e amostra (a₀ > 0) tende a ser
   mais nítida do que na área bruta, facilitando o Otsu

**Fluxo completo com EMSC:**

```python
data  = file.fsm('tecido.fsm')
data  = prep.cut(data, 900, 1800)
data  = prep.snv(data)

# Constrói e aplica o modelo EMSC
model = emsc.create_model(referencia, parafina_loadings, 3, 2)
data  = emsc.emsc_fit(data, model)

# Controle de qualidade baseado no coeficiente EMSC a₀
data  = qt.otsu_emsc(data)

# Agora os dados estão limpos e corrigidos para análise
sh.pc(data, n=1)
```

---

## Comparativo das funções

| Função | Precisa de EMSC? | Parâmetros manuais? | Métrica usada |
|--------|-----------------|---------------------|--------------|
| `area()` | Não | Sim (a, b) | Área integrada (trapézio) |
| `intt()` | Não | Sim (a, b) | Intensidade num ponto |
| `mean()` | Não | Sim (a, b) | Coeficiente α (regressão) |
| `emsc()` | **Sim** | Sim (a, b) | Coeficiente EMSC |
| `otsu_area()` | Não | Não (k opcional) | Área + Otsu automático |
| `otsu_emsc()` | **Sim** | Não | Coeficiente \(a_0\) + Otsu |

---

## Ordem recomendada no pipeline

O controle de qualidade deve ser feito **antes** do pré-processamento
principal, para não desperdiçar tempo processando pixels ruins:

```python
data = file.fsm('tecido.fsm')

# ── 1. Controle de qualidade ─────────────────────────────────────────────
data = qt.otsu_area(data, ini1=900, fim1=1800)
# Remove fundo, bolhas e bordas automaticamente

# ── 2. Pré-processamento ─────────────────────────────────────────────────
data = prep.cut(data, 900, 1800)
data = prep.snv(data)

# ── 3. Análise ───────────────────────────────────────────────────────────
sh.pc(data, n=1)
km.fit(data, k=5)
```

Se você usar EMSC, o `otsu_emsc()` pode ser aplicado depois do EMSC
como um segundo estágio de limpeza mais refinado.
