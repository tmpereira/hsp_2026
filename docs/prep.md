> ⚠️ **Documentação em construção.** Algumas informações podem estar imprecisas ou incompletas.

> ⚠️ **Documentação em construção.** Algumas informações podem estar imprecisas ou incompletas.

# Módulo `prep.py` — Pré-Processamento Espectral

Este módulo contém todas as ferramentas para **preparar** os espectros
antes de qualquer análise química ou estatística. Pular o pré-processamento
é um erro comum entre iniciantes: os dados brutos de um microscópio de FTIR
estão cheios de variações que não têm nada a ver com a composição química
da amostra — diferenças de espessura, espalhamento de luz, ruído eletrônico,
deriva de linha de base. O objetivo do pré-processamento é remover essas
interferências e deixar apenas o sinal químico.

---

## Por que pré-processar?

Imagine dois pixels de uma mesma célula: um está em cima de uma região
ligeiramente mais espessa da seção histológica. Esse pixel vai apresentar
absorbâncias maiores em todas as bandas — não porque tem mais proteína ou
lipídio, mas simplesmente porque o feixe de luz atravessou mais material.
Se você aplicar K-Means ou PCA sem corrigir isso, o algoritmo vai separar
pixels por espessura, não por composição química. Isso é exatamente o que
queremos evitar.

O pipeline típico de pré-processamento é:

```
Dados brutos
    │
    ├─ 1. cut()    — remove regiões sem informação
    ├─ 2. offset() — corrige deslocamento vertical
    ├─ 3. golay()  — suaviza ou calcula derivada
    └─ 4. snv()    — normaliza cada espectro individualmente
         │
         └─► Dados prontos para análise
```

Cada etapa é independente e pode ser combinada de formas diferentes
dependendo do tipo de amostra e dos objetivos científicos.

---

## Rastreamento automático com `log`

Todas as funções de `prep.py` registram automaticamente o que foi feito no
campo `data['log']`. Isso garante rastreabilidade: você sempre saberá exatamente
quais etapas foram aplicadas em um conjunto de dados.

```python
data = file.fsm('amostra.fsm')
data = prep.cut(data, 900, 1800)
data = prep.snv(data)

print(data['log'])
# abrindo arquivo FSM: amostra.fsm
# restrição espectral de 900 cm-1 até 1800 cm-1
# normalização vetorial(SNV) em única região
```

---

## Funções do módulo

### `cut(data, a, b)` — Recorte espectral

```python
data = prep.cut(data, 900, 1800)
```

**O que faz:**

Descarta todos os pontos espectrais fora do intervalo `[a, b]` cm⁻¹.
Após essa operação, `data['r']` terá menos colunas e `data['wn']` será
um vetor menor.

**Por que usar?**

Um espectro de infravermelho médio cobre tipicamente 400–4000 cm⁻¹, mas
a região biológica mais informativa é a **região fingerprint** (900–1800 cm⁻¹),
que contém as bandas de proteínas (Amida I a 1650 cm⁻¹, Amida II a 1550 cm⁻¹),
lipídios (1740 cm⁻¹), ácidos nucleicos (1085 cm⁻¹) e carboidratos.
Trabalhar apenas nessa região:

- Reduz o custo computacional (menos pontos por espectro)
- Evita que regiões com absorção forte de água atmosférica (em torno de
  3400 cm⁻¹ e 1640 cm⁻¹) distorçam os algoritmos de normalização

**Como funciona internamente:**

```python
sel1 = data['wn'] > a    # True acima do limite inferior
sel2 = data['wn'] < b    # True abaixo do limite superior
# Converte para int, soma e subtrai 1: resultado é 1 (dentro) ou 0/negativo (fora)
ver = (sel1.astype(int) + sel2.astype(int)) - 1
sel = ver.astype(bool)   # True apenas onde ver == 1

data['r']  = data['r'][:, sel]
data['wn'] = data['wn'][sel]
```

Essa implementação usa uma soma aritmética como substituto do operador `AND`
lógico: somente quando `sel1` e `sel2` são simultaneamente `True` (ambas valem 1),
a soma menos 1 resulta em 1 (True).

**Exemplo prático:**

```python
# Antes do corte
print(data['wn'][[0, -1]])   # [399.9, 3999.1]
print(data['r'].shape)       # (3600, 1800)

data = prep.cut(data, 900, 1800)

# Depois do corte
print(data['wn'][[0, -1]])   # [900.4, 1799.8]
print(data['r'].shape)       # (3600, 426)  ← muito menos colunas
```

---

### `snv(data)` / `norm(data)` — Normalização vetorial (SNV)

```python
data = prep.snv(data)
# ou, equivalentemente:
data = prep.norm(data)
```

**O que faz:**

Aplica a **SNV** (*Standard Normal Variate*) a cada espectro individualmente.
Para cada espectro \(\mathbf{z}\) com \(p\) pontos:

$$
z_i^{\text{SNV}} = \frac{z_i - \bar{z}}{\sigma_z}, \quad i = 1, \ldots, p
$$

onde \(\bar{z}\) é a média e \(\sigma_z\) é o desvio padrão calculados
**ao longo dos pontos espectrais daquele único espectro** (não entre espectros!).

Após a SNV, cada espectro tem média 0 e desvio padrão 1.

**Por que usar?**

A SNV é a técnica de normalização mais usada em espectroscopia vibracional
porque remove simultaneamente:

| Tipo de variação | Causa | Efeito na SNV |
|-----------------|-------|---------------|
| Variação multiplicativa | Diferença de espessura | Removida pela divisão por σ |
| Variação de offset | Espalhamento de luz | Removida pela subtração da média |
| Diferença de concentração total | Quantidade de material no pixel | Removida |

Depois da SNV, dois pixels com a mesma composição química relativa terão
espectros idênticos, independentemente da espessura ou da quantidade
de material.

**⚠️ Limitação importante:**

A SNV remove informação sobre a quantidade absoluta de material.
Se o objetivo é quantificar concentrações absolutas, a SNV não é adequada.
Para isso, use `emsc.emsc_fit()`, que modela e preserva a escala relativa
através do coeficiente \(a_0\).

**Como funciona internamente:**

```python
media = np.mean(spc, axis=1)   # média de cada espectro → shape (n_espectros,)
std   = np.std(spc, axis=1)    # desvio de cada espectro → shape (n_espectros,)

# [:, None] transforma (n,) em (n, 1) para funcionar com broadcasting do NumPy
data['r'] = (spc - media[:, None]) / std[:, None]
```

O truque do `[:, None]` (equivalente a `reshape(-1, 1)`) é um padrão
fundamental do NumPy: permite subtrair ou dividir um vetor coluna de
uma matriz sem escrever loops.

---

### `golay(data, diff, order, win)` — Filtro de Savitzky-Golay

```python
# Suavização simples
data = prep.golay(data, diff=0, order=2, win=9)

# Primeira derivada
data = prep.golay(data, diff=1, order=2, win=11)

# Segunda derivada
data = prep.golay(data, diff=2, order=2, win=11)
```

**O que faz:**

O filtro de **Savitzky-Golay** é um método de suavização que ajusta um
polinômio local a cada janela de `win` pontos consecutivos do espectro.
O valor filtrado no centro da janela é o valor do polinômio naquele ponto.
Se `diff > 0`, em vez do valor do polinômio, calcula-se sua derivada de
ordem `diff`.

**Por que usar derivadas?**

| `diff` | Operação | Remove | Realça |
|--------|----------|--------|--------|
| 0 | Suavização | Ruído de alta frequência | Nada |
| 1 | 1ª derivada | Offset constante | Flancos de bandas |
| 2 | 2ª derivada | Baseline linear | Posição de pico, bandas sobrepostas |

A segunda derivada é especialmente poderosa em espectroscopia porque
converte picos de absorbância em **vales negativos** e pode resolver
bandas sobrepostas que parecem um único pico largo no espectro original.

**Como escolher os parâmetros?**

- `order` (grau do polinômio): 2 ou 3 são os mais comuns. Valores maiores
  podem distorcer picos estreitos.
- `win` (tamanho da janela): deve ser ímpar e maior que `order`.
  Para resolução espectral de 4 cm⁻¹, janelas de 9–15 pontos funcionam bem.
  Janelas maiores suavizam mais, mas podem alargar ou distorcer bandas.

**Como funciona internamente:**

Em vez de aplicar a convolução ponto a ponto em um loop (lento), a
implementação usa uma **matriz esparsa de convolução** `D` e aplica
tudo de uma vez com um produto matricial:

```
r_filtrado = r · D        (produto matricial: n_pixels × n_pontos)
```

A matriz `D` é construída colocando os coeficientes do filtro SG em
diagonais deslocadas. Cada linha `i` de `D` define como o ponto `i` do
espectro filtrado depende dos pontos vizinhos do espectro original.

> **Nota sobre bordas:** as primeiras e últimas `n = (win-1)/2` colunas
> de `D` são zeradas para evitar artefatos nas bordas do espectro, onde
> a janela não cabe completamente.

**Exemplo visual (segunda derivada):**

```
Espectro original:          │  Após 2ª derivada:
                            │
     ╭─────╮                │      ↓ vales indicam picos
─────╯     ╰─────           │  ────╮     ╭────
    pico largo              │      ╰─────╯
                            │   dois vales = dois picos resolvidos
```

---

### `norm2r(data, ini1, fim1, ini2, fim2)` — SNV em duas regiões

```python
# Normaliza e concatena região fingerprint + região lipídica
data = prep.norm2r(data, 900, 1800, 2800, 3050)
```

**O que faz:**

Versão especializada da SNV para quando se quer trabalhar com **duas
regiões espectrais separadas** ao mesmo tempo.

Cada região é normalizada independentemente pela sua própria média e
desvio padrão, e depois as duas regiões são concatenadas em um único
espectro. O vetor `wn` é igualmente concatenado.

**Por que não usar `cut()` + `snv()` para cada região?**

A SNV sobre uma única região grande pode ser dominada pela região com
intensidade maior, distorcendo a região menor. Ao normalizar cada região
separadamente, garante-se que ambas contribuam igualmente para as análises
multivariadas subsequentes.

**Caso de uso típico:**

Região fingerprint (900–1800 cm⁻¹) + região de estiramento C–H dos
lipídios (2800–3050 cm⁻¹). Essas duas regiões têm intensidades muito
diferentes e fornecer informações complementares sobre a composição
molecular da célula.

```python
data = prep.norm2r(data, 900, 1800, 2800, 3050)
print(data['r'].shape)   # (n_pixels, n_pontos_r1 + n_pontos_r2)
```

---

### `pcares(data, n)` — Denoising por PCA

```python
data = prep.pcares(data, n=10)
```

**O que faz:**

Remove ruído reconstruindo os espectros usando apenas as `n` primeiras
componentes principais (PCs).

A ideia vem da teoria da PCA: os primeiros PCs capturam as direções de
maior variância nos dados (o "sinal"), enquanto os PCs tardios capturam
variações aleatórias de baixa variância (o "ruído" do detector).

**Algoritmo passo a passo:**

```
1. Calcula o espectro médio global: μ = média(r, eixo=pixels)

2. Centraliza os dados: X_c = r - μ

3. PCA decompõe X_c em scores (T) e loadings (P):
   X_c ≈ T · P^T

4. Zera os scores das PCs além da n-ésima:
   T[:, n:] = 0

5. Reconstrói os espectros filtrados:
   r_denoised = μ + T_filtrado · P^T
```

**Como escolher `n`?**

- Valores muito baixos (ex.: `n=3`) → muita suavização, pode remover
  variações químicas reais
- Valores altos (ex.: `n=50`) → pouco benefício, ruído não é removido
- Valores típicos em Micro-FTIR biológico: 5 a 20

Uma forma prática de escolher é plotar o percentual de variância acumulada
pelos primeiros PCs e escolher `n` onde a curva "dobra" (critério do cotovelo).

---

### `napc(dados, noise, npcs)` — Denoising por NAPC

```python
# Requer espectros de ruído puro (ex.: substrato sem amostra)
data  = file.fsm('amostra.fsm')
ruido = file.fsm('substrato_vazio.fsm')

data = prep.napc(data, ruido, npcs=10)
```

**O que faz:**

O **NAPC** (*Noise-Adjusted Principal Components*) é uma versão avançada
da PCA que leva em conta explicitamente o **perfil do ruído instrumental**.

Em vez de maximizar a variância total (como a PCA padrão), o NAPC
maximiza a **razão sinal-ruído** em cada componente. O resultado são
componentes espectrais cujas direções são mais relevantes quimicamente,
porque o ruído instrumental é "subtraído" matematicamente.

**Por que é melhor que `pcares()` em alguns casos?**

O ruído de um microscópio de FTIR não é uniforme em todos os números de
onda — ele tem um perfil espectral próprio que depende da fonte de luz,
do detector e dos espelhos. O `pcares()` trata todos os pontos iguais,
mas o NAPC usa espectros de ruído medidos experimentalmente para levar
em conta essa não-uniformidade.

**Pré-requisito:** você precisa de espectros de "ruído puro",
medidos no mesmo instrumento, nas mesmas condições, mas sem a amostra
(ex.: espectros do substrato de CaF₂ limpo, ou diferença entre duas
aquisições consecutivas da mesma área).

**Algoritmo (simplificado):**

$$
\Sigma_d = \text{cov}(\mathbf{X}_{\text{dados}}) \qquad
\Sigma_n = \text{cov}(\mathbf{X}_{\text{ruído}})
$$

Encontra as direções que maximizam \(\Sigma_d\) com a restrição de que
o ruído ao longo dessas direções seja unitário (matriz de branqueamento
\(\mathbf{F}\) derivada de \(\Sigma_n\)).

---

### `offset(data, ini, fim)` — Remoção de offset

```python
data = prep.offset(data, ini=1800, fim=2000)
```

**O que faz:**

Para cada espectro, encontra o **valor mínimo** na região `[ini, fim]` cm⁻¹
e subtrai esse valor de todo o espectro. Depois da operação, todos os
espectros terão um mínimo próximo de zero naquela região.

**Por que usar?**

É comum que espectros de FTIR apresentem um deslocamento vertical
(offset) causado por:
- Espalhamento de Mie (interferência entre o tamanho das células e o
  comprimento de onda do IR)
- Diferenças no nivelamento do substrato
- Deriva eletrônica do instrumento

A região `[1800, 2000]` cm⁻¹ é frequentemente escolhida como referência
porque não contém bandas de absorção em tecidos biológicos —
qualquer sinal nessa região é basicamente offset ou ruído.

**Como escolher a região de referência:**

Use `sh.intt(data, 1900)` para verificar visualmente se a intensidade
nessa região varia sistematicamente entre pixels. Se variar muito, o
offset está presente.

**Visualização do efeito:**

```
Antes:                     Depois:
    ╭─╮                       ╭─╮
    │ │         ──────         │ │
────╯ ╰──  →            →  ───╯ ╰──
               offset           zero
   offset elevado              removido
```

---

### `binned(data)` — Binagem 2×2 pixels

```python
data = prep.binned(data)
```

**O que faz:**

Agrupa blocos de 2×2 pixels em um único pixel, calculando a média dos
4 espectros. A imagem resultante tem a metade das linhas e a metade das
colunas, mas com melhor razão sinal-ruído.

**Por que funciona?**

Se cada pixel tem ruído aleatório com desvio padrão σ, a média de 4
pixels independentes tem desvio padrão σ/√4 = σ/2. Ou seja, a binagem
2×2 **dobra o SNR** em amplitude (ou quadruplica em potência).

Isso é especialmente útil em regiões de absorção fraca (como os
lipídios na região 2800–3050 cm⁻¹) onde o sinal individual de cada
pixel pode estar no limite do ruído do detector.

**Troca-off:**

Ganho de SNR ×2 em amplitude, mas perde-se resolução espacial: cada
pixel da imagem binada representa uma área 4× maior da amostra.
Para células grandes (> 20 µm) e imagens com alta resolução (pixel de 5 µm),
a binagem é um excelente investimento.

**Dimensões após binagem:**

```
Antes: dx=100, dy=80  →  100×80 = 8000 pixels
Depois: dx=49, dy=39  →  49×39  = 1911 pixels  (≈ 1/4 dos pontos)
```

> **Atenção:** a binagem redefine `dx`, `dy` e `sel`. Use-a antes de
> qualquer função de visualização de imagem (`sh.intt`, `sh.area`, etc.).

---

### `rand(data, k)` — Seleção aleatória de pixels

```python
# Seleciona 500 pixels aleatórios para teste rápido
data_amostra = prep.rand(data, k=500)
```

**O que faz:**

Seleciona `k` pixels aleatórios da imagem e descarta todos os outros.
O campo `sel` é atualizado para refletir apenas os pixels mantidos.

**Casos de uso:**

- **Testes de parâmetros:** antes de rodar um K-Means em 10.000 pixels,
  teste com 500 para verificar se os parâmetros fazem sentido
- **Equalização de classes:** se uma imagem tem muito mais pixels de
  uma classe do que outra, selecione aleatoriamente o mesmo número de
  pixels de cada para balancear os dados de treinamento
- **Redução para EMSC:** a construção do modelo EMSC pode usar um
  subconjunto representativo em vez de todos os pixels

---

### `dsample(data)` — Subamostragem espacial

```python
data = prep.dsample(data)
```

**O que faz:**

Reduz a resolução espacial mantendo apenas pixels em posições alternadas
da grade (passo de 2 em linhas e colunas). Diferente de `binned()`,
**não calcula média** — simplesmente descarta pixels intermediários.

**Quando usar:**

- Quando a resolução espacial original é maior do que a resolução de difração
  do instrumento (o que é comum em imagens de FTIR com pixel de 2–3 µm
  numa fonte de IR clássica com resolução lateral de ~10 µm)
- Para reduzir rapidamente o número de pontos sem o custo computacional
  de calcular médias

**Diferença em relação a `binned()`:**

| | `binned()` | `dsample()` |
|-|------------|-------------|
| Método | Média de 4 pixels | Descarte de pixels alternados |
| Melhora SNR? | Sim (×2) | Não |
| Preserva espectros? | Não (cria espectros médios) | Sim (mantém espectros originais) |
| Velocidade | Mais lento | Mais rápido |

---

## Pipelines de pré-processamento comuns

### Tecido biológico em parafina (seção histológica)

```python
data = file.fsm('tecido_higido.fsm')
data = prep.cut(data, 900, 1800)        # foco na região fingerprint
data = prep.offset(data, 1800, 2000)    # remove offset vertical
data = prep.golay(data, diff=2, order=2, win=11)  # segunda derivada
data = prep.snv(data)                   # normaliza forma espectral
```

### Células em suspensão ou cultivo

```python
data = file.fsm('celulas.fsm')
data = prep.cut(data, 900, 1800)
data = prep.snv(data)
data = prep.pcares(data, n=15)          # denoising por PCA
```

### Análise de lipídios + proteínas simultaneamente

```python
data = file.fsm('tecido.fsm')
# Normaliza e concatena as duas regiões de interesse
data = prep.norm2r(data, 1500, 1800, 2800, 3050)
```

### Imagem grande com poucos recursos computacionais

```python
data = file.fsm('imagem_grande.fsm')
data = prep.cut(data, 900, 1800)
data = prep.binned(data)               # reduz para 1/4 dos pixels
data = prep.snv(data)
```
