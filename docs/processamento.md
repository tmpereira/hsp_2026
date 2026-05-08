> ⚠️ **Documentação em construção.** Algumas informações podem estar imprecisas ou incompletas.

> ⚠️ **Documentação em construção.** Algumas informações podem estar imprecisas ou incompletas.

# Pré-Processamento Espectral

O cubo hiperespectral bruto de Micro-FTIR contém, além do sinal químico de
interesse, artefatos de origem instrumental e física que precisam ser
removidos antes de qualquer análise quantitativa ou multivariada.

Esta página descreve as principais etapas implementadas no módulo `prep` e
no módulo `emsc`.

---

## 1. Recorte Espectral (`prep.cut`)

A região de infravermelho médio útil para tecidos biológicos é geralmente
a **região fingerprint** (900 – 1800 cm⁻¹) e a região de **estiramento C–H**
(2800 – 3050 cm⁻¹).  Remover regiões sem informação reduz o custo
computacional e evita que ruído influencie os modelos.

```python
# Mantém apenas a região fingerprint
data = prep.cut(data, 900, 1800)
```

---

## 2. Normalização por Variação Normal Padrão — SNV

A **SNV** (*Standard Normal Variate*) remove variações multiplicativas e
aditivas causadas por diferenças de espessura e espalhamento de luz entre
pixels.  Para cada espectro \(\mathbf{z}\):

$$
z_i^{\text{SNV}} = \frac{z_i - \bar{z}}{\sigma_z}
$$

onde \(\bar{z}\) é a média e \(\sigma_z\) o desvio padrão do espectro.
Após a SNV, todos os espectros têm média 0 e variância 1, tornando-os
comparáveis independentemente da quantidade de material no caminho óptico.

```python
data = prep.snv(data)
```

!!! warning "SNV vs. normalização por área de banda"
    A SNV usa todo o espectro como referência interna e é preferível quando
    não existe uma banda inerte de referência.  Se a sua amostra contém uma
    banda de intensidade constante (ex.: banda de parafina a 2850 cm⁻¹ em
    seções histológicas embebidas em parafina), considere a normalização por
    área de banda com `prep.norm2r`.

---

## 3. Filtro de Savitzky-Golay e Derivadas

O filtro de Savitzky-Golay ajusta um polinômio local de grau `order` a uma
janela de `win` pontos e calcula a derivada de ordem `diff`.  É amplamente
utilizado em espectroscopia para:

- **diff=0**: suavização — reduz ruído de alta frequência sem distorcer picos
- **diff=1**: primeira derivada — remove baseline linear e acentua picos
- **diff=2**: segunda derivada — remove baseline quadrático e resolve bandas sobrepostas

A implementação usa matrizes esparsas diagonais para eficiência:

$$
\hat{z} = \mathbf{D}^{(\text{diff})} \cdot \mathbf{z}
$$

```python
# Suavização: polinômio de 2º grau, janela de 9 pontos
data = prep.golay(data, diff=0, order=2, win=9)

# Segunda derivada: polinômio de 2º grau, janela de 11 pontos
data = prep.golay(data, diff=2, order=2, win=11)
```

!!! tip "Escolha do tamanho da janela"
    Janelas maiores produzem mais suavização mas podem distorcer bandas
    estreitas.  Para a região fingerprint (resolução típica de 4 cm⁻¹),
    janelas de 9–15 pontos são comuns na literatura.

---

## 4. Correção de Espalhamento Multiplicativo Estendido — EMSC

A **EMSC** (*Extended Multiplicative Scatter Correction*) modela cada
espectro como uma combinação linear de um espectro de referência químico
puro \(\mathbf{r}\), de um baseline polinomial e de espectros interferentes
(ex.: parafina, água atmosférica):

$$
\mathbf{z} = a_0 \cdot \mathbf{r} + \sum_{j=0}^{p} b_j \lambda^j
           + \sum_k c_k \mathbf{p}_k^{\text{par}}
           + \sum_l d_l \mathbf{p}_l^{\text{H_2O}} + \boldsymbol{\varepsilon}
$$

Os coeficientes são estimados por mínimos quadrados e o espectro corrigido é:

$$
\mathbf{z}^{\text{EMSC}} = \frac{\mathbf{z} - \text{baseline} - \text{interferentes}}{a_0}
$$

O coeficiente \(a_0\) é proporcional à espessura/concentração efetiva do
material de interesse.  O mapa de \(a_0\) é uma das imagens mais informativas
para avaliar a homogeneidade da seção histológica.

```python
import hsp_2026.emsc as emsc

# Calcula os PCs da parafina (necessário antes de criar o modelo)
para_loadings, para_var = emsc.pca(parafina_data, k=5)

# Cria o modelo EMSC com baseline de grau 2 e 3 PCs de parafina
model = emsc.create_model(
    target    = espectro_referencia,
    parafin   = para_loadings,
    para_pcs  = 3,
    polyorder = 2
)

# Aplica a correção a todos os espectros da imagem
data = emsc.emsc_fit(data, model)

# Visualiza o coeficiente a0 como imagem
import hsp_2026.sh as sh
sh.emsc(data, b=0)   # b=0 → coeficiente a0
```

---

## 5. Remoção de Offset de Baseline (`prep.offset`)

Remove um deslocamento vertical aditivo usando o mínimo da absorbância numa
região de referência (idealmente uma região sem bandas de absorção):

$$
z_i^{\text{corr}} = z_i - \min_{j \in [a, b]} z_j
$$

```python
# Remove offset usando o mínimo entre 1800–2000 cm⁻¹ (região sem bandas)
data = prep.offset(data, ini=1800, fim=2000)
```

---

## 6. Análise de Componentes Principais para Denoising (`prep.pcares`)

Reconstrói os espectros usando apenas os primeiros `n` componentes
principais, descartando os componentes tardios que capturam majoritariamente
ruído:

$$
\mathbf{X} \approx \mathbf{T}_n \cdot \mathbf{P}_n^T
$$

```python
# Reconstrói com as 10 primeiras PCs (descarta ruído nas PCs > 10)
data = prep.pcares(data, n=10)
```

---

## Fluxo de pré-processamento típico para tecido biológico

```python
import hsp_2026.file as file
import hsp_2026.prep as prep
import hsp_2026.qt   as qt

# 1. Carrega o arquivo
data = file.fsm('tecido.fsm')

# 2. Controle de qualidade: remove pixels com baixa absorção (fundo)
data = qt.otsu_area(data, ini1=900, fim1=1800)

# 3. Recorte espectral
data = prep.cut(data, 900, 1800)

# 4. Suavização
data = prep.golay(data, diff=0, order=2, win=9)

# 5. Normalização SNV
data = prep.snv(data)
```
