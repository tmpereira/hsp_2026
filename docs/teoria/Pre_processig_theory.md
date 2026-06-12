> ⚠️ **Documentação em construção.** Algumas informações podem estar imprecisas ou incompletas.

# Pré-processamento e Análise de Dados em Microespectroscopia FT-IR

> **Para quem é este texto?** Este material foi escrito para estudantes de graduação em Engenharia e Ciências Exatas que estão tendo seu primeiro contato com análise de dados espectrais. O objetivo é construir a intuição passo a passo, do dado bruto até a imagem química final.

---

## 1. O que é um dado de microespectroscopia FT-IR?

Antes de falar em pré-processamento, precisamos entender  **como os dados são organizados** .

Em um experimento de microespectroscopia FT-IR, o equipamento varre uma amostra pixel por pixel — imagine uma câmera fotográfica, mas em vez de registrar cor, cada pixel registra um  **espectro de infravermelho completo** . O resultado é um objeto chamado **hipercubo espectral** (ou  *hyperspectral data cube* ):

```
Dimensões do hipercubo:
  x  →  coordenada espacial horizontal
  y  →  coordenada espacial vertical
  n  →  número de pontos espectrais (comprimentos de onda)
```

Cada posição $(x, y)$ da amostra gera um vetor de intensidades com $n$ pontos. Um experimento típico pode gerar de  **dezenas de milhares a milhões de espectros** , tornando a análise manual completamente inviável. É por isso que métodos multivariados são indispensáveis.

!!! tip "Por que pré-processar?"
Os espectros brutos raramente refletem apenas a composição química da amostra. Eles carregam junto ruído instrumental, variações de espessura, contribuições do ambiente (vapor d'água, CO₂) e artefatos físicos. O pré-processamento tem como objetivo **isolar o sinal químico real** antes de qualquer análise.

---

## 2. Pré-processamento Geral

### 2.1 Conversão de Transmitância para Absorbância

Alguns instrumentos armazenam os espectros em unidades de **transmitância** $T(\tilde{\nu})$, que representa a fração da luz que passou pela amostra. Para análises quantitativas, precisamos trabalhar em **absorbância** $A(\tilde{\nu})$, que é linearmente proporcional à concentração dos componentes (Lei de Lambert-Beer).

A conversão é simples:

$$
A(\tilde{\nu}) = -\log_{10} \left[ T(\tilde{\nu}) \right]
$$

!!! note "Intuição"
Transmitância de 100% significa que nada foi absorvido → absorbância = 0. Transmitância de 1% significa absorção quase total → absorbância = 2. A escala logarítmica lineariza a relação com a concentração.

---

### 2.2 Normalização

A intensidade de um espectro de absorbância FT-IR depende diretamente da **espessura da amostra** — amostras mais espessas absorvem mais luz em todos os comprimentos de onda. Se a espessura varia de pixel para pixel (o que é muito comum em amostras biológicas preparadas com micrótomo), os espectros ficam incomparáveis entre si.

A **normalização** remove esse efeito, colocando todos os espectros em uma escala comum.

#### 2.2.1 Normalização Min-Max

O valor mínimo do espectro é mapeado para 0 e o máximo para 1:

$$
A_{\text{norm}}(\tilde{\nu}) = \frac{A(\tilde{\nu}) - A_{\min}}{A_{\max} - A_{\min}}
$$

É simples, mas sensível a picos muito intensos ou outliers.

#### 2.2.2 Normalização Vetorial *(Vector Normalization)*

É a abordagem mais utilizada em espectroscopia de células e tecidos. O espectro é dividido por um fator de normalização obtido como a raiz quadrada da soma dos quadrados de todas as intensidades:

$$
A_{\text{norm}}(\tilde{\nu}) = \frac{A(\tilde{\nu})}{\sqrt{\sum_{\tilde{\nu}} A(\tilde{\nu})^2}}
$$

O resultado é um vetor com  **norma euclidiana igual a 1** . Frequentemente, o espectro é centralizado na média ( *mean-centered* ) antes da normalização.

!!! tip "Quando usar cada uma?"

- **Min-Max** : útil para visualização e comparação qualitativa rápida.
- **Vetorial** : preferida para análise multivariada (PCA, HCA, etc.), pois remove diferenças de escala sem distorcer as relações relativas entre bandas.

---

### 2.3 Redução de Ruído

O ruído é inevitável em qualquer medição instrumental. Em microespectroscopia FT-IR, o ruído é especialmente problemático na região de **alto número de onda** (acima de 3000 cm⁻¹), onde detectores de estado sólido têm sensibilidade reduzida.

#### 2.3.1 Suavização de Savitzky-Golay

O método de Savitzky-Golay é o algoritmo de suavização mais utilizado em espectroscopia. Ele aplica uma **janela deslizante** ao longo do espectro: em cada ponto, ajusta um polinômio de baixo grau aos pontos vizinhos e substitui o valor central pelo valor estimado pelo polinômio.

Matematicamente:

$$
Y^* *j = \frac{\sum* {i=-m}^{m} C_i \cdot Y_{j+i}}{N}
$$

onde:

| Símbolo  | Significado                        |
| --------- | ---------------------------------- |
| $Y^*_j$ | Valor suavizado no ponto$j$      |
| $C_i$   | Coeficientes da janela (tabelados) |
| $N$     | Constante de normalização        |
| $m$     | Metade da largura da janela        |

!!! warning "Cuidado com janelas muito largas"
Uma janela de suavização muito ampla remove ruído, mas também começa a distorcer os picos químicos, alargando e achatando bandas reais. Sempre verifique visualmente o espectro suavizado em relação ao original.

!!! info "Pré-requisito importante"
O espectro deve estar representado com **espaçamento linear em número de onda** antes de aplicar Savitzky-Golay. Caso contrário, a janela produzirá distorções.

#### 2.3.2 Análise de Componentes Principais Ajustada pelo Ruído (NA-PCA)


---

##### 2.3.2.1 O problema: ruído não é uniforme em espectros FT-IR

Antes de entender a NA-PCA, precisamos entender **por que a PCA convencional não é suficiente** para dados de microespectroscopia.

Em experimentos de FT-IR com detectores de estado sólido (como o HgCdTe), o ruído instrumental **não é constante ao longo do espectro**. Ele varia significativamente com o número de onda:

- Na região de **alto número de onda** (~4000 cm⁻¹): ruído muito maior, pois o detector tem sensibilidade reduzida próxima ao seu limiar de detecção (*band gap*).
- Na região de **baixo número de onda** (~1000 cm⁻¹): ruído menor, onde o detector opera com melhor eficiência.

Esse perfil de ruído dependente do comprimento de onda é chamado de **ruído colorido** (*colored noise*), em oposição ao **ruído branco** (*white noise*), que seria uniforme em todas as frequências.

!!! warning "Por que isso é um problema para a PCA convencional?"
    A PCA convencional não distingue variância química de variância instrumental. Se o ruído é muito maior em certas regiões espectrais, a PCA irá "desperdiçar" componentes principais tentando explicar esse ruído ao invés de capturar diferenças bioquímicas reais. Em casos extremos, os primeiros PCs podem ser dominados pelo perfil de ruído do instrumento, e não pela química da amostra.

A **NA-PCA (Noise-Adjusted Principal Component Analysis)** resolve exatamente esse problema: ela **mede o ruído do instrumento diretamente dos dados** e usa essa informação para reorganizar os componentes principais em ordem de **relação sinal-ruído decrescente**, não apenas de variância decrescente.

---

##### 2.3.2.2 Revisão rápida: PCA convencional

Para entender a NA-PCA, é útil revisar brevemente como a PCA convencional funciona.

Dado um conjunto de $n$ espectros organizados em uma matriz $\mathbf{S}_D$ (onde cada coluna é um espectro de $m$ pontos), a PCA calcula a **matriz de covariância do sinal**:

$$
\mathbf{C}_D = \mathbf{S}_D \mathbf{S}_D^T
$$

A diagonalização desta matriz fornece os **autovetores** (vetores de loading, ou componentes principais) e **autovalores** (variância explicada por cada PC):

$$
\mathbf{P}^T \mathbf{C}_D \mathbf{P} = \mathbf{\Lambda}_D
$$

Os PCs são ordenados do maior para o menor autovalor — ou seja, do que explica mais variância para o que explica menos. O problema é que **variância ≠ sinal químico**: um PC com alto autovalor pode estar descrevendo ruído instrumental intenso, não química real.

---

##### 2.3.2.3 A ideia central da NA-PCA

A NA-PCA foi originalmente desenvolvida pela comunidade de sensoriamento remoto (*remote sensing*) e adaptada para dados de microespectroscopia FT-IR por Bhargava e colaboradores.

A ideia central é elegante: **usar os próprios dados para medir o ruído do instrumento**.

Em um hipercubo espectral de microespectroscopia, nem todos os pixels contêm amostra. Alguns pixels correspondem a **regiões vazias** — fundo aquoso, espaço sem célula ou tecido, substrato puro. Esses pixels não carregam informação química da amostra, mas carregam **todo o ruído instrumental** presente no experimento.

Se separarmos os pixels em dois grupos:

- **Pixels de ruído** ($\mathbf{S}_N$): coletados em regiões sem amostra → representam o ruído puro do instrumento.
- **Pixels de sinal** ($\mathbf{S}_D$): coletados em regiões com célula/tecido → representam sinal + ruído.

Então podemos construir uma **matriz de covariância do ruído** $\mathbf{C}_N$ separadamente da **matriz de covariância do sinal** $\mathbf{C}_D$, e usar $\mathbf{C}_N$ para "calibrar" a PCA do sinal.

!!! tip "A intuição geométrica"
    Pense assim: a PCA convencional encontra as direções de maior variância total no espaço espectral. A NA-PCA encontra as direções de maior **relação sinal/ruído** — ou seja, as direções onde o sinal é grande *em relação ao* ruído. Essas são as direções mais informativas bioquimicamente.

---

##### 2.3.2.4 Derivação matemática passo a passo

###### Passo 1: Separação dos pixels

O primeiro passo é identificar quais pixels pertencem à amostra e quais pertencem ao fundo (ruído). Em dados de células vivas ou tecido, isso pode ser feito usando a intensidade integrada da **banda Amida I** (região 1626–1676 cm⁻¹):

- Pixels com alta intensidade Amida I → contêm célula → matriz $\mathbf{S}_D$
- Pixels com baixa intensidade Amida I → fundo/ruído → matriz $\mathbf{S}_N$

###### Passo 2: Construção das matrizes de covariância

Com os dois grupos separados, calculamos duas matrizes de covariância, ambas baseadas em vetores espectrais centralizados na média:

$$
\mathbf{C}_N = \mathbf{S}_N \mathbf{S}_N^T \quad \text{(covariância do ruído)}
$$

$$
\mathbf{C}_D = \mathbf{S}_D \mathbf{S}_D^T \quad \text{(covariância do sinal)}
$$

Ambas são matrizes $m \times m$, onde $m$ é o número de pontos espectrais.

###### Passo 3: Diagonalização da matriz de ruído

A matriz de covariância do ruído é diagonalizada:

$$
\mathbf{E}^T \mathbf{C}_N \mathbf{E} = \mathbf{\Delta}_N
$$

onde $\mathbf{E}$ é a matriz de autovetores do ruído e $\mathbf{\Delta}_N$ é a matriz diagonal de autovalores do ruído.

###### Passo 4: Transformação de branqueamento do ruído (*noise whitening*)

A ideia agora é construir uma transformação $\mathbf{F}$ que torne o ruído **uniforme em todas as direções** (ruído branco). Essa transformação é:

$$
\mathbf{F} = \mathbf{E} \cdot \mathbf{\Delta}_N^{-1/2}
$$

Pode-se verificar que ao aplicar $\mathbf{F}$, a matriz de covariância do ruído se torna a identidade:

$$
\mathbf{F}^T \mathbf{C}_N \mathbf{F} = \mathbf{\Delta}_N^{-1/2} \mathbf{E}^T \mathbf{\Delta}_N \mathbf{E} \mathbf{\Delta}_N^{-1/2} = \mathbf{I}
$$

!!! note "O que significa 'branquear o ruído'?"
    Após esta transformação, o ruído é matematicamente equivalente em todas as direções do espaço espectral. Isso significa que, quando aplicarmos a PCA a seguir, ela não será mais "atraída" pelas direções de alto ruído — pois todas as direções têm a mesma quantidade de ruído. O que sobrar de variância após o branqueamento será **sinal**.

###### Passo 5: PCA ajustada pelo ruído

Agora calculamos a **matriz de covariância do sinal ajustada pelo ruído**:

$$
\mathbf{C}_{\text{adj}} = \mathbf{F}^T \mathbf{C}_D \mathbf{F}
$$

Esta matriz é então diagonalizada:

$$
\mathbf{G}^T \mathbf{C}_{\text{adj}} \mathbf{G} = \mathbf{\Lambda}_{\text{adj}}
$$

Os autovalores $\mathbf{\Lambda}_{\text{adj}}$ agora expressam **relação sinal-ruído**, não apenas variância bruta.

###### Passo 6: Combinação das transformações

As matrizes $\mathbf{F}$ e $\mathbf{G}$ são combinadas em uma única matriz de transformação:

$$
\mathbf{H} = \mathbf{F} \cdot \mathbf{G}
$$

###### Passo 7: Cálculo dos componentes principais ajustados

Os componentes principais ajustados pelo ruído são calculados como:

$$
\mathbf{Z}'_{\text{adj}} = \mathbf{H}^T \mathbf{S}_D
$$

Esses componentes estão ordenados em ordem **decrescente de relação sinal-ruído**. Os primeiros componentes carregam o máximo de informação química com o mínimo de ruído; os últimos componentes são essencialmente ruído puro.

###### Passo 8: Reconstrução com redução de ruído

Finalmente, os espectros são reconstruídos usando apenas os **primeiros $p$ componentes ajustados** (tipicamente $p = 30$ a $50$):

$$
\hat{\mathbf{S}}_D = \sum_{k=1}^{p} z'_k \cdot h_k^T
$$

Os componentes de alta ordem (ruído) são simplesmente descartados. Os espectros reconstruídos têm **relação sinal-ruído significativamente melhorada**.

---

##### 2.3.2.5 Resumo esquemático do algoritmo

```
ENTRADA: Hipercubo espectral (x, y, n)
         │
         ▼
┌─────────────────────────────────────────────┐
│  SEGMENTAÇÃO DOS PIXELS                     │
│                                             │
│  Critério: intensidade integrada Amida I    │
│                                             │
│  → Pixels de RUÍDO    → matriz S_N          │
│  → Pixels de SINAL    → matriz S_D          │
└─────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────┐
│  CONSTRUÇÃO DAS COVARIÂNCIAS                │
│                                             │
│  C_N = S_N · S_N^T   (covariância do ruído) │
│  C_D = S_D · S_D^T   (covariância do sinal) │
└─────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────┐
│  BRANQUEAMENTO DO RUÍDO                     │
│                                             │
│  Diagonaliza C_N → E, Δ_N                  │
│  Calcula F = E · Δ_N^(-1/2)               │
│  (agora F^T · C_N · F = I)                 │
└─────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────┐
│  PCA AJUSTADA PELO RUÍDO                    │
│                                             │
│  C_adj = F^T · C_D · F                     │
│  Diagonaliza C_adj → G, Λ_adj              │
│  H = F · G                                 │
│  Z'_adj = H^T · S_D                        │
└─────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────┐
│  RECONSTRUÇÃO COM p COMPONENTES             │
│                                             │
│  Tipicamente p = 30 a 50                   │
│  Descarta componentes de alta ordem (ruído)│
│  Reconstrói espectros com S/N melhorado    │
└─────────────────────────────────────────────┘
         │
         ▼
  Espectros reconstruídos com
  alta relação sinal-ruído
```

---

##### 2.3.2.6 NA-PCA vs. suavização convencional

Por que usar NA-PCA em vez de simplesmente aplicar um filtro de suavização (como Savitzky-Golay)?

A tabela abaixo compara os dois métodos:

| Característica | Savitzky-Golay | NA-PCA |
|---|---|---|
| **Considera o perfil de ruído do instrumento?** | Não | Sim |
| **Trata ruído uniforme ao longo do espectro?** | Sim (aplica a mesma janela em todo o espectro) | Não necessário (adapta-se ao perfil real) |
| **Risco de distorcer picos químicos?** | Sim (janelas largas alargam e achatam bandas) | Baixo (quando $p$ é escolhido com cuidado) |
| **Melhoria típica em S/N** | Fator ~3 (janela de 9 pontos) | Fator ~3 a 10 |
| **Requer pixels sem amostra?** | Não | Sim |
| **Aplicável a espectros individuais?** | Sim | Não (requer conjunto de espectros) |

O ponto crítico é este: em instrumentos FT-IR como o Perkin Elmer Spotlight 400, **o ruído na região de 4000 cm⁻¹ é muito maior do que na região de 1000 cm⁻¹**. Um filtro de Savitzky-Golay com janela suficientemente larga para reduzir o ruído em 4000 cm⁻¹ irá **supersuavizar** a região de 1000 cm⁻¹ (que já tem bom S/N), distorcendo as bandas químicas importantes. A NA-PCA aplica, de forma implícita, uma suavização **diferente em cada região espectral**, proporcional ao ruído presente naquela região.

!!! example "Resultado prático (Marcsisin et al., 2012)"
    Usando 30 componentes principais ajustados na reconstrução, a NA-PCA reduziu o nível de ruído em células HeLa vivas por um **fator de 3** — equivalente a uma suavização Savitzky-Golay de 9 pontos — mas sem distorcer a forma das bandas.

    A S/N média dos espectros celulares saltou de **475 (imagem bruta)** para **1254 (após NA-PCA)**, superando inclusive os **647** obtidos com o modo de ponto com 16× mais interferogramas co-adicionados.

---

##### 2.3.2.7 Um detalhe importante: os "pixels vazios" são valiosos

Em microespectroscopia de células vivas ou tecido, é comum que uma fração significativa da imagem corresponda a **regiões sem amostra** — fundo aquoso, espaço entre células, substrato de CaF₂.

Intuitivamente, poderíamos pensar que esses pixels são um "desperdício" de tempo de aquisição e espaço de armazenamento. A NA-PCA inverte essa perspectiva:

> **Os pixels vazios não são desperdício — eles são a medição do ruído do instrumento.**

Esses pixels fornecem uma estimativa precisa e contemporânea (coletada durante o mesmo experimento, nas mesmas condições instrumentais) do perfil de ruído do detector. Essa informação é usada para construir a matriz $\mathbf{C}_N$ e, consequentemente, para realizar o branqueamento do ruído.

!!! tip "Implicação prática"
    Ao planejar experimentos de microespectroscopia FT-IR para uso com NA-PCA, é vantajoso garantir que a imagem coletada inclua uma **área de fundo limpo e representativo** — sem contaminação por amostra — suficientemente grande para estimar o perfil de ruído com precisão estatística.

---

##### 2.3.2.8 Escolha do número de componentes $p$

A qualidade da reconstrução depende criticamente de quantos componentes ajustados são retidos. Algumas diretrizes práticas:

- **Poucos componentes** ($p$ muito pequeno): reconstrução suavizada demais, perde informação química real, especialmente de bandas menores.
- **Muitos componentes** ($p$ muito grande): inclui componentes de ruído, deteriorando a S/N.
- **Valor típico**: $p = 30$ a $50$ componentes costuma ser suficiente para capturar mais de 99% da variância química em hipercubos de células.

Uma forma objetiva de escolher $p$ é inspecionar os **autovalores ajustados** $\mathbf{\Lambda}_{\text{adj}}$: valores maiores que 1 indicam que aquele componente tem mais sinal do que ruído; valores próximos a 1 indicam que o componente é dominado por ruído.

---

##### 2.3.2.9 Aplicação em células vivas: um exemplo concreto

O trabalho de Marcsisin et al. (2012) ilustra perfeitamente o impacto da NA-PCA em um cenário de alta relevância biomédica: monitoramento de células HeLa vivas por FT-IR em ambiente aquoso.

###### O desafio

Células vivas precisam ser medidas em solução aquosa, o que impõe:

- **Caminho óptico muito curto** (6–10 μm) para minimizar a absorção da água líquida.
- **Poucas varreduras co-adicionadas** por pixel (apenas 8) para completar a imagem em tempo hábil (~30 minutos), antes que as células sofram alterações biológicas.
- **Resultado**: espectros com S/N mais baixa do que em amostras fixas/secas.

###### A solução com NA-PCA

Com o modo de imagem + NA-PCA, o fluxo de trabalho foi:

1. Coletar um mapa de 1 mm × 1 mm (25.600 espectros, 8 interferogramas/pixel) em ~30 minutos.
2. Identificar pixels de ruído (fundo aquoso) e pixels de sinal (células) pela intensidade Amida I.
3. Aplicar NA-PCA com $p = 30$ componentes → redução do ruído por fator 3.
4. Aplicar o algoritmo PapMap para co-adicionar pixels pertencentes à mesma célula → espectro médio de cada célula individual.
5. Converter para segunda derivada (Savitzky-Golay 9 pontos) + normalização vetorial.
6. Análise por PCA para comparar células tratadas vs. não tratadas.

###### O resultado

Com este fluxo, foi possível **detectar diferenças espectrais entre células HeLa não tratadas e células tratadas com ciclofosfamida** (um agente quimioterápico alquilante). As diferenças foram observadas principalmente:

- Na posição da banda **Amida I** (~1652 cm⁻¹): variações em proteínas.
- No estiramento simétrico **PO₂⁻** (~1080 cm⁻¹): variações no DNA/RNA.

Esses resultados são consistentes com o mecanismo de ação conhecido da ciclofosfamida: danos ao DNA por alquilação das bases guanina.

---


---

### 2.4 Cálculo de Derivadas Espectrais

Uma das transformações mais importantes em pré-processamento de espectros FT-IR é o cálculo da  **segunda derivada** . Essa operação oferece três benefícios simultâneos:

1. **Remove linhas de base inclinadas** — derivadas de ordens 1 e 2 eliminam offsets constantes e lineares.
2. **Aumenta a resolução aparente** — bandas sobrepostas que parecem uma só tornam-se distinguíveis na derivada.
3. **Reduz o efeito de onda estacionária** — efeito físico que distorce espectros em substratos metálicos.

#### Como interpretar a segunda derivada?

Considere um pico de absorbância com formato gaussiano. Sua segunda derivada é uma função com um **mínimo** na posição do pico original. Onde havia um máximo positivo, agora há um mínimo negativo — os picos ficam "invertidos".

Quando duas bandas se sobrepõem fortemente, a absorbância bruta mostra apenas um único pico largo. A segunda derivada, por ser sensível à **curvatura** do espectro, consegue revelar os dois mínimos separados, mesmo que a separação seja menor que a largura das bandas.

!!! warning "Custo da derivada: relação sinal-ruído"
A segunda derivada amplifica o ruído de alta frequência. Por isso, ela deve sempre ser combinada com uma etapa de suavização (ex.: Savitzky-Golay) aplicada antes ou durante o cálculo da derivada. Há um compromisso intrínseco: mais suavização preserva melhor o sinal, mas pode mascarar detalhes espectrais finos.

As derivadas podem ser calculadas pelos coeficientes específicos do algoritmo de Savitzky-Golay para derivadas, ou no espaço de Fourier multiplicando o interferograma por uma função parabólica centrada no ponto de diferença de caminho zero (ZPD).

---

## 3. Efeitos Confundidores em FT-IR

Além do ruído instrumental, existem contribuições físicas e químicas que se superpõem ao sinal de interesse e precisam ser tratadas especificamente. Dois dos mais importantes em amostras biológicas são o **vapor d'água** e a  **parafina** .

### 3.1 Vapor d'Água

#### Por que é um problema?

A molécula de água ($\mathrm{H_2O}$) possui modos rotacionais de baixa energia que, à temperatura ambiente, geram um **espectro roto-vibracional** muito complexo. Esse espectro ocupa duas regiões amplas:

* **3200 – 4000 cm⁻¹** : região do estiramento O-H
* **1250 – 2000 cm⁻¹** : região da deformação H-O-H

O problema é grave porque **a região de 1250–2000 cm⁻¹ coincide exatamente com as bandas Amida I (~1650 cm⁻¹) e Amida II (~1550 cm⁻¹)** de proteínas — as bandas mais informativas em tecidos biológicos.

#### Como surge a contaminação?

Os espectrômetros FT-IR são purgados com ar seco ou operados em vácuo para minimizar o vapor d'água. Contudo, em experimentos de longa duração (como aquisição de um mapa de tecido que pode levar horas), a quantidade de vapor d'água no instrumento pode variar ligeiramente ao longo do tempo. Quando isso ocorre:

* Espectros coletados no início da medição têm uma quantidade de vapor d'água.
* Espectros coletados no final têm uma quantidade diferente.
* Ao dividir pelo espectro de fundo ( *background* ), surgem **contribuições positivas ou negativas** de vapor d'água nos espectros da amostra.

Mesmo quando essas contribuições são invisíveis a olho nu, elas dominam os  **componentes de alta ordem em PCA** , mascarando diferenças químicas reais.

#### Estratégias de correção

**Estratégia 1 — Subtração simples escalada:**
Medir um espectro puro de vapor d'água e subtraí-lo de cada espectro da amostra com um fator de escala. A dificuldade é que o espectro de vapor d'água é **sensível à temperatura** — variações de temperatura deslocam os máximos dos ramos P e R do contorno rotacional — e o fator de escala correto é difícil de determinar na presença de ruído.

**Estratégia 2 — Correção multivariada por EMSC:**
A abordagem mais robusta foi desenvolvida por Bruun et al. (2006). Consiste em:

1. Coletar espectros de vapor d'água em diferentes condições de temperatura e concentração.
2. Aplicar PCA sobre esses espectros para extrair as **componentes principais de variação** do vapor d'água: o espectro médio e o primeiro espectro de variância.
3. Incluir esses dois vetores no modelo EMSC como espectros constituintes de interferência.
4. O modelo estima e remove a contribuição do vapor d'água de cada espectro da amostra automaticamente.

Este procedimento reduz os resíduos de vapor d'água em mais de uma ordem de grandeza.

!!! example "Impacto prático"
Antes da correção, a análise de componentes principais de espectros de tecido frequentemente mostra que os primeiros PCs de alta ordem são dominados por padrões de linhas finas característicos do vapor d'água. Após a correção, esses PCs passam a refletir diferenças bioquímicas reais entre as células.

---

### 3.2 Parafina

#### Por que é um problema?

Em histologia clínica, tecidos biológicos são rotineiramente **embebidos em parafina** (FFPE —  *Formalin-Fixed Paraffin-Embedded* ) para facilitar o corte em micrótomo e o armazenamento a longo prazo. Isso significa que a grande maioria dos espécimes de biópsia disponíveis nos bancos de patologia hospitalar estão em parafina.

A parafina é uma mistura de alcanos de cadeia longa ($\mathrm{C_nH_{2n+2}}$, com $n$ entre 20 e 40). Seus modos vibracionais geram bandas de absorbância **muito intensas** nas seguintes regiões:

| Região (cm⁻¹) | Atribuição                                                    |
| ---------------- | --------------------------------------------------------------- |
| 2920 e 2850      | Estiramento C-H assimétrico e simétrico ($-\mathrm{CH_2}-$) |
| 1465             | Deformação C-H no plano                                       |
| 720              | Deformação C-H fora do plano (rocking)                        |

O problema crítico é que as bandas em **~2920 e ~2850 cm⁻¹ sobrepõem-se às bandas de ácidos graxos** dos lipídios celulares, e a banda em  **1465 cm⁻¹ interfere com bandas de proteínas e lipídios** . Isso impede a análise direta dos espectros sem remoção prévia da parafina.

#### Estratégias de remoção

**Estratégia 1 — Deparafinização química:**
O tecido é submetido a banhos sucessivos em xileno (ou substitutos como Histo-Clear), que dissolve a parafina. Embora eficaz, esse processo é demorado, requer manuseio de solventes orgânicos e pode alterar a composição lipídica da amostra.

**Estratégia 2 — Subtração espectral digital:**
Um espectro de referência de parafina pura é medido separadamente. Para cada espectro do mapa, a contribuição da parafina é estimada (normalmente usando a intensidade das bandas em 2920 ou 2850 cm⁻¹, que são exclusivas da parafina) e subtraída com o fator de escala adequado:

$$
A_{\text{corr}}(\tilde{\nu}) = A_{\text{amostra}}(\tilde{\nu}) - \alpha \cdot A_{\text{parafina}}(\tilde{\nu})
$$

onde $\alpha$ é determinado minimizando o resíduo na região das bandas de parafina.

**Estratégia 3 — Correção por EMSC:**
Analogamente ao vapor d'água, o espectro de diferença da parafina pode ser incluído como constituinte no modelo EMSC, permitindo estimar e remover sua contribuição de forma mais robusta e simultânea a outras correções.

!!! warning "Atenção à região do C-H"
Mesmo após a remoção computacional da parafina, a região 2800–3000 cm⁻¹ pode apresentar resíduos. Muitos estudos optam por **excluir essa região** da análise multivariada quando a remoção não é perfeitamente satisfatória, focando apenas na região de impressão digital (800–1800 cm⁻¹).

---

## 4. Métodos Não-Supervisionados de Segmentação

Após o pré-processamento, o objetivo é extrair informação química e estrutural do hipercubo espectral. Os **métodos não-supervisionados** fazem isso sem nenhuma informação prévia sobre o resultado esperado — apenas a variância intrínseca dos espectros é utilizada.

!!! info "Supervisionado vs. Não-supervisionado"

- **Não-supervisionado** : o algoritmo agrupa espectros apenas pela sua similaridade matemática. Não precisa de "respostas certas" para treinar. Ideal para exploração e descoberta.
- **Supervisionado** : o algoritmo é treinado com espectros de classes conhecidas e depois classifica espectros desconhecidos. Requer dados rotulados.

Nesta seção, focamos nos três principais métodos não-supervisionados:  **PCA** , **HCA** e  **KMCA** .

---

### 4.1 Análise de Componentes Principais (PCA)

#### A ideia central

Imagine um conjunto de 10.000 espectros, cada um com 500 pontos de intensidade. Matematicamente, cada espectro é um  **ponto em um espaço de 500 dimensões** . Isso é impossível de visualizar, mas a PCA resolve isso encontrando as **direções de maior variância** nesse espaço e projetando os dados nessas direções.

O resultado são novos eixos chamados **componentes principais (PCs)** ou  **vetores de loading** , ordenados do mais ao menos informativo. Tipicamente, os primeiros 10-15 PCs capturam mais de 99% de toda a variância do conjunto de dados.

#### Construção matemática

O conjunto completo de $n$ espectros é organizado em uma matriz $\mathbf{S}$, onde cada coluna é um espectro. A **matriz de correlação** é calculada como:

$$
\mathbf{C} = \mathbf{S} \cdot \mathbf{S}^T
$$

Os elementos $C_{kl}$ representam a correlação entre as intensidades nos números de onda $\tilde{\nu}_k$ e $\tilde{\nu}_l$, somada sobre todos os espectros.

A **diagonalização** da matriz de correlação fornece os autovetores (vetores de loading) e autovalores:

$$
\mathbf{P}^T \mathbf{C} \mathbf{P} = \mathbf{\Lambda}
$$

Os autovalores $\Lambda$ expressam a **variância explicada** por cada PC. Os **scores** $\alpha$ descrevem quanto cada PC contribui para cada espectro:

$$
S' *i(\tilde{\nu}) = \sum* {j=1}^{p} \alpha_{ij} \cdot Z_j(\tilde{\nu})
$$

#### O que os PCs representam?

* Os **primeiros PCs** capturam as maiores fontes de variação no conjunto de dados. Em espectros de tecido, isso geralmente corresponde a diferenças entre tipos celulares.
* **PCs intermediários** podem capturar interferências como vapor d'água ou parafina residual.
* **PCs de alta ordem** geralmente contêm apenas ruído não correlacionado.

#### Visualização: o *scores plot*

Plotando os scores de dois PCs entre si (ex.: PC1 × PC2), cada ponto representa um espectro. Espectros quimicamente semelhantes aparecem próximos; espectros diferentes, afastados. Essa visualização em espaço de **dimensionalidade muito reduzida** (2-3 dimensões vs. 500 originais) é uma das ferramentas mais poderosas para exploração inicial de dados.

#### PCA como imagem

Os scores de cada PC podem ser mapeados de volta às coordenadas espaciais $(x, y)$, gerando uma  **imagem de scores** . Regiões da amostra com composição química similar aparecem com valores de score similares. Até três imagens de scores podem ser sobrepostas em canais R, G, B para gerar um mapa pseudocor.

!!! tip "PCA como redução de ruído"
Os primeiros 15-30 PCs contêm essencialmente todo o sinal relevante. Reconstruindo os espectros usando apenas esses PCs (descartando os de alta ordem), obtém-se uma **redução significativa de ruído** sem perda apreciável de informação química.

---

### 4.2 Análise de Agrupamento Hierárquico (HCA)

#### A ideia central

Enquanto a PCA decompõe os dados em componentes matemáticos, a HCA tem um objetivo diferente: **agrupar espectros por similaridade** e gerar uma imagem pseudocor onde pixels com espectros parecidos recebem a mesma cor.

A HCA é  **completamente não-supervisionada** : nenhuma informação sobre tipos celulares ou estrutura tecidual é fornecida ao algoritmo. A imagem emerge puramente da matemática.

#### Medidas de similaridade

A HCA começa calculando a  **similaridade entre todos os pares de espectros** . As métricas mais usadas são:

**Coeficiente de correlação de Pearson:**

$$
C' *{ij} = \sum* {N=1}^{n} S^i(\tilde{\nu}_N) \cdot S^j(\tilde{\nu}_N)
$$

Quando $C' *{ij} \approx 1$: espectros muito similares. Quando $C'* {ij} \approx 0$: espectros muito diferentes.

**Distância Euclidiana:**

$$
D_{ij} = \sqrt{\sum_{N=1}^{n} \left[ S^i(\tilde{\nu}_N) - S^j(\tilde{\nu}_N) \right]^2}
$$

Distância pequena = espectros parecidos. Distância grande = espectros diferentes.

!!! warning "Custo computacional"
Para um mapa de 200 × 200 pixels (40.000 espectros), a matriz de correlação tem dimensão 40.000 × 40.000, com 1,6 × 10⁹ entradas e requer cerca de  **6 GB de memória** . Por isso, HCA é prático apenas para conjuntos de até ~100.000 espectros.

#### O algoritmo de agrupamento

Uma vez calculadas as similaridades:

1. **Inicialização** : cada espectro é tratado como um cluster individual.
2. **Fusão** : os dois espectros (ou clusters) mais similares são fundidos em um novo cluster.
3. **Atualização** : a similaridade do novo cluster com todos os demais é recalculada.
4. **Repetição** : o processo se repete até que todos os espectros pertençam a um único cluster.

O método de **Ward** é o mais utilizado para a fusão, pois minimiza o aumento da variância intra-cluster e produz grupos mais homogêneos.

#### Dendrograma e imagem pseudocor

O processo de fusão gera um **dendrograma** — uma árvore hierárquica que mostra em que nível de similaridade cada fusão ocorreu. Cortando o dendrograma em um determinado nível, obtemos um número desejado de clusters (tipicamente 5-8 para dados de tecido).

Cada cluster recebe uma cor arbitrária, e as coordenadas espaciais dos espectros são mapeadas nessa cor, gerando a  **imagem pseudocor de HCA** . O espectro médio de cada cluster pode ser calculado para representar a composição química daquela região.

!!! example "Resultado típico"
Em tecido de tumor, os clusters do HCA geralmente correspondem a: epitélio tumoral, estroma, músculo, necrose e células inflamatórias — tudo isso sem qualquer informação prévia fornecida ao algoritmo.

---

### 4.3 Análise de Agrupamento K-Médias (KMCA)

#### A ideia central

A KMCA ( *K-Means Cluster Analysis* ) é uma alternativa computacionalmente mais eficiente ao HCA, especialmente adequada para **grandes conjuntos de dados** (>100.000 espectros). Em vez de construir uma hierarquia, ela parte diretamente do número $k$ de clusters desejado e distribui os espectros entre eles.

#### O algoritmo passo a passo

**Passo 1 — Inicialização aleatória:**
$k$ espectros são selecionados aleatoriamente como **centroides iniciais** de cada cluster.

**Passo 2 — Atribuição:**
Cada espectro é atribuído ao cluster cujo centroide é o mais próximo (menor distância):

$$
\min_k \sum_j \left( S^j_i - m_k \right)^2
$$

**Passo 3 — Atualização:**
Os centroides são recalculados como a **média de todos os espectros** atribuídos a cada cluster.

**Passo 4 — Convergência:**
Os passos 2 e 3 são repetidos até que nenhum espectro mude de cluster entre iterações consecutivas, ou até que a melhoria na homogeneidade intra-cluster caia abaixo de um limite pré-definido.

#### Vantagens e limitações

| Aspecto                          | KMCA                                               | HCA                                   |
| -------------------------------- | -------------------------------------------------- | ------------------------------------- |
| **Velocidade**             | Rápida                                            | Lenta (cresce quadraticamente)        |
| **Tamanho do dado**        | Grandes conjuntos (>100k)                          | Limitada (~100k espectros)            |
| **Número de clusters**    | Deve ser definido pelo usuário                    | Determinado pelo dendrograma          |
| **Determinismo**           | Resultado pode variar (inicialização aleatória) | Resultado é determinístico          |
| **Qualidade dos clusters** | Boa                                                | Excelente (clusters mais homogêneos) |

!!! tip "Qual usar na prática?"

- **HCA** : preferida para dados FT-IR de alta relação S/N, onde espectros de segunda derivada podem ser utilizados e o conjunto tem até ~100.000 pixels.
- **KMCA** : preferida para conjuntos muito grandes ou quando o tempo de computação é limitante.
- Ambas produzem **espectros médios de cluster** de altíssima qualidade, que podem ser interpretados bioquimicamente.

---

## 5. Fluxo de Trabalho Completo

O diagrama abaixo resume a sequência típica de processamento de um hipercubo FT-IR de tecido biológico:

```
Dados brutos (hipercubo x, y, n)
         │
         ▼
┌─────────────────────────────────────────┐
│         PRÉ-PROCESSAMENTO               │
│                                         │
│  1. Conversão Transmitância → Absorbância│
│  2. Remoção de vapor d'água (EMSC)      │
│  3. Remoção de parafina (subtração)     │
│  4. Normalização vetorial               │
│  5. Suavização Savitzky-Golay           │
│  6. Cálculo da 2ª derivada             │
└─────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────┐
│     PRÉ-SEGMENTAÇÃO (PCA)               │
│                                         │
│  • Redução de dimensionalidade          │
│  • Identificação de fontes de variância │
│  • Redução de ruído                     │
│  • Scores plots para exploração         │
└─────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────┐
│    SEGMENTAÇÃO NÃO-SUPERVISIONADA       │
│                                         │
│  HCA (dados menores, alta qualidade)    │
│       ou                                │
│  KMCA (dados grandes, alta velocidade)  │
│                                         │
│  → Imagem pseudocor                     │
│  → Espectros médios de cluster          │
└─────────────────────────────────────────┘
         │
         ▼
   Interpretação bioquímica
   e análise supervisionada
```

---

## 6. Referências

> Diem, M. (2015).  **Data Preprocessing and Data Processing in Microspectral Analysis** . In:  *Modern Vibrational Spectroscopy and Micro-Spectroscopy: Theory, Instrumentation and Biomedical Applications* . John Wiley & Sons, Ltd. Chapter 12, pp. 251–282.

> Bruun, S.W., Kohler, A., Adt, I., et al. (2006). Correcting attenuated total reflection-Fourier transform infrared spectra for water vapor and carbon dioxide.  *Applied Spectroscopy* , 60(9), 1029–1039.

> Savitzky, A., Golay, M.J.E. (1964). Smoothing and differentiation of data by simplified least-squares procedures.  *Analytical Chemistry* , 36(8), 1627–1639.

> Kohler, A., et al. (2005). Extended multiplicative signal correction as a tool for separation and characterization of physical and chemical information in FT-IR microscopy images.  *Applied Spectroscopy* , 59, 707–716.

> Ward, J.H. (1963). Hierarchical grouping to optimize an objective function.  *Journal of the American Statistical Association* , 58(301), 236–244.

> Marcsisin, E.J., Uttero, C.M., Mazur, A.I., Miljkovic, M., Bird, B., Diem, M. (2012). **Noise Adjusted Principal Component reconstruction to optimize infrared microspectroscopy of individual live cells**. *Analyst*, 137, 2958–2964. DOI: 10.1039/c2an15868j

> Bhargava, R., Ribar, T., Koenig, J.L. (1999). Towards faster FT-IR imaging by reducing noise. *Applied Spectroscopy*, 53(11), 1313–1322.

> Green, A.A., Berman, M., Switzer, P., Craig, M.D. (1988). A transformation for ordering multispectral data in terms of image quality with implications for noise removal. *IEEE Transactions on Geoscience and Remote Sensing*, 26(1), 65–74.

> Diem, M. (2015). Data Preprocessing and Data Processing in Microspectral Analysis. In: *Modern Vibrational Spectroscopy and Micro-Spectroscopy*. John Wiley & Sons. Chapter 12.
