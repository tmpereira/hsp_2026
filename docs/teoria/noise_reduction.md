> ⚠️ **Documentação em construção.** Algumas informações podem estar imprecisas ou incompletas.

> **Para quem é este texto?** Este material foi escrito para estudantes de graduação em Engenharia e Ciências Exatas que estão tendo seu primeiro contato com análise de dados espectrais. O objetivo é construir a intuição passo a passo, do dado bruto até a imagem química final.

---

## 1. Redução de Ruído

O ruído é inevitável em qualquer medição instrumental. Em microespectroscopia FT-IR, o ruído é especialmente problemático na região de **alto número de onda** (acima de 3000 cm⁻¹), onde detectores de estado sólido têm sensibilidade reduzida.

### 1.1 Suavização de Savitzky-Golay

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

### 1.2 Análise de Componentes Principais Ajustada pelo Ruído (NA-PCA)


---

#### 1.2.1 O problema: ruído não é uniforme em espectros FT-IR

Antes de entender a NA-PCA, precisamos entender **por que a PCA convencional não é suficiente** para dados de microespectroscopia.

Em experimentos de FT-IR com detectores de estado sólido (como o HgCdTe), o ruído instrumental **não é constante ao longo do espectro**. Ele varia significativamente com o número de onda:

- Na região de **alto número de onda** (~4000 cm⁻¹): ruído muito maior, pois o detector tem sensibilidade reduzida próxima ao seu limiar de detecção (*band gap*).
- Na região de **baixo número de onda** (~1000 cm⁻¹): ruído menor, onde o detector opera com melhor eficiência.

Esse perfil de ruído dependente do comprimento de onda é chamado de **ruído colorido** (*colored noise*), em oposição ao **ruído branco** (*white noise*), que seria uniforme em todas as frequências.

!!! warning "Por que isso é um problema para a PCA convencional?"
    A PCA convencional não distingue variância química de variância instrumental. Se o ruído é muito maior em certas regiões espectrais, a PCA irá "desperdiçar" componentes principais tentando explicar esse ruído ao invés de capturar diferenças bioquímicas reais. Em casos extremos, os primeiros PCs podem ser dominados pelo perfil de ruído do instrumento, e não pela química da amostra.

A **NA-PCA (Noise-Adjusted Principal Component Analysis)** resolve exatamente esse problema: ela **mede o ruído do instrumento diretamente dos dados** e usa essa informação para reorganizar os componentes principais em ordem de **relação sinal-ruído decrescente**, não apenas de variância decrescente.

---

#### 1.2.2 Revisão rápida: PCA convencional

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

#### 1.2.3 A ideia central da NA-PCA

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

#### 1.2.4 Derivação matemática passo a passo

##### 1.2.4.1 Passo 1: Separação dos pixels

O primeiro passo é identificar quais pixels pertencem à amostra e quais pertencem ao fundo (ruído). Em dados de células vivas ou tecido, isso pode ser feito usando a intensidade integrada da **banda Amida I** (região 1626–1676 cm⁻¹):

- Pixels com alta intensidade Amida I → contêm célula → matriz $\mathbf{S}_D$
- Pixels com baixa intensidade Amida I → fundo/ruído → matriz $\mathbf{S}_N$

##### 1.2.4.2 Passo 2: Construção das matrizes de covariância

Com os dois grupos separados, calculamos duas matrizes de covariância, ambas baseadas em vetores espectrais centralizados na média:

$$
\mathbf{C}_N = \mathbf{S}_N \mathbf{S}_N^T \quad \text{(covariância do ruído)}
$$

$$
\mathbf{C}_D = \mathbf{S}_D \mathbf{S}_D^T \quad \text{(covariância do sinal)}
$$

Ambas são matrizes $m \times m$, onde $m$ é o número de pontos espectrais.

##### 1.2.4.3 Passo 3: Diagonalização da matriz de ruído

A matriz de covariância do ruído é diagonalizada:

$$
\mathbf{E}^T \mathbf{C}_N \mathbf{E} = \mathbf{\Delta}_N
$$

onde $\mathbf{E}$ é a matriz de autovetores do ruído e $\mathbf{\Delta}_N$ é a matriz diagonal de autovalores do ruído.

##### 1.2.4.4 Passo 4: Transformação de branqueamento do ruído (*noise whitening*)

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

##### 1.2.4.5 Passo 5: PCA ajustada pelo ruído

Agora calculamos a **matriz de covariância do sinal ajustada pelo ruído**:

$$
\mathbf{C}_{\text{adj}} = \mathbf{F}^T \mathbf{C}_D \mathbf{F}
$$

Esta matriz é então diagonalizada:

$$
\mathbf{G}^T \mathbf{C}_{\text{adj}} \mathbf{G} = \mathbf{\Lambda}_{\text{adj}}
$$

Os autovalores $\mathbf{\Lambda}_{\text{adj}}$ agora expressam **relação sinal-ruído**, não apenas variância bruta.

##### 1.2.4.6 Passo 6: Combinação das transformações

As matrizes $\mathbf{F}$ e $\mathbf{G}$ são combinadas em uma única matriz de transformação:

$$
\mathbf{H} = \mathbf{F} \cdot \mathbf{G}
$$

##### 1.2.4.7 Passo 7: Cálculo dos componentes principais ajustados

Os componentes principais ajustados pelo ruído são calculados como:

$$
\mathbf{Z}'_{\text{adj}} = \mathbf{H}^T \mathbf{S}_D
$$

Esses componentes estão ordenados em ordem **decrescente de relação sinal-ruído**. Os primeiros componentes carregam o máximo de informação química com o mínimo de ruído; os últimos componentes são essencialmente ruído puro.

##### 1.2.4.8 Passo 8: Reconstrução com redução de ruído

Finalmente, os espectros são reconstruídos usando apenas os **primeiros $p$ componentes ajustados** (tipicamente $p = 30$ a $50$):

$$
\hat{\mathbf{S}}_D = \sum_{k=1}^{p} z'_k \cdot h_k^T
$$

Os componentes de alta ordem (ruído) são simplesmente descartados. Os espectros reconstruídos têm **relação sinal-ruído significativamente melhorada**.

---

#### 1.2.5 Resumo esquemático do algoritmo

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

#### 1.2.6 NA-PCA vs. suavização convencional

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

#### 1.2.7 Um detalhe importante: os "pixels vazios" são valiosos

Em microespectroscopia de células vivas ou tecido, é comum que uma fração significativa da imagem corresponda a **regiões sem amostra** — fundo aquoso, espaço entre células, substrato de CaF₂.

Intuitivamente, poderíamos pensar que esses pixels são um "desperdício" de tempo de aquisição e espaço de armazenamento. A NA-PCA inverte essa perspectiva:

> **Os pixels vazios não são desperdício — eles são a medição do ruído do instrumento.**

Esses pixels fornecem uma estimativa precisa e contemporânea (coletada durante o mesmo experimento, nas mesmas condições instrumentais) do perfil de ruído do detector. Essa informação é usada para construir a matriz $\mathbf{C}_N$ e, consequentemente, para realizar o branqueamento do ruído.

!!! tip "Implicação prática"
    Ao planejar experimentos de microespectroscopia FT-IR para uso com NA-PCA, é vantajoso garantir que a imagem coletada inclua uma **área de fundo limpo e representativo** — sem contaminação por amostra — suficientemente grande para estimar o perfil de ruído com precisão estatística.

---

#### 1.2.8 Escolha do número de componentes $p$

A qualidade da reconstrução depende criticamente de quantos componentes ajustados são retidos. Algumas diretrizes práticas:

- **Poucos componentes** ($p$ muito pequeno): reconstrução suavizada demais, perde informação química real, especialmente de bandas menores.
- **Muitos componentes** ($p$ muito grande): inclui componentes de ruído, deteriorando a S/N.
- **Valor típico**: $p = 30$ a $50$ componentes costuma ser suficiente para capturar mais de 99% da variância química em hipercubos de células.

Uma forma objetiva de escolher $p$ é inspecionar os **autovalores ajustados** $\mathbf{\Lambda}_{\text{adj}}$: valores maiores que 1 indicam que aquele componente tem mais sinal do que ruído; valores próximos a 1 indicam que o componente é dominado por ruído.

---

#### 1.2.9 Aplicação em células vivas: um exemplo concreto

O trabalho de Marcsisin et al. (2012) ilustra perfeitamente o impacto da NA-PCA em um cenário de alta relevância biomédica: monitoramento de células HeLa vivas por FT-IR em ambiente aquoso.

##### 1.2.9.1 O desafio

Células vivas precisam ser medidas em solução aquosa, o que impõe:

- **Caminho óptico muito curto** (6–10 μm) para minimizar a absorção da água líquida.
- **Poucas varreduras co-adicionadas** por pixel (apenas 8) para completar a imagem em tempo hábil (~30 minutos), antes que as células sofram alterações biológicas.
- **Resultado**: espectros com S/N mais baixa do que em amostras fixas/secas.

##### 1.2.9.2 A solução com NA-PCA

Com o modo de imagem + NA-PCA, o fluxo de trabalho foi:

1. Coletar um mapa de 1 mm × 1 mm (25.600 espectros, 8 interferogramas/pixel) em ~30 minutos.
2. Identificar pixels de ruído (fundo aquoso) e pixels de sinal (células) pela intensidade Amida I.
3. Aplicar NA-PCA com $p = 30$ componentes → redução do ruído por fator 3.
4. Aplicar o algoritmo PapMap para co-adicionar pixels pertencentes à mesma célula → espectro médio de cada célula individual.
5. Converter para segunda derivada (Savitzky-Golay 9 pontos) + normalização vetorial.
6. Análise por PCA para comparar células tratadas vs. não tratadas.

##### 1.2.9.3 O resultado

Com este fluxo, foi possível **detectar diferenças espectrais entre células HeLa não tratadas e células tratadas com ciclofosfamida** (um agente quimioterápico alquilante). As diferenças foram observadas principalmente:

- Na posição da banda **Amida I** (~1652 cm⁻¹): variações em proteínas.
- No estiramento simétrico **PO₂⁻** (~1080 cm⁻¹): variações no DNA/RNA.

Esses resultados são consistentes com o mecanismo de ação conhecido da ciclofosfamida: danos ao DNA por alquilação das bases guanina.

---

