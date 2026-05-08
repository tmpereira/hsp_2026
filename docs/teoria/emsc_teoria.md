# Correção Multiplicativa de Sinal Estendida (EMSC)

> **Para quem é este texto?** Este material foi escrito para estudantes de graduação em Engenharia e Ciências Exatas que estão tendo seu primeiro contato com técnicas de pré-processamento de espectros vibracionais. Nenhum conhecimento prévio de espectroscopia avançada é necessário — apenas álgebra linear básica e um pouco de curiosidade!

---

## 1. Por que precisamos pré-processar espectros?

Imagine que você está tentando medir a concentração de proteínas em uma amostra de tecido biológico usando espectroscopia infravermelha (FT-IR) ou Raman. Você posiciona a amostra no equipamento, coleta o espectro e... ele não parece exatamente o que você esperava.

Por quê? Porque o espectro medido raramente contém **apenas** a informação química que você quer. Ele carrega junto uma série de outros efeitos indesejados:

* **Variações de linha de base** ( *baseline* ): um "fundo" que eleva ou inclina todo o espectro, sem relação com a química da amostra.
* **Efeitos multiplicativos** : toda a intensidade do espectro escala para cima ou para baixo dependendo, por exemplo, da espessura da amostra ou da potência do laser.
* **Interferências químicas** : sinais de substâncias que não são de interesse — vapor d'água, CO₂ do ambiente, fluorescência.
* **Variação entre réplicas** : mesmo medindo a mesma amostra duas vezes, há diferenças instrumentais.

O **pré-processamento** é o conjunto de técnicas que nos permite *limpar* o espectro, separando o sinal de interesse de todos esses efeitos indesejados.

!!! tip "Dois tipos de pré-processamento"

- **Métodos de filtragem** : simplesmente transformam o espectro (ex.: derivada segunda, normalização vetorial). São rápidos, mas jogam fora informação.
- **Métodos baseados em modelo** : *quantificam* os efeitos indesejados e os removem de forma controlada. A informação filtrada não é perdida — fica armazenada nos parâmetros do modelo.

A **EMSC (Extended Multiplicative Signal Correction)** pertence à segunda categoria e é o foco deste texto.

---

## 2. Ponto de partida: a Lei de Lambert-Beer

Antes de entender a EMSC, precisamos relembrar a  **Lei de Lambert-Beer** , que é a base da espectroscopia de absorbância.

Para uma amostra transparente com uma única espécie química absorvedora, a absorbância $A(\tilde{\nu})$ em um número de onda $\tilde{\nu}$ é:

$$
A(\tilde{\nu}) = k(\tilde{\nu}) \cdot c \cdot b
$$

onde:

| Símbolo           | Significado                                           |
| ------------------ | ----------------------------------------------------- |
| $k(\tilde{\nu})$ | Absortividade (característica da substância)        |
| $c$              | Concentração da espécie química                   |
| $b$              | Comprimento do caminho óptico (espessura da amostra) |

!!! note "Intuição física"
Quanto mais espessa for a amostra ($b$ maior) ou mais concentrada ($c$ maior), maior será a absorbância. Se $b$ varia de medição para medição, todo o espectro escala proporcionalmente — isso é o chamado  **efeito multiplicativo** .

Para amostras biológicas complexas, com $J$ espécies absorvedoras ao mesmo tempo, a lei se generaliza como uma  **soma de contribuições** :

$$
A(\tilde{\nu}) = \left( \sum_{j=1}^{J} c_j \cdot k_j(\tilde{\nu}) \right) \cdot b
$$

---

## 3. Da Lei de Lambert-Beer ao modelo MSC

### 3.1 A ideia do espectro médio como referência

Em conjuntos de espectros biológicos, há uma observação prática muito importante:  **os espectros de amostras da mesma classe se parecem muito entre si** . Isso significa que o espectro médio $\bar{x}(\tilde{\nu})$ do conjunto é uma boa aproximação para cada espectro individual.

Matematicamente, podemos escrever o espectro de cada constituinte como:

$$
k_j(\tilde{\nu}) = \bar{x}(\tilde{\nu}) + \Delta k_j(\tilde{\nu})
$$

onde $\Delta k_j(\tilde{\nu})$ é o **desvio** daquele constituinte em relação à média. Substituindo na equação anterior e aplicando a condição $\sum_j c_j = 1$ (as concentrações normalizam para 100%), chegamos a:

$$
A(\tilde{\nu}) \approx \bar{x}(\tilde{\nu}) \cdot b + e(\tilde{\nu})
$$

onde $e(\tilde{\nu})$ é um resíduo que captura as diferenças químicas entre as amostras.

### 3.2 O modelo MSC básico

O modelo **MSC (Multiplicative Signal Correction)** adiciona a esse resultado um termo de linha de base constante $a$:

$$
\boxed{A(\tilde{\nu}) = a + \bar{x}(\tilde{\nu}) \cdot b + e(\tilde{\nu})}
$$

Os parâmetros $a$ e $b$ são estimados por **regressão de mínimos quadrados** (ordinary least squares — OLS). Uma vez estimados, o espectro corrigido é calculado como:

$$
A_{\text{corr}}(\tilde{\nu}) = \frac{A(\tilde{\nu}) - a}{b}
$$

!!! info "O que a correção faz?"

- Subtraímos $a$ para remover o deslocamento de linha de base constante.
- Dividimos por $b$ para normalizar o espectro, removendo o efeito multiplicativo (variação de espessura, por exemplo).
- O resíduo $e(\tilde{\nu})$ contém as **diferenças químicas reais** entre as amostras — exatamente o que queremos estudar!

---

## 4. Estendendo o modelo: a EMSC

O modelo MSC funciona bem quando a linha de base é aproximadamente constante. Mas e quando ela tem uma curvatura? Isso é muito comum em espectros  **Raman** , onde o fenômeno de **fluorescência** pode criar linhas de base fortemente curvas.

### 4.1 Extensão polinomial

A **EMSC** resolve esse problema adicionando termos polinomiais ao modelo:

$$
\boxed{A(\tilde{\nu}) = a + \bar{x}(\tilde{\nu}) \cdot b + d_1\tilde{\nu} + d_2\tilde{\nu}^2 + \cdots + d_n\tilde{\nu}^n + e(\tilde{\nu})}
$$

E a correção correspondente remove todos esses termos:

$$
A_{\text{corr}}(\tilde{\nu}) = \frac{A(\tilde{\nu}) - a - d_1\tilde{\nu} - d_2\tilde{\nu}^2 - \cdots - d_n\tilde{\nu}^n}{b}
$$

!!! tip "Qual grau do polinômio usar?"
Na prática, para espectros Raman com fluorescência, polinômios de **grau 6 ou 7** costumam ser suficientes. Polinômios de baixo grau descrevem apenas variações lentas e amplas, enquanto os picos químicos Raman são muito estreitos — portanto não há risco de confundir sinal químico com linha de base ao usar polinômios nessa faixa.

**Por que isso funciona sem superajuste ( *overfitting* )?** Porque os termos polinomiais descrevem variações *amplas e suaves* no espectro, enquanto os picos Raman de interesse são  *estreitos e agudos* . Esses dois tipos de estrutura são matematicamente muito diferentes (linearmente independentes), então o modelo os separa naturalmente.

---

## 5. Incorporando espectros constituintes

Até agora, as diferenças químicas entre amostras ficavam apenas no resíduo $e(\tilde{\nu})$. Em alguns casos, conhecemos *a priori* o espectro de uma interferência — e podemos incluí-la diretamente no modelo.

### 5.1 O modelo EMSC com espectros de diferença

Se soubermos o espectro de diferença $\Delta k_j(\tilde{\nu})$ de um constituinte interferente (ex.: vapor d'água), podemos incluí-lo explicitamente:

$$
A(\tilde{\nu}) = b \cdot \bar{x}(\tilde{\nu}) + \sum_{j=1}^{I} h_j \cdot \Delta k_j(\tilde{\nu}) + e(\tilde{\nu})
$$

onde $h_j = b \cdot c_j$ é o coeficiente que quantifica **o quanto** daquele constituinte está presente naquele espectro.

!!! example "Exemplo prático: água em tecido muscular"
Imagine espectros de FT-IR de seções de tecido de carne coletados em dias diferentes. A quantidade de **água ligada** no tecido varia com a umidade do ar no laboratório, afetando fortemente a região de estiramento O-H (3500–3400 cm⁻¹) e a região Amida I (~1640 cm⁻¹).

```
**Como corrigir isso com EMSC?**

1. Mede-se o espectro da mesma seção de tecido sob **alta umidade** e **baixa umidade**.
2. Calcula-se o **espectro de diferença**: alta umidade − baixa umidade. Esse é o $\Delta k_{\text{água}}(\tilde{\nu})$.
3. Inclui-se esse espectro de diferença no modelo EMSC.
4. O modelo estima automaticamente $h_j$ para cada espectro e remove a contribuição da água.

Resultado: a variação dia-a-dia desaparece completamente!
```

---

## 6. Usando subespaços ortogonais: a correção por réplicas

E se não soubermos exatamente qual é a interferência, mas tivermos **réplicas técnicas** das amostras? A EMSC permite usar essa informação.

### 6.1 A lógica da abordagem

Se medimos a mesma amostra química várias vezes (réplicas), qualquer diferença entre as réplicas é, por definição, **artefato físico** — não pode ser informação química. Podemos então:

1. Aplicar o modelo EMSC básico a todas as réplicas de cada amostra.
2. Centralizar na média ( *mean-center* ) as réplicas corrigidas de cada amostra.
3. Reunir todos os espectros de variação entre réplicas em uma única tabela.
4. Aplicar **PCA (Análise de Componentes Principais)** para extrair os padrões dominantes dessa variação.

Os vetores de loadings $p_k(\tilde{\nu})$ obtidos representam os **principais modos de variação física** no equipamento.

### 6.2 O modelo EMSC com correção por réplicas

Esses componentes são então adicionados ao modelo global:

$$
A(\tilde{\nu}) = a + \bar{x}(\tilde{\nu}) \cdot b + d_1\tilde{\nu} + d_2\tilde{\nu}^2 + \sum_{k=1}^{K} g_k \cdot p_k(\tilde{\nu}) + e(\tilde{\nu})
$$

E a correção remove todos os termos de interferência:

$$
A_{\text{corr}}(\tilde{\nu}) = \frac{A(\tilde{\nu}) - a - d_1\tilde{\nu} - d_2\tilde{\nu}^2 - \sum_{k=1}^{K} g_k \cdot p_k(\tilde{\nu})}{b}
$$

!!! warning "Atenção: seja restritivo no número de componentes!"
Como os modos de variação física são  **estimados a partir dos dados** , incluir componentes demais pode remover informação química real. É necessário inspecionar visualmente os loadings e usar validação para decidir quantos componentes incluir. Em geral, os loadings associados a artefatos físicos têm estrutura  **larga e suave** , enquanto sinais químicos têm picos  **estreitos e agudos** .

---

## 7. Resumo: a hierarquia dos modelos EMSC

A tabela abaixo resume a evolução dos modelos, do mais simples ao mais completo:

| Modelo                         | Equação                                                           | Quando usar                                                |
| ------------------------------ | ------------------------------------------------------------------- | ---------------------------------------------------------- |
| **MSC**                  | $A = a + \bar{x} \cdot b + e$                                     | Linhas de base constantes, efeitos multiplicativos simples |
| **EMSC básico**         | $A = a + \bar{x} \cdot b + d_1\tilde{\nu} + d_2\tilde{\nu}^2 + e$ | Linhas de base com curvatura moderada                      |
| **EMSC polinomial**      | $A = a + \bar{x} \cdot b + \sum d_n\tilde{\nu}^n + e$             | Fluorescência forte em Raman (grau 6–7)                  |
| **EMSC + constituintes** | $A = b\bar{x} + \sum h_j\Delta k_j + e$                           | Interferência química conhecida (ex.: vapor d'água)     |
| **EMSC + subespaços**   | $A = a + \bar{x} \cdot b + \sum g_k p_k + e$                      | Variação física entre réplicas desconhecida            |

---

## 8. O significado dos parâmetros estimados

Uma das grandes vantagens da EMSC sobre métodos de filtragem é que ela **não apenas remove** as interferências — ela as  **quantifica** .

Cada parâmetro estimado pelo modelo tem um significado físico ou químico:

* **$b$** (parâmetro multiplicativo): proporcional à espessura efetiva da amostra ou ao volume de amostragem (em Raman).
* **$a$, $d_1$, $d_2$, ...** (parâmetros de linha de base): descrevem o comportamento do fundo espectral.
* **$h_j$** (coeficientes de constituintes): quantificam a contribuição de cada interferência química conhecida.
* **$g_k$** (coeficientes de subespaço): medem o quanto cada modo de variação física afetou aquela medição específica.

Esses parâmetros podem ser usados como **variáveis auxiliares** em análises posteriores. Por exemplo, Kohler et al. (2005) relacionaram os parâmetros EMSC de espectros FT-IR de seções de músculo diretamente a propriedades texturais do músculo após tratamento térmico.

---

## 9. Considerações práticas

### Rank do modelo

Para que a estimativa por mínimos quadrados funcione corretamente, a **matriz de espectros do modelo** deve ter posto ( *rank* ) completo. Isso significa que os espectros incluídos no modelo (referência, polinômios, espectros de diferença, loadings) não podem ser linearmente dependentes entre si. Na prática:

* Os polinômios são sempre independentes entre si e da referência.
* Espectros de diferença $\Delta k_j$ são independentes por definição.
* Os loadings de PCA são ortogonais por construção.

### Escolha do espectro de referência

Usualmente, o espectro médio $\bar{x}(\tilde{\nu})$ do conjunto de dados é usado como referência. Uma referência diferente apenas desloca o nível do espectro corrigido, mas **não afeta** a análise multivariada subsequente (PCA, PLS, etc.).


---

## 10. Referências

> Afseth, N.K., Kohler, A. (2012).  **Extended multiplicative signal correction in vibrational spectroscopy, a tutorial** .  *Chemometrics and Intelligent Laboratory Systems* , 117, 92–99. https://doi.org/10.1016/j.chemolab.2012.03.004

> Martens, H., Stark, E. (1991). Extended multiplicative signal correction and spectral interference subtraction — new preprocessing methods for near-infrared spectroscopy.  *Journal of Pharmaceutical and Biomedical Analysis* , 9, 625–635.

> Kohler, A., et al. (2005). Extended multiplicative signal correction as a tool for separation and characterization of physical and chemical information in FT-IR microscopy images.  *Applied Spectroscopy* , 59, 707–716.
>
