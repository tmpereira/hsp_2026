> Módulo `rg.py` — Histogramas para Definição de Limiares

O módulo `rg` (*range*) é um **auxiliar de diagnóstico visual**.
Ele gera histogramas das métricas usadas pelos filtros do `qt.py`,
permitindo que você escolha os limiares `a` e `b` de forma informada
**antes** de aplicar qualquer remoção de pixels.

---

## Por que inspecionar antes de filtrar?

Os filtros do `qt.py` pedem parâmetros `a` (mínimo) e `b` (máximo):

```python
data = qt.area(data, 900, 1800, a=50, b=500)
```

Se você escolher `a` e `b` sem olhar os dados, corre o risco de:

- **`a` muito alto** → remover pixels válidos de amostra (sub-segmentação)
- **`a` muito baixo** → deixar pixels de fundo passarem (sobre-segmentação)
- **`b` muito baixo** → eliminar pixels de amostras espessas que são legítimos

O histograma exibe a **distribuição real** da métrica para os seus dados
e torna a escolha dos limiares visual e intuitiva.

### Anatomia típica de um histograma de área

```
Frequência
    │
    │   ▐██▌                              ← pico 1: fundo (substrato, bolhas)
    │   ▐██▌▌                                área muito baixa → remover
    │   ▐██████▌
    │         ▐████████████▌             ← pico 2: amostra real
    │              ▐████████████████▌       amplitude normal → manter
    │                          ▐▌
    └─────────────────────────────────── Área integrada
          ↑ a (mínimo)        ↑ b (máximo)
```

O limiar `a` deve ser colocado no **vale entre os dois picos**.
O limiar `b` deve cobrir a cauda direita da amostra real, excluindo
apenas eventuais outliers de saturação.

---

## Workflow completo com `rg.py` + `qt.py`

O ciclo de uso é sempre o mesmo, independente da métrica:

```
rg.area()  →  inspeciona histograma  →  qt.area(a=…, b=…)
rg.intt()  →  inspeciona histograma  →  qt.intt(a=…, b=…)
rg.emsc()  →  inspeciona histograma  →  qt.emsc(a=…, b=…)
```

Exemplo completo:

```python
import hsp_2026.file as file
import hsp_2026.rg   as rg
import hsp_2026.qt   as qt

data = file.fsm('tecido.fsm')

# ── Etapa 1: diagnóstico visual ──────────────────────────────────────────
rg.area(data, 900, 1800)
# → Identifica visualmente: fundo abaixo de 40, amostra entre 40 e 450

# ── Etapa 2: aplica o filtro com os valores identificados ────────────────
data = qt.area(data, 900, 1800, a=40, b=450)
```

Quando os dois grupos se sobrepõem muito no histograma (amostras
heterogêneas ou muito espessas), prefira `qt.otsu_area()` ou
`qt.otsu_emsc()`, que determinam o limiar automaticamente.

---

## Funções do módulo

### `area(data, ini, fim)` — Histograma das áreas integradas

```python
rg.area(data, ini=900, fim=1800)
```

**O que faz:**

Para cada pixel, calcula a **área integrada** do espectro na faixa
`[ini, fim]` cm⁻¹ usando a regra do trapézio:

$$
A_i = \int_{\text{ini}}^{\text{fim}} z_i(\tilde{\nu})\,d\tilde{\nu}
\;\approx\;
\sum_{k} \frac{z_i(\tilde{\nu}_k) + z_i(\tilde{\nu}_{k+1})}{2}
\,\Delta\tilde{\nu}
$$

onde \(z_i\) é o espectro do pixel \(i\) e \(\tilde{\nu}\) é o número
de onda. O histograma com 300 bins é gerado para todos os \(A_i\).

**Parâmetros:**

| Parâmetro | Tipo  | Descrição                              |
| ---------- | ----- | ---------------------------------------- |
| `data`   | dict  | Dicionário hsp com `'r'` e `'wn'`   |
| `ini`    | float | Limite inferior da integração (cm⁻¹) |
| `fim`    | float | Limite superior da integração (cm⁻¹) |

**Qual região usar?**

A escolha de `ini` e `fim` determina quais bandas entram na integral.
Algumas opções comuns:

| Região              | ini  | fim  | Conteúdo principal         |
| -------------------- | ---- | ---- | --------------------------- |
| Fingerprint completa | 900  | 1800 | Todas as bandas biológicas |
| Região Amida        | 1550 | 1700 | Proteínas (Amida I e II)   |
| Lipídios C–H       | 2800 | 3000 | Ácidos graxos              |
| Região Amida        | 900  | 1800 | Análise geral de tecido    |

A região fingerprint completa (900–1800 cm⁻¹) é a escolha mais segura
para separar fundo de amostra biológica porque inclui múltiplas bandas.

**Exemplo de saída esperada:**

```
Título do gráfico: "histograma de área entre 900 até 1800"
Eixo x: Área integrada
Eixo y: Frequência (contagem de pixels)
```

---

### `intt(data, b)` — Histograma da intensidade num ponto espectral

```python
rg.intt(data, b=1650)
```

**O que faz:**

Extrai a **intensidade de absorção no número de onda mais próximo de `b`**
para cada pixel e plota o histograma dessas intensidades.

Internamente, a função seleciona `data['wn'] > b` e pega a primeira
coluna resultante — ou seja, o ponto imediatamente acima de `b` no vetor
de números de onda.

**Parâmetros:**

| Parâmetro | Tipo  | Descrição                            |
| ---------- | ----- | -------------------------------------- |
| `data`   | dict  | Dicionário hsp com `'r'` e `'wn'` |
| `b`      | float | Número de onda de interesse (cm⁻¹)  |

**Quando usar no lugar de `rg.area()`?**

`rg.intt()` é útil quando você quer avaliar a intensidade de uma **banda
muito específica e característica**. Exemplos práticos:

| Banda             | b (cm⁻¹) | Molécula  | Aplicação                    |
| ----------------- | ---------- | ---------- | ------------------------------ |
| Amida I           | 1650       | Proteínas | Filtrar pixels sem proteína   |
| Amida II          | 1546       | Proteínas | Complementar à Amida I        |
| C=O éster        | 1740       | Lipídios  | Filtrar em estudos de membrana |
| CH₂ assimétrico | 2924       | Lipídios  | Lipídios insaturados          |

**Limitação:**

Como a intensidade num único ponto é mais sensível a ruído do que a área
integrada numa região, `rg.intt()` pode mostrar histogramas mais
dispersos. Se os dois grupos (fundo e amostra) não aparecerem bem
separados, use `rg.area()` com uma janela estreita em torno da banda
de interesse.

**⚠️ Detalhe de implementação:**

A função busca `data['wn'] > b`, portanto usa o primeiro ponto **acima**
de `b`, não o ponto mais próximo. Para um vetor de números de onda
espaçado de 4 cm⁻¹, a diferença é pequena; mas esteja ciente disso ao
interpretar o título do histograma.

---

### `emsc(data, a)` — Histograma de um coeficiente do modelo EMSC

```python
# Visualiza o coeficiente a₀ (escala da referência)
rg.emsc(data, a=0)

# Visualiza o coeficiente de baseline linear (índice 2)
rg.emsc(data, a=2)
```

**O que faz:**

Plota o histograma do coeficiente de índice `a` em `data['EMSC_coeff']`
para todos os pixels. Requer que `emsc.emsc_fit()` tenha sido executado
antes.

**Parâmetros:**

| Parâmetro | Tipo | Descrição                          |
| ---------- | ---- | ------------------------------------ |
| `data`   | dict | Dicionário hsp com `'EMSC_coeff'` |
| `a`      | int  | Índice do coeficiente (0-based)     |

**Entendendo os índices dos coeficientes:**

O modelo EMSC ajusta cada espectro como combinação linear de uma
referência, polinômios de baseline e espectros de interferentes:

$$
z = a_0 \cdot \mathbf{r}_{\text{ref}} + b_0 + b_1\lambda + b_2\lambda^2 + \ldots + \sum_k c_k\mathbf{p}_k + \varepsilon
$$

O vetor `data['EMSC_coeff']` armazena os coeficientes ajustados para
cada pixel na seguinte ordem:

| Índice `a` | Coeficiente | Interpretação física                                           |
| ------------- | ----------- | ----------------------------------------------------------------- |
| 0             | \(a_0\)     | Escala da referência — quanto de "sinal biológico" tem o pixel |
| 1             | \(b_0\)     | Offset (baseline constante)                                       |
| 2             | \(b_1\)     | Inclinação linear da baseline                                   |
| 3             | \(b_2\)     | Curvatura da baseline                                             |
| 4, 5, …      | \(c_k\)     | Contribuição de cada interferente (parafina, água, etc.)       |

**Por que o coeficiente \(a_0\) é especialmente útil?**

O \(a_0\) representa diretamente **o quanto aquele pixel "se parece"
com o espectro de referência da amostra**, já descontados os efeitos
de baseline e interferentes. Por isso:

- Pixels de fundo (substrato, bolhas): \(a_0 \approx 0\) ou negativo
- Pixels de amostra real: \(a_0 > 0\), com valores típicos entre 0.5 e 2.5

O histograma do \(a_0\) tende a mostrar uma separação muito mais nítida
entre fundo e amostra do que o histograma da área bruta, especialmente
em amostras embebidas em parafina ou hidratadas.

**Exemplo de uso completo:**

```python
import hsp_2026.emsc as emsc
import hsp_2026.rg   as rg
import hsp_2026.qt   as qt

# 1. Constrói e aplica o modelo EMSC
model = emsc.create_model(ref_espectro, parafina_loadings, n_poly=3, n_pc=2)
data  = emsc.emsc_fit(data, model)

# 2. Inspeciona a distribuição do coeficiente a₀
rg.emsc(data, a=0)
# → Dois grupos claros: fundo (a₀ < 0.2) e amostra (a₀ > 0.2)

# 3. Aplica o filtro com os limites identificados
data = qt.emsc(data, ini=0, a=0.2, b=3.0)
```

---

## Tabela comparativa das funções

| Função                    | Usa `emsc_fit()`? | Métrica                    | Companion em `qt.py` |
| --------------------------- | ------------------- | --------------------------- | ---------------------- |
| `rg.area(data, ini, fim)` | Não                | Área integrada (trapézio) | `qt.area()`          |
| `rg.intt(data, b)`        | Não                | Intensidade num ponto       | `qt.intt()`          |
| `rg.emsc(data, a)`        | **Sim**       | Coeficiente do modelo EMSC  | `qt.emsc()`          |

---

## Posição no pipeline geral

O `rg.py` é exclusivamente uma ferramenta de **diagnóstico e exploração**.
Ele não modifica `data` e deve ser chamado antes dos filtros do `qt.py`:

```
file.fsm()          → carrega dados brutos
     │
     ├─ rg.area()   → histograma → escolha de a, b
     │
qt.area()           → remove pixels defeituosos (modifica data)
     │
prep.cut()
prep.snv()
     │
sh.pc() / km.fit()  → análise multivariada
```

Use `rg.py` **quantas vezes precisar** — ele não consome os dados, apenas
plota. É comum rodar `rg.area()`, ajustar `a` e `b`, aplicar `qt.area()`,
e depois rodar `rg.area()` novamente para confirmar que o histograma
resultante está adequado.
