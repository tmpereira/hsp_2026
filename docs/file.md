# Módulo `file.py` — Leitura de Arquivos de Micro-FTIR

Este módulo é responsável por **abrir e converter** os arquivos gerados
pelos microscópios de FTIR para o formato interno da biblioteca (o
"dicionário hsp"). É o ponto de entrada de qualquer análise: sem este
módulo, não há dados para processar.

---

## O que é um arquivo de Micro-FTIR?

Quando você faz uma medição com um microscópio de FTIR, o equipamento gera
um **cubo de dados**: imagine uma pilha de imagens, em que cada imagem
corresponde a um comprimento de onda diferente. Outro jeito de visualizar
é pensar que cada pixel da imagem contém um **espectro completo** de
infravermelho.

```
     Imagem (espaço)
     ┌──────────────┐
     │  pixel (1,1) │──► [ absorção em 900 cm⁻¹, 901 cm⁻¹, ..., 1800 cm⁻¹ ]
     │  pixel (1,2) │──► [ absorção em 900 cm⁻¹, 901 cm⁻¹, ..., 1800 cm⁻¹ ]
     │  ...         │
     └──────────────┘
```

Cada fabricante salva esse cubo em um formato binário proprietário.
O módulo `file.py` sabe como decodificar os dois formatos mais comuns nos
laboratórios brasileiros:

| Fabricante | Arquivo | Função |
|------------|---------|--------|
| Perkin Elmer Spotlight | `.fsm` | `fsm()` |
| Agilent Cary 620/670   | `.dmd` + `.dmt` | `age()` |

---

## O dicionário hsp — a moeda da biblioteca

Todas as funções de leitura retornam um **dicionário Python** com sempre
as mesmas chaves. Esse padrão garante que qualquer função de processamento
(`prep.py`, `emsc.py`, etc.) saiba exatamente onde encontrar os dados.

| Chave | Tipo | O que armazena |
|-------|------|----------------|
| `r` | `ndarray float32` (n_pixels × n_pontos) | Absorbâncias de todos os pixels |
| `wn` | `ndarray float64` (n_pontos,) | Números de onda em cm⁻¹ |
| `dx` | `int` | Número de **linhas** da imagem |
| `dy` | `int` | Número de **colunas** da imagem |
| `sel` | `ndarray bool` (n_pixels,) | Máscara: `True` = pixel válido |
| `filename` | `str` | Nome do arquivo de origem |
| `log` | `str` | Histórico de operações aplicadas |

> **Dica para iniciantes:** o número de linhas em `r` é sempre
> `dx × dy` (total de pixels). Para reconstruir a imagem 2D, basta
> fazer `r.reshape(dx, dy, -1)`.

---

## Funções de uso direto (API pública)

### `fsm(arq)` — Ler arquivo Perkin Elmer Spotlight

```python
data = file.fsm('minha_amostra.fsm')
```

**O que faz passo a passo:**

1. Abre o arquivo `.fsm` em modo binário e valida se os primeiros 4 bytes
   são `PEPE` (a "assinatura" que identifica o formato Perkin Elmer).
   Se não for, lança um erro explicando o problema.

2. Percorre todos os **blocos** de dados do arquivo. Um bloco é como um
   capítulo de um livro: cada um tem um número de identificação (ID) e
   um conteúdo. Os blocos relevantes são:
   - **Bloco 5100:** dimensões da imagem, número de pontos espectrais,
     intervalo de números de onda
   - **Bloco 5104:** metadados do instrumento (analista, data, modelo)
   - **Bloco 5105:** os espectros brutos em transmitância

3. **Converte transmitância → absorbância** usando a Lei de Beer-Lambert:

$$
A = -\log_{10}\!\left(\frac{T}{100}\right)
$$

   Por quê isso é necessário? O equipamento mede a **transmitância** (fração
   da luz que passou pela amostra), mas para análise química usamos
   **absorbância**, que é proporcional à concentração do composto.

4. **Inverte o eixo espectral** com `np.fliplr` e `np.flipud`, pois o
   arquivo FSM armazena os espectros do maior para o menor número de onda,
   e a convenção da biblioteca é do menor para o maior.

5. Retorna o dicionário hsp com todos os pixels marcados como válidos
   (`sel = True` para todos).

**Exemplo de uso:**

```python
import hsp_2026.file as file
import hsp_2026.sh as sh

data = file.fsm('tecido_higido.fsm')
print(data['r'].shape)     # ex.: (3600, 426) → 3600 pixels, 426 pontos
print(data['wn'][[0, -1]]) # ex.: [899.5, 1799.8] → intervalo espectral
print(data['dx'], data['dy'])  # ex.: 60 60 → imagem 60×60 pixels

# Visualiza a intensidade em 1650 cm⁻¹ (banda de proteínas - Amida I)
sh.intt(data, 1650)
```

---

### `get_fsm_files(path)` — Listar arquivos FSM em um diretório

```python
lista = file.get_fsm_files('/dados/experimento_01/')
```

**O que faz:**

Varre o diretório `path` e retorna uma lista com o caminho completo de
todos os arquivos `.fsm` encontrados. Útil quando você tem vários arquivos
de uma série experimental e quer processá-los em loop.

> ⚠️ **Atenção:** esta função usa `os.chdir(path)` internamente, o que
> muda o diretório de trabalho atual do Python. Se você estiver rodando
> um notebook Jupyter, saiba que isso pode afetar onde outros arquivos são
> salvos. Prefira usar `os.chdir` explicitamente ou trabalhar com caminhos
> absolutos quando possível.

**Exemplo de uso:**

```python
import hsp_2026.file as file
import hsp_2026.prep as prep

arquivos = file.get_fsm_files('/dados/experimento_01/')
for arq in arquivos:
    data = file.fsm(arq)
    data = prep.cut(data, 900, 1800)  # corta a região fingerprint
    # ... demais processamentos
```

---

### `age(path)` — Ler mosaico Agilent (`.dmd` + `.dmt`)

```python
data = file.age('/dados/mosaico_tecido/')
```

**O que faz passo a passo:**

O microscópio Agilent divide a imagem em **pedaços menores chamados tiles**,
cada um com 32 × 32 pixels. Cada tile é salvo como um arquivo `.dmd`
separado. O arquivo `.dmt` contém os metadados (número de pontos espectrais,
número de onda inicial e o passo espectral).

1. **Descobre as dimensões do mosaico** lendo os nomes dos arquivos `.dmd`.
   O nome segue o padrão `XXXX_YYYY.dmd`, onde XXXX e YYYY são os índices
   do tile nas direções x e y.

2. **Lê os metadados espectrais** do arquivo `.dmt`:
   - `npoints`: quantos pontos espectrais existem por espectro
   - `start`: o número de onda inicial (em unidades internas)
   - `steps`: o passo espectral em cm⁻¹
   - Reconstrói o vetor `wn` com `np.linspace`

3. **Lê cada tile `.dmd`**, reorienta os eixos com `np.transpose` e
   `[:, ::-1, :]` (inversão necessária pela convenção do Agilent), e
   **encaixa cada tile na posição correta** do array global.

4. Achata o array 3D `(linhas, colunas, n_pontos)` para 2D
   `(n_pixels, n_pontos)` e retorna o dicionário hsp.

**Exemplo de uso:**

```python
data = file.age('/dados/agilent/amostra_001/')
print(data['dx'], data['dy'])   # ex.: 160 128 → imagem de 160×128 pixels
```

---

### `npz_save(arq, data)` — Salvar em formato NumPy

```python
file.npz_save('amostra_processada', data)
# cria o arquivo: amostra_processada.npz
```

**O que faz:**

Salva o dicionário hsp inteiro (com todos os seus arrays e metadados) em
um único arquivo `.npz`, que é o formato de arquivo comprimido do NumPy.

**Por que usar?**

- Abrir um arquivo `.fsm` grande pode levar vários segundos (parsing binário)
- Após pré-processar os dados (corte, SNV, EMSC), você provavelmente quer
  salvar o estado processado
- Recarregar um `.npz` é muito mais rápido do que reprocessar do zero

---

### `npz_load(arq)` — Carregar formato NumPy

```python
data = file.npz_load('amostra_processada.npz')
```

**O que faz:**

Carrega um arquivo `.npz` previamente salvo por `npz_save` e reconstrói
o dicionário hsp com todas as chaves e arrays originais.

O parâmetro `allow_pickle=True` é necessário porque o dicionário contém
strings (`filename`, `log`) além de arrays numéricos, e o NumPy por
segurança desabilita o carregamento de objetos Python arbitrários por
padrão.

**Fluxo de trabalho recomendado:**

```python
# ── Sessão 1: processamento pesado ──────────────────────
data = file.fsm('tecido.fsm')          # abre arquivo binário
data = prep.cut(data, 900, 1800)       # pré-processa
data = prep.snv(data)
emsc_model = emsc.create_model(...)
data = emsc.emsc_fit(data, emsc_model)
file.npz_save('tecido_processado', data)  # salva estado processado

# ── Sessão 2: análise rápida ─────────────────────────────
data = file.npz_load('tecido_processado.npz')  # carrega em < 1 segundo
sh.pc(data, n=1)   # já pode ir direto para a análise
```

---

## Funções internas (prefixo `_`)

Estas funções não são chamadas diretamente pelo usuário. Elas são usadas
internamente por `fsm()` para decodificar o arquivo binário FSM peça por peça.
O prefixo `_` é uma convenção Python que sinaliza: "esta função é interna,
não é parte da API pública".

> **Para fins didáticos**, entender estas funções ajuda a compreender como
> um arquivo binário proprietário é decodificado "do zero", sem depender
> de bibliotecas externas.

---

### `_fsm_block_info(data)` — Lê o cabeçalho de um bloco

```python
block_id, block_size = _fsm_block_info(content[pos:pos + 6])
```

Recebe 6 bytes e os interpreta como:
- **2 bytes** → `block_id` (número inteiro sem sinal): identifica o tipo do bloco
- **4 bytes** → `block_size` (número inteiro com sinal): tamanho do bloco em bytes

Usa `struct.unpack('<Hi', data)`:
- `<` = little-endian (o byte menos significativo vem primeiro, padrão Intel/Windows)
- `H` = unsigned short (2 bytes, inteiro sem sinal)
- `i` = signed int (4 bytes, inteiro com sinal)

---

### `_fsm_decode_5100(data)` — Decodifica os metadados do mapeamento

Este é o bloco mais importante do arquivo FSM. Ele contém:

| Campo | Significado |
|-------|-------------|
| `n_x` | Número de colunas da imagem (pixels no eixo horizontal) |
| `n_y` | Número de linhas da imagem (pixels no eixo vertical) |
| `n_z` | Número de pontos espectrais por pixel |
| `z_start` | Número de onda inicial em cm⁻¹ |
| `z_end` | Número de onda final em cm⁻¹ |
| `z_delta` | Espaçamento entre pontos espectrais em cm⁻¹ |
| `x_delta`, `y_delta` | Tamanho de cada pixel em µm |

A decodificação usa `struct.unpack` com o formato
`'<ddddddddddiiihBhBhBhB'`, onde cada letra representa um tipo de dado:
`d` = double (8 bytes), `i` = int (4 bytes), `h` = short (2 bytes),
`B` = unsigned char (1 byte).

---

### `_fsm_decode_5104(data)` — Decodifica as informações do instrumento

Este bloco contém texto com informações operacionais:
analista, data da medição, modelo e número de série do equipamento,
versão do software.

O formato é uma sequência de registros com **tags de 2 bytes**:
- `#u` → uma string UTF-8 (precedida do seu tamanho em 2 bytes)
- `$u` e `,u` → um número inteiro curto (2 bytes)

A função percorre o bloco byte a byte, identifica a tag, extrai o valor
correspondente e armazena em uma lista. No final, mapeia os valores
extraídos para nomes de campos conhecidos (`analyst`, `date`,
`instrument_model`, etc.).

---

### `_fsm_decode_5105(data)` — Decodifica os espectros brutos

Este é o maior bloco do arquivo: contém todos os valores de transmitância
de todos os pixels concatenados em sequência, como `float32` (ponto
flutuante de 32 bits, precisão simples).

A decodificação é feita com uma única linha:

```python
np.frombuffer(data, dtype=np.float32)
```

`np.frombuffer` interpreta os bytes como um array NumPy sem copiar os
dados na memória (operação eficiente). O reshape para
`(n_espectros, n_pontos)` é feito em `_read_fsm_binary`.

---

### `_read_fsm_binary(path)` — Orquestra a leitura completa do FSM

Esta é a função que coordena tudo: abre o arquivo, valida a assinatura
`PEPE`, e percorre os blocos em um loop, chamando o decoder adequado para
cada ID encontrado. Os blocos com IDs desconhecidos são simplesmente
ignorados.

O fluxo é:

```
Abre arquivo
    │
    ├── Valida assinatura 'PEPE'
    │
    └── Loop pelos blocos:
            │
            ├── Bloco 5100 ──► _fsm_decode_5100() ──► atualiza meta{}
            ├── Bloco 5104 ──► _fsm_decode_5104() ──► atualiza meta{}
            ├── Bloco 5105 ──► _fsm_decode_5105() ──► acumula espectros
            └── outros IDs ──► ignora

Monta amplitudes (ndarray 2D) e wavelength (vetor de cm⁻¹)
Retorna (amplitudes, wavelength, meta)
```

A separação entre `_read_fsm_binary` e `fsm` existe porque a primeira faz
apenas a leitura binária pura, e a segunda aplica as transformações
científicas (conversão de unidades, orientação dos eixos) e monta o
dicionário hsp.
