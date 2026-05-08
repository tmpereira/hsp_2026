# Leitura de Dados

A biblioteca HSP suporta nativamente os principais formatos de arquivo
gerados por microscópios de FTIR comerciais, além de um formato interno
baseado em NumPy para armazenamento eficiente.

---

## Perkin Elmer Spotlight — `.fsm`

```python
import hsp_2026.file as file

# Carrega um arquivo .fsm
data = file.fsm('amostra.fsm')

print(data['r'].shape)    # (n_pixels, n_pontos_espectrais)
print(data['wn'][[0,-1]]) # primeiro e último número de onda
print(data['dx'], data['dy'])  # dimensões da imagem
```

Para listar todos os arquivos `.fsm` em um diretório:

```python
arquivos = file.get_fsm_files('/caminho/para/pasta/')
for arq in arquivos:
    data = file.fsm(arq)
```

### Conversão de transmitância para absorbância

O leitor `.fsm` converte automaticamente os dados de **transmitância** (%)
para **absorbância** usando a Lei de Beer-Lambert:

$$
A = -\log_{10}\!\left(\frac{T}{100}\right)
$$

---

## Agilent Cary 620/670 — `.dmd` + `.dmt`

O formato Agilent armazena a imagem em mosaico de tiles de 32 × 32 pixels.
O leitor reconstrói automaticamente a grade completa:

```python
data = file.age('/caminho/para/pasta_com_tiles/')
```

!!! note
    A pasta deve conter os pares de arquivos `.dmd` (dados espectrais) e
    `.dmt` (metadados), organizados pelo software Agilent MicroLab.

---

## Formato interno NumPy — `.npz`

Para evitar re-processar arquivos pesados a cada sessão, salve o dicionário
HSP já processado:

```python
# Salva após pré-processamento
file.npz_save('amostra_processada', data)

# Recarrega numa sessão futura
data = file.npz_load('amostra_processada.npz')
```
