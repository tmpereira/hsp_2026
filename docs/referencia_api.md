# Referência da API

Documentação completa de todas as funções públicas da biblioteca HSP.
Os docstrings são renderizados automaticamente a partir do código-fonte
pelo plugin `mkdocstrings`.

---

## `hsp_2026.file` — Leitura de Arquivos

::: hsp_2026.file
    options:
      members:
        - fsm
        - get_fsm_files
        - age
        - npz_save
        - npz_load

---

## `hsp_2026.prep` — Pré-Processamento

::: hsp_2026.prep
    options:
      members:
        - cut
        - snv
        - norm
        - golay
        - norm2r
        - pcares
        - napc
        - offset
        - binned
        - rand
        - dsample

---

## `hsp_2026.emsc` — Correção de Espalhamento

::: hsp_2026.emsc
    options:
      members:
        - pca
        - create_model
        - create_model_h2o
        - emsc_model_view
        - emsc_fit

---

## `hsp_2026.qt` — Controle de Qualidade

::: hsp_2026.qt
    options:
      members:
        - area
        - intt
        - emsc
        - mean
        - otsu_area
        - otsu_emsc

---

## `hsp_2026.rg` — Histogramas

::: hsp_2026.rg
    options:
      members:
        - area
        - intt
        - emsc

---

## `hsp_2026.sh` — Visualização de Imagens

::: hsp_2026.sh
    options:
      members:
        - intt
        - area
        - mean
        - pplot
        - emsc
        - pc
        - int_plt

---

## `hsp_2026.km` — Segmentação K-Means

::: hsp_2026.km
    options:
      members:
        - fit
        - sh
        - spc
        - fit_common
        - sh2

---

## Estrutura do dicionário HSP

Todas as funções da biblioteca utilizam e retornam um dicionário padrão com
as seguintes chaves:

| Chave      | Tipo                     | Descrição                                          |
|------------|--------------------------|----------------------------------------------------|
| `r`        | `ndarray float32 (n, p)` | Absorbâncias: n espectros × p pontos espectrais    |
| `wn`       | `ndarray float64 (p,)`   | Números de onda em cm⁻¹                            |
| `dx`       | `int`                    | Número de linhas da imagem                         |
| `dy`       | `int`                    | Número de colunas da imagem                        |
| `sel`      | `ndarray bool (n,)`      | Máscara de pixels válidos (`True` = válido)        |
| `filename` | `str`                    | Nome do arquivo de origem                          |
| `log`      | `str`                    | Histórico das operações aplicadas                  |

Chaves adicionais criadas pelo módulo `emsc`:

| Chave          | Tipo                       | Descrição                              |
|----------------|----------------------------|----------------------------------------|
| `EMSC_model`   | `ndarray (p, m)`           | Matriz do modelo EMSC                  |
| `EMSC_coeff`   | `ndarray (n_total, m)`     | Coeficientes do ajuste por pixel       |
