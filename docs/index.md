# HSP — Hyperspectral Spectral Processing

**HSP** é uma biblioteca Python de código aberto para aquisição, pré-processamento e análise
quantitativa de imagens hiperespectrais obtidas por **Micro-FTIR**
(Microscopia de Infravermelho por Transformada de Fourier).

---

## Motivação científica

A espectroscopia de infravermelho por transformada de Fourier no modo de
microimagem (Micro-FTIR) permite mapear a distribuição espacial de compostos
químicos em tecidos, células e materiais com resolução micrométrica, sem
necessidade de marcadores exógenos.  Cada pixel da imagem contém um espectro
completo de absorbância na região do IR médio (400 – 4000 cm⁻¹), gerando um
cubo hiperespectral:

$$
\mathbf{X} \in \mathbb{R}^{n_x \times n_y \times n_\lambda}
$$

onde \(n_x \times n_y\) é o número de pixels e \(n_\lambda\) é o número de
pontos espectrais.

O processamento desse cubo envolve etapas como remoção de variações físicas
de espalhamento (EMSC), normalização, derivação espectral e análise
multivariada (PCA, K-Means), que são implementadas e documentadas nesta
biblioteca.

---

## Visão geral dos módulos

| Módulo  | Função principal                                                           |
| -------- | ---------------------------------------------------------------------------- |
| `file` | Leitura de formatos binários FSM (Perkin Elmer) e Agilent DMD/DMT           |
| `prep` | Pré-processamento espectral: corte, normalização SNV, Savitzky-Golay, PCA |
| `emsc` | Correção de espalhamento multiplicativo estendido (EMSC)                   |
| `qt`   | Controle de qualidade: remoção de pixels defeituosos                       |
| `rg`   | Histogramas para inspeção das métricas de qualidade                       |
| `sh`   | Geração de imagens hiperespectrais (intensidade, área, PCA, EMSC)         |
| `km`   | Segmentação por K-Means com visualização de clusters                     |

---

## Início rápido

```python
import hsp_2026.file as file
import hsp_2026.prep as prep
import hsp_2026.sh   as sh

# 1. Carrega um arquivo Perkin Elmer Spotlight (.fsm)
data = file.fsm('minha_amostra.fsm')

# 2. Recorta a região fingerprint do IR biológico
data = prep.cut(data, 900, 1800)

# 3. Normalização SNV
data = prep.snv(data)

# 4. Imagem da intensidade na banda da Amida I (1650 cm⁻¹)
sh.intt(data, 1650)

# 5. Imagem da PC1 (análise de componentes principais)
scores, loadings, var = sh.pc(data, n=1)
```

---

## Instalação

```bash

```

```bash
git clone https://github.com/seu_usuario/hsp_2026
cd hsp_2026
pip install -e .
```

Consulte a página [Instalação](instalacao.md) para detalhes sobre dependências.

---


## Licença

MIT License — veja o arquivo `LICENSE` no repositório.

---

*Desenvolvido pelo Prof. Thiago Martini Pereira — Universidade Federal de São Paulo (UNIFESP)*
