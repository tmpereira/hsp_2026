# Instalação

## Requisitos do sistema

| Requisito | Versão mínima |
|-----------|---------------|
| Python    | 3.8           |
| NumPy     | 1.21          |
| Matplotlib| 3.5           |
| SciPy     | 1.7           |
| scikit-learn | 1.0        |
| scikit-image | 0.19       |
| pandas    | 1.4           |

---

## Instalação via pip

```bash
pip install hsp_2026
```

---

## Instalação para desenvolvimento

Clone o repositório e instale em modo editável (`-e`), que vincula diretamente
ao código-fonte sem necessidade de reinstalar a cada alteração:

```bash
git clone https://github.com/seu_usuario/hsp_2026
cd hsp_2026
pip install -e .
```

---

## Instalação em ambiente virtual (recomendado)

O uso de ambientes virtuais isola as dependências do projeto e evita conflitos
com outros pacotes do sistema.

=== "venv (padrão Python)"

    ```bash
    python -m venv .venv
    source .venv/bin/activate        # Linux / macOS
    .venv\Scripts\activate           # Windows

    pip install hsp_2026
    ```

=== "conda"

    ```bash
    conda create -n hsp_2026 python=3.10
    conda activate hsp_2026

    pip install hsp_2026
    ```

---

## Verificação da instalação

Após instalar, execute o seguinte no terminal Python para confirmar que todos
os módulos estão acessíveis:

```python
import hsp_2026.file as file
import hsp_2026.prep as prep
import hsp_2026.emsc as emsc
import hsp_2026.sh   as sh
import hsp_2026.km   as km
import hsp_2026.qt   as qt
import hsp_2026.rg   as rg

print("HSP instalado com sucesso!")
```

---

## Dependências para a documentação

Se desejar compilar e servir esta documentação localmente:

```bash
pip install mkdocs-material "mkdocstrings[python]"
mkdocs serve
```

A documentação ficará disponível em `http://127.0.0.1:8000`.
