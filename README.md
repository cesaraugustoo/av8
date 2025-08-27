# **Implementação do Método de Hückel**

## **Resumo**

Este repositório contém um pipeline em Python que aplica o **Método de Hückel** para resolver problemas de química quântica, permitindo:

* Cálculo dos **níveis de energia** dos orbitais moleculares (OMs).
* Determinação dos **coeficientes dos OMs**.
* Visualização de propriedades eletrônicas fundamentais, como:

  * **HOMO** (*Highest Occupied Molecular Orbital* – Orbital Molecular Ocupado de Mais Alta Energia).
  * **LUMO** (*Lowest Unoccupied Molecular Orbital* – Orbital Molecular Desocupado de Mais Baixa Energia).

O código constrói a **matriz de Hückel**, resolve o **problema de autovalores** para obter energias e orbitais e gera **diagramas e relatórios detalhados**.

---

## **Problemas Resolvidos**

O pipeline implementado no script `moleculas.py` aborda diretamente as tarefas propostas:

### **Tarefa 1: Análise do Oxepino (C₆H₆O)**

O script realiza uma análise completa de Hückel para o Oxepino. Ele calcula os níveis de energia dos elétrons π, identifica os orbitais **HOMO** e **LUMO**, e gera visualizações dos orbitais de fronteira. O relatório final fornece os dados necessários para discutir a natureza eletrônica da molécula, incluindo seu gap de energia de **2.077 eV**.

### **Tarefa 2: Estudo Comparativo do Trifenileno e Aza-derivados**

Para a segunda tarefa, o pipeline analisa sistematicamente o Trifenileno e seus análogos substituídos por nitrogênio. O método de Hückel avançado ajusta os parâmetros de Coulomb ($\alpha$) e ressonância ($\beta$) para os heteroátomos, permitindo uma comparação direta e significativa. Os relatórios gerados permitem analisar como o número de átomos de nitrogênio afeta a distribuição eletrônica, as ordens de ligação e, crucialmente, os orbitais de fronteira.

---

## **Pré-requisitos**

Para executar este projeto, é necessário ter **Python 3.x** instalado, juntamente com as seguintes bibliotecas:

* **Python 3.10+**
* **NumPy**: cálculos numéricos e manipulação de matrizes.
* **Matplotlib**: geração de gráficos e diagramas de energia.
* **RDKit**: manipulação de estruturas moleculares e identificação de sistemas π.

---

## **Instalação das Dependências com Conda**

A instalação se dá via **conda**:

### **1. Criar e ativar um novo ambiente Conda**

```bash
conda create -c conda-forge -n huckel_env python=3.10 rdkit numpy matplotlib
```

### **2. Ativar o ambiente**

```bash
conda activate huckel_env
```

---

## **Como Executar**

O projeto foi desenvolvido como um **pipeline autocontido** no script `huckel_pipeline.py`.
O script principal para resolver as tarefas propostas é o `moleculas.py`. Para executá-lo, utilize o seguinte comando no seu terminal (com o ambiente Conda ativado):

```bash
python moleculas.py
```

Este comando irá executar a análise de Hückel para as quatro moléculas definidas no script e gerar todos os arquivos de saída no diretório `output/`.

---

## **Arquivos Gerados**

Após a execução, será criado um diretório chamado **`demo_output/`**, contendo:

* **Relatório HTML (.html)**: relatório completo com diagramas, tabelas de energia, ordens de ligação e mapas de orbitais de fronteira.
* **Imagens (.png)**:
  * Diagrama de níveis de energia dos orbitais moleculares.
  * Visualização dos coeficientes dos orbitais HOMO e LUMO sobre a estrutura molecular.
  * Mapa de ordens de ligação π.
  * Estrutura molecular com a numeração dos átomos do sistema π.
* **Dados JSON (.json)**: arquivo com resultados numéricos (energias, coeficientes, etc.) para processamento.

---

## **Estrutura do Projeto**

O projeto está concentrado em um único arquivo principal para simplicidade:

* **`huckel_pipeline.py`**:

  * **Classes principais**:

    * `MoleculeConfig`, `Results`: manipulação de dados.
    * `HuckelCalculator`: lógica do método (construção e diagonalização da matriz).
    * `Visualizer`: geração de gráficos e imagens.
    * `ReportGenerator`: criação de relatórios em HTML e JSON.
    * `QuantumChemistryPipeline`: orquestra o fluxo completo.

---

## **Funcionalidades Principais**

  * **Construção da Matriz de Hückel**: gerada a partir da topologia molecular identificada via RDKit, com parâmetros ajustáveis para heteroátomos.
  * **Cálculo de Autovalores e Autovetores**: diagonalização da matriz para obter energias e coeficientes dos OMs.
  * **Análise dos Orbitais de Fronteira**: identificação de **HOMO**, **LUMO** e cálculo do **gap de energia**.
  * **Propriedades Eletrônicas**: cálculo de ordens de ligação π e densidades de carga π.
  * **Visualizações e Relatórios**: gráficos em PNG e relatórios completos em HTML para cada molécula analisada.

---

## **Diagrama de Fluxo**

```mermaid
flowchart TD

%% ============================
%% SEÇÃO 1 - NÍVEL ALTO
%% ============================
subgraph HighLevel["Fluxo Geral do Pipeline Huckel"]
    A1["1. Entrada da Molécula (SMILES/arquivo)"]
    A2["2. Geração da Estrutura Molecular (RDKit)"]
    A3["3. Construção da Matriz de Adjacência"]
    A4["4. Montagem da Matriz Huckel"]
    A5["5. Cálculo de Autovalores e Autovetores"]
    A6["6. Identificação de HOMO e LUMO"]
    A7["7. Visualizações e Resultados"]
end

A1 --> A2 --> A3 --> A4 --> A5 --> A6 --> A7

%% ============================
%% SEÇÃO 2 - DETALHES DO CÓDIGO
%% ============================
subgraph Details["Detalhes Internos do Código"]
    B1["Importação das Bibliotecas: numpy, rdkit, matplotlib"]
    B2["Leitura da molécula via SMILES (Chem.MolFromSmiles)"]
    B3["Conversão para grafo -> Matriz de adjacência"]
    B4["Definição de parâmetros α (alpha) e β (beta) do Huckel"]
    B5["Construção da Matriz Hamiltoniana H = αI + βA"]
    B6["Diagonalização: numpy.linalg.eigh"]
    B7["Ordenação dos níveis de energia"]
    B8["Determinação de HOMO/LUMO pela contagem de elétrons"]
    B9["Plot de:
        - Estrutura molecular
        - Diagrama de níveis de energia
        - Orbitais moleculares"]
    B10["Geração dos arquivos de saída (.png, .csv)"]
end

%% RELAÇÕES ENTRE ALTO E DETALHE
A1 -->|Implementado por| B2
A2 -->|Visualização com RDKit| B9
A3 -->|Implementado por| B3
A4 -->|Baseado em| B4 & B5
A5 -->|Implementado por| B6 & B7
A6 -->|Baseado em| B8
A7 -->|Gerado por| B9 & B10
```

---

## **Etapas de Visualização**

* **Mapa de Ordem de Ligação π** (`plot_bond_order_overlay`): sobrepõe as ordens de ligação π na estrutura 2D.
* **Diagrama de Níveis de Energia** (`plot_orbital_energy_diagram`): exibe níveis, HOMO, LUMO e gap.
* **Mapas de Coeficientes HOMO/LUMO** (`plot_huckel_mo_on_molecule`): desenha orbitais sobre a molécula.

---

## **Exemplo de Uso**

```python
from huckel_pipeline import QuantumChemistryPipeline, MoleculeConfig

# 1. Inicializar o pipeline
pipeline = QuantumChemistryPipeline(output_dir="meus_calculos")

# 2. Definir a molécula (SMILES do benzeno)
molecula_benzeno = MoleculeConfig(smiles="C1=CC=CC=C1", name="benzeno")

# 3. Executar o cálculo
resultados = pipeline.run_single_calculation(molecula_benzeno)

# 4. Exibir resultados principais
if resultados.success:
    print(f"Cálculo para {resultados.molecule_name} concluído com sucesso!")
    print(f"Energia do HOMO: {resultados.homo_energy:.3f} eV")
    print(f"Energia do LUMO: {resultados.lumo_energy:.3f} eV")
    print(f"Gap HOMO-LUMO: {resultados.homo_lumo_gap:.3f} eV")
```

## **Exemplo de Uso (Adicionando uma Nova Molécula)**

Para analisar uma nova molécula, basta adicioná-la à lista `molecules_to_run` no arquivo `moleculas.py`:

```python
from huckel_pipeline import QuantumChemistryPipeline, MoleculeConfig

# 1. Inicializar o pipeline
pipeline = QuantumChemistryPipeline(output_dir="output")

# 2. Definir a lista de moléculas
molecules_to_run = [
    # ... moléculas existentes
    
    # Adicione uma nova molécula aqui
    MoleculeConfig(
        smiles="C1=CC2=C(C=C1)C=C3C=CC=C3C2", 
        name="Anthracene"
    )
]

# 3. Executar o cálculo em lote
resultados = pipeline.run_batch_calculations(molecules_to_run)
```