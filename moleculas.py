from huckel_pipeline import QuantumChemistryPipeline, MoleculeConfig

# --- Inicializa o Pipeline ---
pipeline = QuantumChemistryPipeline("output")

# --- Define as Moléculas para Análise ---
molecules_to_run = [
    MoleculeConfig(
        smiles="C1=CC=COC=C1", 
        name="Oxepin"
    ),
    MoleculeConfig(
        smiles="C1=CC=C2C(=C1)C3=CC=CC=C3C4=CC=CC=C24", 
        name="Triphenylene"
    ),
    MoleculeConfig(
        smiles="C1=CC=C2C(=C1)C3=CC=CC=C3C4=NN=NC=C24", 
        name="Triazatriphenylene"
    ),
    MoleculeConfig(
        smiles="C1=CN=C2C(=N1)C3=CN=CN=C3C4=NC=CN=C24", 
        name="Hexaazatriphenylene"
    )
]

# --- Executa os Cálculos ---
print("="*50)
print("🔬 Iniciando Análise pelo Método de Hückel para Múltiplas Moléculas")
print("="*50)

all_results = pipeline.run_batch_calculations(molecules_to_run)

# --- Resumo Final ---
print("\n" + "="*50)
print("Processamento concluído.")
print(f"Verifique o diretório '{pipeline.output_dir}' para todos os relatórios e imagens gerados.")

# Exibe uma tabela de resumo no console para referência rápida.
print("\nResumo dos Cálculos:")
print("-" * 80)
print(f"{'Molécula':<25} {'Elétrons π':<12} {'HOMO (eV)':<12} {'LUMO (eV)':<12} {'Gap (eV)':<10} {'Status'}")
print("-" * 80)

for result in all_results:
    homo = f"{result.homo_energy:.3f}" if result.homo_energy is not None else "N/D"
    lumo = f"{result.lumo_energy:.3f}" if result.lumo_energy is not None else "N/D"
    gap = f"{result.homo_lumo_gap:.3f}" if result.homo_lumo_gap is not None else "N/D"
    status = "✓" if result.success else "✗"
    
    print(f"{result.molecule_name:<25} {str(result.num_pi_electrons):<12} {homo:<12} {lumo:<12} {gap:<10} {status}")