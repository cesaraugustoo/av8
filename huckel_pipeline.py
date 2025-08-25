import os
import sys
import json
import logging
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union, Any
from dataclasses import dataclass, asdict
import traceback

# Bibliotecas científicas principais
import numpy as np
import pandas as pd
from scipy import linalg
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import matplotlib.patches as patches
import seaborn as sns
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

# Manipulação molecular
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
    from rdkit.Chem import Draw
    from rdkit.Chem.Draw import rdMolDraw2D
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    warnings.warn("RDKit não está disponível. Algumas funcionalidades serão limitadas.")

# Configura o logging
logging.basicConfig(level=logging.INFO, 
                   format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# ============================================================================
# Classes de Configuração e Dados
# ============================================================================

@dataclass
class MoleculeConfig:
    """Configuração para o sistema molecular."""
    smiles: str = ""
    name: str = "desconhecido"
    charge: int = 0

@dataclass
class HuckelConfig:
    """Configuração para o método de Hückel baseada em parâmetros de referência."""
    alpha_0: float = 0.0   # Valor base de alfa para Carbono (eV)
    beta_0: float = -2.5   # Valor base de beta para ligação aromática C-C (eV)

@dataclass
class Results:
    """Contêiner para os resultados do cálculo."""
    molecule_name: str
    input_smiles: str
    energy: Optional[float] = None
    homo_energy: Optional[float] = None
    lumo_energy: Optional[float] = None
    homo_lumo_gap: Optional[float] = None
    molecular_orbitals: Optional[np.ndarray] = None
    orbital_energies: Optional[np.ndarray] = None
    calculation_method: str = "desconhecido"
    success: bool = False
    error_message: str = ""

    # Resultados específicos de Hückel
    connectivity: Optional[Dict[int, List[int]]] = None
    atom_populations: Optional[np.ndarray] = None
    pi_charges: Optional[np.ndarray] = None
    bond_orders: Optional[Dict[Tuple[int, int], float]] = None
    pi_atoms: Optional[List[int]] = None
    num_pi_electrons: Optional[int] = None
    num_occupied_orbitals: Optional[int] = None
    homo_index: Optional[int] = None
    lumo_index: Optional[int] = None
    somo_index: Optional[int] = None
    somo_energy: Optional[float] = None
    huckel_matrix: Optional[np.ndarray] = None

    # Resultados de planaridade
    planarity_status: Optional[str] = None
    planarity_rms: Optional[float] = None
    planarity_max_dev: Optional[float] = None
    planarity_method: Optional[str] = None

# ============================================================================
# Módulo de Manipulação de Entrada
# ============================================================================

class InputHandler:
    """Manipula a entrada molecular e a geração de estruturas."""
    
    def __init__(self):
        self.logger = logging.getLogger(f"{__name__}.InputHandler")
        if not RDKIT_AVAILABLE:
            raise ImportError("RDKit é necessário para o InputHandler")
    
    def smiles_to_mol(self, smiles: str, name: str = "molecule") -> Chem.Mol:
        """Converte uma string SMILES em um objeto Mol do RDKit."""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError(f"SMILES inválido: {smiles}")
            
            mol.SetProp("_Name", name)
            self.logger.info(f"SMILES processado com sucesso: {smiles}")
            return mol
            
        except Exception as e:
            self.logger.error(f"Erro ao processar o SMILES {smiles}: {e}")
            raise
    
    def get_molecular_properties(self, mol: Chem.Mol) -> Dict[str, Any]:
        """Calcula propriedades moleculares básicas."""
        try:
            props = {
                "massa_molecular": Descriptors.MolWt(mol),
                "num_atomos": mol.GetNumAtoms(),
                "num_atomos_pesados": mol.GetNumHeavyAtoms(),
                "num_ligacoes_rotaveis": Descriptors.NumRotatableBonds(mol),
                "tpsa": Descriptors.TPSA(mol),
                "logp": Descriptors.MolLogP(mol),
                "num_aneis_aromaticos": Descriptors.NumAromaticRings(mol),
                "formula": rdMolDescriptors.CalcMolFormula(mol)
            }
            
            self.logger.info("Propriedades moleculares calculadas")
            return props
            
        except Exception as e:
            self.logger.error(f"Erro ao calcular propriedades: {e}")
            return {}

# ============================================================================
# Implementação do Método de Hückel
# ============================================================================

class HuckelCalculator:
    """
    Implementação do método de Hückel usando parâmetros específicos do ambiente químico.
    """
    
    def __init__(self, config: HuckelConfig = None):
        self.config = config or HuckelConfig()
        self.logger = logging.getLogger(f"{__name__}.HuckelCalculator")

        if self.config.beta_0 > 0:
            self.logger.warning(f"beta_0 é positivo ({self.config.beta_0}). "
                                f"A convenção padrão é beta_0 < 0. Os níveis de energia podem estar invertidos.")

        self.h_params = {
            'B': -1.0, 'C': 0.0, 'N_pyridine': 0.5, 'N_pyrrole': 1.5, 'N_positive': 2.0,
            'O_carbonyl': 1.0, 'O_furan': 2.0, 'O_positive': 2.5, 'F': 3.0, 'Cl': 2.0, 'Br': 1.5,
        }
        self.k_params = {
            'C-C': 0.9, 'C=C': 1.1, 'C_aromatic': 1.0, 'C-B': 0.7, 'C-N': 0.8,
            'C_N_aromatic': 1.0, 'C-O': 0.8, 'C=O': 1.0, 'N-O': 0.7, 'C-F': 0.7,
            'C-Cl': 0.4, 'C-Br': 0.3,
        }

    def _get_atom_type(self, atom: Chem.Atom) -> str:
        symbol = atom.GetSymbol()
        if symbol in ['C', 'B', 'F', 'Cl', 'Br']: return symbol
        if symbol == 'N':
            if atom.GetFormalCharge() == 1: return 'N_positive'
            if atom.GetIsAromatic(): return 'N_pyrrole' if atom.GetTotalNumHs() > 0 else 'N_pyridine'
            if any(b.GetBondType() == Chem.BondType.DOUBLE for b in atom.GetBonds()): return 'N_pyridine'
            return 'N_pyrrole'
        if symbol == 'O':
            if atom.GetFormalCharge() == 1: return 'O_positive'
            return 'O_carbonyl' if len(atom.GetBonds()) == 1 else 'O_furan'
        return symbol

    def _get_bond_type(self, bond: Chem.Bond) -> str:
        atom1, atom2 = bond.GetBeginAtom(), bond.GetEndAtom()
        symbols = sorted([atom1.GetSymbol(), atom2.GetSymbol()])
        bond_key = f"{symbols[0]}-{symbols[1]}"
        if bond.GetIsAromatic():
            if bond_key == "C-C": return "C_aromatic"
            if bond_key == "C-N": return "C_N_aromatic"
            return bond_key
        bond_type = bond.GetBondType()
        if bond_type == Chem.BondType.SINGLE: return bond_key
        if bond_type == Chem.BondType.DOUBLE: return f"{symbols[0]}={symbols[1]}"
        return "unknown"

    def identify_pi_system(self, mol: Chem.Mol) -> List[int]:
        pi_atoms = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic() or atom.GetHybridization() == Chem.HybridizationType.SP2}
        self.logger.info(f"Identificados {len(pi_atoms)} átomos no sistema π")
        return sorted(list(pi_atoms))

    def build_huckel_matrix(self, mol: Chem.Mol, pi_atoms: List[int]) -> np.ndarray:
        n = len(pi_atoms)
        if n == 0: return np.array([[]])
        H = np.zeros((n, n))
        atom_map = {atom_idx: i for i, atom_idx in enumerate(pi_atoms)}
        for i, atom_idx in enumerate(pi_atoms):
            atom = mol.GetAtomWithIdx(atom_idx)
            atom_type = self._get_atom_type(atom)
            h_val = self.h_params.get(atom_type, 0.0)
            H[i, i] = self.config.alpha_0 + h_val * self.config.beta_0
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() in atom_map:
                    j = atom_map[neighbor.GetIdx()]
                    if i < j:
                        bond = mol.GetBondBetweenAtoms(atom_idx, neighbor.GetIdx())
                        bond_type = self._get_bond_type(bond)
                        k_val = self.k_params.get(bond_type, self.k_params.get(bond_type.replace('=', '-'), 0.0))
                        H[i, j] = H[j, i] = k_val * self.config.beta_0
        self.logger.info(f"Matriz de Hückel {n}x{n} construída.")
        return H

    def solve_huckel(self, mol: Chem.Mol, mol_config: MoleculeConfig) -> Results:
        results = Results(molecule_name=mol_config.name, input_smiles=mol_config.smiles, calculation_method="Hückel (Avançado)")
        try:
            pi_atoms = self.identify_pi_system(mol)
            if not pi_atoms:
                results.error_message = "Nenhum sistema π encontrado na molécula"
                return results
            H = self.build_huckel_matrix(mol, pi_atoms)
            try:
                results.huckel_matrix = H.copy() if isinstance(H, np.ndarray) else np.asarray(H)
            except Exception:
                results.huckel_matrix = None
            if H.size == 0:
                results.error_message = "Falha ao construir a matriz de Hückel"
                return results
            
            eigenvalues, eigenvectors = linalg.eigh(H)
            idx = np.argsort(eigenvalues)
            results.orbital_energies, results.molecular_orbitals = eigenvalues[idx], eigenvectors[:, idx]
            results.pi_atoms, results.success = pi_atoms, True
            
            pi_electron_supply = np.zeros(len(pi_atoms))
            for i, atom_idx in enumerate(pi_atoms):
                atom = mol.GetAtomWithIdx(atom_idx)
                symbol, charge = atom.GetSymbol(), atom.GetFormalCharge()
                electrons = 1 - charge if symbol == 'C' else (2 if self._get_atom_type(atom) in ['N_pyrrole', 'O_furan'] else (0 if self._get_atom_type(atom) == 'B' else 1))
                pi_electron_supply[i] = max(0, min(2, electrons))
            
            num_pi_electrons = int(np.sum(pi_electron_supply))
            results.num_pi_electrons = num_pi_electrons
            
            conn, pops, charges, orders = self._calculate_huckel_properties(mol, pi_atoms, eigenvectors[:, idx], num_pi_electrons, pi_electron_supply)
            results.connectivity, results.atom_populations, results.pi_charges, results.bond_orders = conn, pops, charges, orders
            
            is_odd_electron = num_pi_electrons % 2 != 0
            occupied_orbitals = num_pi_electrons // 2
            results.num_occupied_orbitals = occupied_orbitals

            if is_odd_electron:
                somo_idx = occupied_orbitals
                results.somo_index, results.somo_energy = somo_idx, eigenvalues[idx[somo_idx]]
                results.energy = 2 * np.sum(eigenvalues[idx[:occupied_orbitals]]) + eigenvalues[idx[somo_idx]]
                if somo_idx > 0: results.homo_index, results.homo_energy = somo_idx - 1, eigenvalues[idx[somo_idx - 1]]
                if somo_idx + 1 < len(eigenvalues):
                    results.lumo_index, results.lumo_energy = somo_idx + 1, eigenvalues[idx[somo_idx + 1]]
                    results.homo_lumo_gap = results.lumo_energy - results.somo_energy
            else:
                if occupied_orbitals > 0:
                    results.homo_index, results.homo_energy = occupied_orbitals - 1, eigenvalues[idx[occupied_orbitals - 1]]
                    results.energy = 2 * np.sum(eigenvalues[idx[:occupied_orbitals]])
                    if occupied_orbitals < len(eigenvalues):
                        results.lumo_index, results.lumo_energy = occupied_orbitals, eigenvalues[idx[occupied_orbitals]]
                        results.homo_lumo_gap = results.lumo_energy - results.homo_energy
            
            self.logger.info(f"Cálculo de Hückel concluído para {num_pi_electrons} elétrons π.")
        except Exception as e:
            error_msg = f"Cálculo de Hückel falhou: {str(e)}"
            self.logger.error(f"{error_msg}\n{traceback.format_exc()}")
            results.error_message, results.success = error_msg, False
        return results
    
    def _calculate_huckel_properties(self, mol, pi_atoms, eigenvectors, num_pi_electrons, pi_electron_supply):
        atom_map = {atom_idx: i for i, atom_idx in enumerate(pi_atoms)}
        connectivity = {i: [atom_map[n.GetIdx()] for n in mol.GetAtomWithIdx(pi_atoms[i]).GetNeighbors() if n.GetIdx() in atom_map] for i in range(len(pi_atoms))}
        num_occupied = num_pi_electrons // 2
        is_odd = num_pi_electrons % 2 != 0
        populations = 2 * np.sum(eigenvectors[:, :num_occupied] ** 2, axis=1) if num_occupied > 0 else np.zeros(len(pi_atoms))
        if is_odd: populations += eigenvectors[:, num_occupied] ** 2
        bond_orders = {}
        for i in range(len(pi_atoms)):
            for j in range(i + 1, len(pi_atoms)):
                if j in connectivity[i]:
                    order = 2 * np.sum(eigenvectors[i, :num_occupied] * eigenvectors[j, :num_occupied]) if num_occupied > 0 else 0.0
                    if is_odd: order += eigenvectors[i, num_occupied] * eigenvectors[j, num_occupied]
                    bond_orders[(i, j)] = order
        pi_charges = pi_electron_supply - populations
        return connectivity, populations, pi_charges, bond_orders

# ============================================================================
# Módulo de Visualização
# ============================================================================

class Visualizer:
    """Cria visualizações científicas dos resultados."""
    
    def __init__(self, output_dir: str = "output"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.logger = logging.getLogger(f"{__name__}.Visualizer")
        plt.style.use('seaborn-v0_8-whitegrid')
        sns.set_palette("husl")
    
    def plot_molecular_structure(self, mol: Chem.Mol, title: str = "Molecule") -> str:
        if not RDKIT_AVAILABLE: return ""
        try:
            AllChem.Compute2DCoords(mol)
            drawer = rdMolDraw2D.MolDraw2DCairo(800, 600)
            drawer.DrawMolecule(mol)
            drawer.FinishDrawing()
            filename = self.output_dir / f"{title.replace(' ', '_')}_estrutura.png"
            with open(filename, 'wb') as f: f.write(drawer.GetDrawingText())
            self.logger.info(f"Estrutura molecular salva: {filename}")
            return str(filename)
        except Exception as e:
            self.logger.error(f"Erro ao plotar a estrutura molecular: {e}"); return ""
    
    def plot_orbital_energy_diagram(self, results: Results, save_name: str) -> str:
        if results.orbital_energies is None: return ""
        try:
            fig, ax = plt.subplots(figsize=(6, 8))
            energies = results.orbital_energies
            num_occupied = results.num_occupied_orbitals
            is_odd_electron = results.somo_index is not None
            for i, energy in enumerate(energies):
                color = 'blue' if i < num_occupied else ('green' if is_odd_electron and i == results.somo_index else 'red')
                linewidth = 3 if is_odd_electron and i == results.somo_index else 2
                ax.hlines(energy, -0.4, 0.4, colors=color, linewidth=linewidth)
                ax.text(0.5, energy, f'{energy:.3f}', ha='left', va='center', fontsize=9)
            if results.homo_energy is not None and results.lumo_energy is not None:
                start_e = results.somo_energy if is_odd_electron else results.homo_energy
                mid_point = (start_e + results.lumo_energy) / 2
                gap_label = f'Gap S-L:\n{results.homo_lumo_gap:.3f} eV' if is_odd_electron else f'Gap:\n{results.homo_lumo_gap:.3f} eV'
                ax.annotate("", xy=(0, start_e), xytext=(0, results.lumo_energy), arrowprops=dict(arrowstyle='<->', color='purple', lw=2))
                ax.text(-0.1, mid_point, gap_label, ha='right', va='center', bbox=dict(boxstyle='round,pad=0.5', facecolor='lavender', alpha=0.7))
            ax.set_ylabel('Energia (eV)'); ax.set_title(f'Diagrama de OM - {results.molecule_name}'); ax.set_xticks([])
            legend_handles = [patches.Patch(color='blue', label='Ocupado')]
            if is_odd_electron: legend_handles.append(patches.Patch(color='green', label='SOMO'))
            legend_handles.append(patches.Patch(color='red', label='Virtual'))
            ncols = len(legend_handles)
            ax.legend(
               handles=legend_handles,
               loc='upper center',
               bbox_to_anchor=(0.5, -0.08),
               ncol=ncols,
               frameon=True,
               fontsize=10,
               bbox_transform=fig.transFigure
            )
            plt.tight_layout()
            plt.subplots_adjust(bottom=0.01)
            filename = self.output_dir / f"{save_name.replace(' ', '_')}.png"
            plt.savefig(filename, dpi=300, bbox_inches='tight'); plt.close()
            self.logger.info(f"Diagrama de energia dos orbitais salvo: {filename}")
            return str(filename)
        except Exception as e:
            self.logger.error(f"Erro ao criar o diagrama de energia dos orbitais: {e}"); return ""

    def plot_orbital_coefficients(self, results: Results, orbital_idx: int = None, save_name: str = None) -> str:
        if results.molecular_orbitals is None: return ""
        try:
            orbitals = results.molecular_orbitals
            energies = results.orbital_energies
            if orbital_idx is None:
                orbital_idx = results.somo_index if results.somo_index is not None else (results.homo_index if results.homo_index is not None else orbitals.shape[1] // 2 - 1)
            coefficients = orbitals[:, orbital_idx]
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
            atom_indices = range(len(coefficients))
            colors = ['red' if c < 0 else 'blue' for c in coefficients]
            ax1.bar(atom_indices, coefficients, color=colors, alpha=0.7)
            ax1.set_title(f'Coeficientes do OM {orbital_idx + 1}\nEnergia: {energies[orbital_idx]:.3f} eV')
            ax2.plot(atom_indices, np.abs(coefficients), 'o-')
            ax2.set_title('Magnitude dos Coeficientes')
            if save_name is None: save_name = f"{results.molecule_name}_om_{orbital_idx + 1}_coeficientes"
            filename = self.output_dir / f"{save_name.replace(' ', '_')}.png"
            plt.savefig(filename, dpi=300); plt.close()
            return str(filename)
        except Exception as e:
            self.logger.error(f"Erro ao plotar os coeficientes dos orbitais: {e}"); return ""
    
    def create_interactive_plot(self, results_list: List[Results], property_name: str = "homo_lumo_gap") -> str:
        try:
            data = [{'molecula': r.molecule_name, 'valor': getattr(r, property_name), 'metodo': r.calculation_method} for r in results_list if r.success and hasattr(r, property_name) and getattr(r, property_name) is not None]
            if not data: return ""
            df = pd.DataFrame(data)
            fig = px.bar(df, x='molecula', y='valor', color='valor', text='valor', title=f'Comparação de {property_name.replace("_", " ").title()}', labels={'molecula': 'Molécula', 'valor': property_name.replace('_', ' ').title()})
            filename = self.output_dir / f"comparacao_interativa_{property_name}.html"
            fig.write_html(filename)
            return str(filename)
        except Exception as e:
            self.logger.error(f"Erro ao criar o gráfico interativo: {e}"); return ""
        
    def plot_huckel_mo_on_molecule(self, results: Results, mol: Chem.Mol, orbital_idx: int, pi_atoms: List[int], save_name: str = None) -> str:
        if not RDKIT_AVAILABLE or results.molecular_orbitals is None: return ""
        try:
            mol_to_draw = Chem.Mol(mol); AllChem.Compute2DCoords(mol_to_draw)
            coeffs = results.molecular_orbitals[:, orbital_idx]
            drawer = rdMolDraw2D.MolDraw2DCairo(500, 500)
            drawer.drawOptions().addAtomIndices = True
            atom_colors = {atom_idx: ((0, 0, 1) if coeffs[i] > 0.05 else (1, 0, 0)) for i, atom_idx in enumerate(pi_atoms) if abs(coeffs[i]) > 0.05}
            max_abs_coeff = np.max(np.abs(coeffs))
            radii = {atom_idx: 0.1 + 0.4 * (abs(coeffs[i]) / max_abs_coeff) if max_abs_coeff > 1e-6 else 0.1 for i, atom_idx in enumerate(pi_atoms)}
            drawer.DrawMolecule(mol_to_draw, highlightAtoms=pi_atoms, highlightAtomColors=atom_colors, highlightAtomRadii=radii)
            drawer.FinishDrawing()
            if save_name is None:
                orbital_name = "HOMO" if orbital_idx == results.homo_index else "LUMO" if orbital_idx == results.lumo_index else f"OM_{orbital_idx+1}"
                save_name = f"{results.molecule_name}_{orbital_name}_coeffs"
            filename = self.output_dir / f"{save_name.replace(' ', '_')}.png"
            with open(filename, 'wb') as f: f.write(drawer.GetDrawingText())
            self.logger.info(f"Gráfico do OM de Hückel salvo: {filename}")
            return str(filename)
        except Exception as e:
            self.logger.error(f"Erro ao plotar o OM de Hückel na molécula: {e}"); return ""
        
    def plot_numbered_molecule(self, mol: Chem.Mol, pi_atoms: List[int], save_name: str) -> str:
        if not RDKIT_AVAILABLE: return ""
        try:
            mol_to_draw = Chem.Mol(mol); AllChem.Compute2DCoords(mol_to_draw)
            drawer = rdMolDraw2D.MolDraw2DCairo(500, 500)
            drawer.drawOptions().addAtomIndices = True
            drawer.drawOptions().highlightColour = (0.85, 0.85, 0.85) 
            drawer.DrawMolecule(mol_to_draw, highlightAtoms=pi_atoms)
            drawer.FinishDrawing()
            filename = self.output_dir / f"{save_name.replace(' ', '_')}_numerada.png"
            with open(filename, 'wb') as f: f.write(drawer.GetDrawingText())
            self.logger.info(f"Esquema da molécula numerada salvo: {filename}")
            return str(filename)
        except Exception as e:
            self.logger.error(f"Erro ao plotar a molécula numerada: {e}"); return ""

    def compute_pi_planarity(self, mol: Chem.Mol, pi_atoms: List[int]) -> Dict[str, Any]:
        """Calcula uma métrica de planaridade para o sistema π."""
        coords, method = None, "desconhecido"
        mol_copy = Chem.Mol(mol)
        if mol_copy.GetNumConformers() > 0:
            conf = mol_copy.GetConformer()
            method = "3D-existente"
        else:
            try:
                AllChem.EmbedMolecule(mol_copy, randomSeed=42)
                conf = mol_copy.GetConformer()
                method = "3D-embutido"
            except Exception:
                AllChem.Compute2DCoords(mol_copy)
                method = "2D-assumido-planar"
        
        if method == "2D-assumido-planar":
            return dict(rms_dev=0.0, max_dev=0.0, is_planar="planar", method=method)

        conf = mol_copy.GetConformer()
        coords = np.array([conf.GetAtomPosition(i) for i in pi_atoms])
        centroid = coords.mean(axis=0)
        X = coords - centroid
        _, _, vh = np.linalg.svd(X, full_matrices=False)
        normal = vh[-1, :]
        dists = np.abs(X @ normal)
        rms_dev, max_dev = np.sqrt((dists**2).mean()), dists.max()
        
        status = "nao-planar"
        if rms_dev <= 0.10 and max_dev <= 0.20: status = "planar"
        elif rms_dev <= 0.20 and max_dev <= 0.35: status = "quase-planar"
        
        self.logger.info(f"Verificação de planaridade ({method}): status={status}, Desvio RMS={rms_dev:.4f} Å, Desvio Máx={max_dev:.4f} Å")
        return dict(rms_dev=rms_dev, max_dev=max_dev, is_planar=status, method=method)

    def plot_bond_order_overlay(self, results: Results, mol: Chem.Mol, save_name: str, min_label: float = 0.05) -> str:
        """Renderiza uma sobreposição 2D das ordens de ligação π na estrutura molecular."""
        if not results.bond_orders or not results.pi_atoms:
            self.logger.warning("Nenhum dado de ordem de ligação disponível para gerar a sobreposição.")
            return ""
        try:
            mol_copy = Chem.Mol(mol)
            AllChem.Compute2DCoords(mol_copy)
            conf = mol_copy.GetConformer()
            
            fig, ax = plt.subplots(figsize=(8, 8))
            
            segments, values = [], []
            for (i, j), v in results.bond_orders.items():
                a, b = results.pi_atoms[i], results.pi_atoms[j]
                if mol_copy.GetBondBetweenAtoms(a, b) is None: continue
                p1, p2 = conf.GetAtomPosition(a), conf.GetAtomPosition(b)
                segments.append(((p1.x, p1.y), (p2.x, p2.y)))
                values.append(v)

            if not values: return ""
            
            vmin, vmax = min(values), max(values)
            norm = Normalize(vmin=vmin, vmax=vmax)
            cmap = plt.cm.viridis

            for ((x1, y1), (x2, y2)), v in zip(segments, values):
                t = 0.5 if vmax == vmin else norm(v)
                width = 2 + 6 * t
                color = cmap(t)
                ax.plot([x1, x2], [y1, y2], color=color, linewidth=width, solid_capstyle='round', zorder=1)
                if v >= min_label:
                    xm, ym = (x1 + x2) / 2, (y1 + y2) / 2
                    ax.text(xm, ym, f'{v:.2f}', fontsize=8, ha='center', va='center', color='black',
                            bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7), zorder=3)

            for atom in mol_copy.GetAtoms():
                pos = conf.GetAtomPosition(atom.GetIdx())
                ax.text(pos.x, pos.y, atom.GetSymbol(), fontsize=10, ha='center', va='center',
                        bbox=dict(boxstyle='circle,pad=0.3', facecolor='white', edgecolor='none'), zorder=2)

            ax.set_aspect('equal'); ax.axis('off')
            fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax, label='Ordem de Ligação π', shrink=0.8)
            plt.title(f"Mapa de Ordem de Ligação π - {results.molecule_name}")
            plt.tight_layout()
            
            filename = self.output_dir / f"{save_name.replace(' ', '_')}_ordens_ligacao.png"
            plt.savefig(filename, dpi=300)
            plt.close()
            
            self.logger.info(f"Sobreposição de ordem de ligação salva: {filename}")
            return str(filename)
        except Exception as e:
            self.logger.error(f"Erro ao plotar a sobreposição de ordem de ligação: {e}"); return ""
        
    # def plot_huckel_matrix_heatmap(self, H_matrix, pi_atoms, save_name: str) -> str:
    #     """
    #     Visualiza a matriz de Hückel como um mapa de calor.
    #     """
    #     try:
    #         if H_matrix is None:
    #             self.logger.warning("plot_huckel_matrix_heatmap: received None H_matrix, skipping.")
    #             return ""

    #         arr = np.asarray(H_matrix)

    #         if arr.ndim == 0:
    #             self.logger.warning(f"plot_huckel_matrix_heatmap: H_matrix is 0-d (scalar). shape={arr.shape}; skipping.")
    #             return ""

    #         if arr.size == 0 or arr.shape[0] == 0 or arr.shape[1] == 0:
    #             self.logger.warning(f"plot_huckel_matrix_heatmap: H_matrix is empty or has a zero dimension. shape={arr.shape}; skipping.")
    #             return ""

    #         if arr.ndim != 2:
    #             self.logger.warning(f"plot_huckel_matrix_heatmap: H_matrix is not 2-D (ndim={arr.ndim}); skipping.")
    #             return ""

    #         fig, ax = plt.subplots(figsize=(7, 6))
    #         sns.heatmap(
    #             arr,
    #             annot=True,
    #             fmt=".3f",
    #             cmap="RdBu_r",
    #             center=0,
    #             cbar_kws={"label": "Elemento da matriz (eV)"},
    #             ax=ax,
    #             square=True
    #         )

    #         n = arr.shape[0]
    #         labels = [str(i) for i in range(n)]
    #         ax.set_xticklabels(labels, rotation=45)
    #         ax.set_yticklabels(labels, rotation=0)
    #         ax.set_title(f"Matriz de Hückel")
    #         ax.set_xlabel("Índice")
    #         ax.set_ylabel("Índice")
    #         plt.tight_layout()

    #         filename = self.output_dir / f"{save_name.replace(' ', '_')}_huckel_matrix.png"
    #         plt.savefig(filename, dpi=300, bbox_inches="tight")
    #         plt.close()
    #         self.logger.info(f"Matriz de Hückel salva: {filename}")
    #         return str(filename)
    #     except Exception as e:
    #         self.logger.error(f"Erro em plot_huckel_matrix_heatmap: {e}")
    #         return ""

    def plot_huckel_matrix_heatmap(self, H_matrix, pi_atoms, save_name: str) -> str:
        try:
            if H_matrix is None:
                self.logger.warning("plot_huckel_matrix_heatmap: received None H_matrix, skipping.")
                return ""

            arr = np.asarray(H_matrix, dtype=float)

            # Basic validity checks
            if arr.ndim != 2 or arr.size == 0:
                self.logger.warning(f"plot_huckel_matrix_heatmap: invalid matrix shape={arr.shape}; skipping.")
                return ""

            n = arr.shape[0]

            annot = True if n <= 30 else False
            tick_target = 20
            tick_step = max(1, n // tick_target)

            max_abs = np.max(np.abs(arr)) if arr.size > 0 else 0.0
            rel_thresh = 0.03
            thresh = max_abs * rel_thresh
            mask = np.abs(arr) < thresh

            base = 6
            scale = min(0.30 * n, 18) 
            figsize = (base + scale, base + scale)

            fig, ax = plt.subplots(figsize=figsize)

            sns.heatmap(
                arr,
                mask=mask,
                annot=annot,
                fmt=".2f" if annot else None,
                cmap="RdBu_r",
                center=0,
                cbar=False,
                ax=ax,
                square=True,
                linewidths=0.2,
            )

            if n <= 40:
                xticks = np.arange(n)
                yticks = np.arange(n)
                xtick_labels = [str(i) for i in xticks]
                ytick_labels = [str(i) for i in yticks]
            else:
                ticks = np.arange(0, n, tick_step)
                xticks = ticks
                yticks = ticks
                xtick_labels = [str(int(t)) for t in ticks]
                ytick_labels = [str(int(t)) for t in ticks]

                ax.set_xticks(xticks + 0.5)
                ax.set_yticks(yticks + 0.5)

            ax.set_xticklabels(xtick_labels, rotation=45, ha="right", fontsize=8)
            ax.set_yticklabels(ytick_labels, rotation=0, fontsize=8)

            ax.set_title(f"Matriz de Hückel")
            ax.set_xlabel("Índice")
            ax.set_ylabel("Índice")
            plt.tight_layout()

            filename = self.output_dir / f"{save_name.replace(' ', '_')}_huckel_matrix.png"
            plt.savefig(filename, dpi=300, bbox_inches="tight")
            plt.close()
            self.logger.info(f"Hückel matrix heatmap saved: {filename}")
            return str(filename)
        except Exception as e:
            self.logger.error(f"Erro em plot_huckel_matrix_heatmap: {e}")
            return ""

# ============================================================================
# Módulo Gerador de Relatórios
# ============================================================================

class ReportGenerator:
    """Gera relatórios abrangentes dos cálculos."""
    
    def __init__(self, output_dir: str = "output"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.logger = logging.getLogger(f"{__name__}.ReportGenerator")
    
    def generate_html_report(self, results: Results, mol: Chem.Mol = None, images: Dict[str, str] = None) -> str:
        try:
            html_content = self._create_html_template(results, mol, images or {})
            method_slug = results.calculation_method.replace(' ', '_').replace('/', '-').replace('(', '').replace(')', '')
            filename = self.output_dir / f"{results.molecule_name}_{method_slug}_relatorio.html"
            with open(filename, 'w', encoding='utf-8') as f: f.write(html_content)
            self.logger.info(f"Relatório HTML gerado: {filename}")
            return str(filename)
        except Exception as e:
            self.logger.error(f"Erro ao gerar o relatório HTML: {e}"); return ""
    
    def _create_html_template(self, results: Results, mol: Chem.Mol = None, images: Dict[str, str] = None) -> str:
        energy_str = f"{results.energy:.6f} eV" if results.energy is not None else "N/D"
        homo_str = f"{results.homo_energy:.6f} eV" if results.homo_energy is not None else "N/D"
        lumo_str = f"{results.lumo_energy:.6f} eV" if results.lumo_energy is not None else "N/D"
        gap_str = f"{results.homo_lumo_gap:.6f} eV" if results.homo_lumo_gap is not None else "N/D"

        return f"""
        <!DOCTYPE html>
        <html lang="pt-br"><head><meta charset="UTF-8"><title>Relatório Hückel - {results.molecule_name}</title>
        <style>
            body {{ font-family: 'Segoe UI', sans-serif; line-height: 1.6; margin: 20px; background-color: #f5f5f5; }}
            .container {{ max-width: 1200px; margin: auto; background: white; padding: 30px; border-radius: 10px; box-shadow: 0 0 20px rgba(0,0,0,0.1); }}
            .header {{ text-align: center; color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 20px; margin-bottom: 30px; }}
            .section {{ margin-bottom: 30px; padding: 20px; border-left: 4px solid #3498db; background-color: #f8f9fa; text-align: center; }}
            .section h3 {{ color: #2c3e50; margin-top: 0; }}
            .property-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 20px; margin-top: 20px; }}
            .property-item {{ background: white; padding: 15px; border-radius: 5px; border-left: 3px solid #e74c3c; }}
            .image-container img {{ max-width: 100%; height: auto; border: 1px solid #ddd; border-radius: 5px; text-align: center; }}
            table {{ width: 100%; border-collapse: collapse; margin-top: 15px; margin: 0 auto; }}
            th, td {{ padding: 12px; text-align: left; border-bottom: 1px solid #ddd; }}
            th {{ background-color: #3498db; color: white; }}
            .status-planar {{ color: #27ae60; font-weight: bold; }} .status-quase-planar {{ color: #f39c12; font-weight: bold; }} .status-nao-planar {{ color: #e74c3c; font-weight: bold; }}
            .warning {{ border: 1px solid #e74c3c; background-color: #fbeae5; padding: 10px; border-radius: 5px; color: #c0392b; }}
        </style></head><body><div class="container">
            <div class="header"><h1>Relatório do Método de Hückel</h1><h2>{results.molecule_name}</h2><p><strong>SMILES:</strong> {results.input_smiles}</p></div>
            <div class="section"><h3>Resumo do Cálculo</h3><div class="property-grid">
                <div class="property-item"><strong>Energia dos Elétrons π:</strong><br>{energy_str}</div>
                <div class="property-item"><strong>Energia do HOMO:</strong><br>{homo_str}</div>
                <div class="property-item"><strong>Energia do LUMO:</strong><br>{lumo_str}</div>
                <div class="property-item"><strong>Diferença de Energia HOMO-LUMO:</strong><br>{gap_str}</div>
            </div></div>
            {self._generate_planarity_section(results)}
            {self._generate_image_section('Estrutura Molecular', images, 'structure')}
            {self._generate_image_section('Mapa de Ordem de Ligação π', images, 'bond_orders')}
            {self._generate_huckel_scheme_section(results, mol, images)}
            {self._generate_image_section('Matriz de Hückel', images, 'huckel_matrix')}
            {self._generate_image_section('Diagrama de Energia dos Orbitais', images, 'energy_diagram')}
            {self._generate_orbital_table(results)}
            {self._generate_population_table(results)}
            {self._generate_bond_order_table(results, mol)}
            {self._generate_frontier_maps_section(images)}
            {self._generate_molecular_properties_section(mol)}
        </div></body></html>"""

    def _generate_image_section(self, title, images, key):
        return f'<div class="section"><h3>{title}</h3><div class="image-container"><img src="{os.path.basename(images[key])}" alt="{title}"></div></div>' if key in images else ""

    def _generate_huckel_scheme_section(self, results, mol, images):
        if 'numbered_molecule' not in images: return ""
        return f"""<div class="section"><h3>Esquema do Sistema π de Hückel</h3>
            <div class="image-container"><img src="{os.path.basename(images['numbered_molecule'])}" alt="Molécula Numerada"></div>
            {self._generate_huckel_index_map_table(results, mol)}
            {self._generate_connectivity_table(results)}
        </div>"""

    def _generate_planarity_section(self, results: Results) -> str:
        if not results.planarity_status: return ""
        status_class = f"status-{results.planarity_status.replace('-', '_')}"
        warning_html = ""
        if results.planarity_status == "nao-planar":
            warning_html = '<p class="warning"><strong>Aviso:</strong> Sistema π não planar detectado. A teoria de Hückel assume planaridade, portanto estes resultados podem não ser confiáveis.</p>'
        elif results.planarity_status == "quase-planar":
            warning_html = '<p style="color: #d35400;"><strong>Nota:</strong> Desvios moderados da planaridade foram detectados.</p>'
        
        return f"""
        <div class="section">
            <h3>Verificação de Planaridade</h3>
            {warning_html}
            <div class="property-grid">
                <div class="property-item"><strong>Status:</strong><br><span class="{status_class}">{results.planarity_status.replace('-', ' ').title()}</span></div>
                <div class="property-item"><strong>Desvio RMS:</strong><br>{results.planarity_rms:.4f} Å</div>
                <div class="property-item"><strong>Desvio Máximo:</strong><br>{results.planarity_max_dev:.4f} Å</div>
                <div class="property-item"><strong>Método:</strong><br>{results.planarity_method.replace('_', ' ').title()}</div>
            </div>
        </div>"""

    def _generate_frontier_maps_section(self, images):
        if 'homo_map' not in images and 'lumo_map' not in images: return ""
        homo_html = f'<div class="image-container"><h4>Mapa do HOMO</h4><img src="{os.path.basename(images["homo_map"])}" alt="Mapa do HOMO"></div>' if 'homo_map' in images else ""
        lumo_html = f'<div class="image-container"><h4>Mapa do LUMO</h4><img src="{os.path.basename(images["lumo_map"])}" alt="Mapa do LUMO"></div>' if 'lumo_map' in images else ""
        return f'<div class="section"><h3>Mapas dos Orbitais Moleculares de Fronteira</h3><div class="property-grid">{homo_html}{lumo_html}</div></div>'

    def _generate_orbital_table(self, results: Results) -> str:
        if results.orbital_energies is None: return ""
        rows = ""
        for i, energy in enumerate(results.orbital_energies):
            orbital_type = "Ocupado"
            if results.somo_index is not None and i == results.somo_index: orbital_type = "SOMO"
            elif i >= results.num_occupied_orbitals: orbital_type = "Virtual"
            special = " (HOMO)" if i == results.homo_index else " (LUMO)" if i == results.lumo_index else " (SOMO)" if i == results.somo_index else ""
            rows += f"<tr><td>{i + 1}{special}</td><td>{energy:.6f}</td><td>{orbital_type}</td></tr>"
        return f'<div class="section"><h3>Energias dos Orbitais Moleculares</h3><table><thead><tr><th>Nº do OM</th><th>Energia (eV)</th><th>Tipo</th></tr></thead><tbody>{rows}</tbody></table></div>'

    def _generate_molecular_properties_section(self, mol: Chem.Mol) -> str:
        if not RDKIT_AVAILABLE or mol is None: return ""
        try:
            props = InputHandler().get_molecular_properties(mol)
            prop_items = "".join([f'<div class="property-item"><strong>{k.replace("_", " ").title()}:</strong><br>{v}</div>' for k, v in props.items()])
            return f'<div class="section"><h3>Propriedades Moleculares</h3><div class="property-grid">{prop_items}</div></div>'
        except: return ""
          
    def _generate_huckel_index_map_table(self, results: Results, mol: Chem.Mol) -> str:
        if results.pi_atoms is None: return ""
        rows = "".join([f"<tr><td>{h_idx}</td><td>{r_idx}</td><td>{mol.GetAtomWithIdx(r_idx).GetSymbol()}</td></tr>" for h_idx, r_idx in enumerate(results.pi_atoms)])
        return f'<h4>Mapeamento de Índices</h4><table><thead><tr><th>Índice Hückel</th><th>Índice RDKit</th><th>Átomo</th></tr></thead><tbody>{rows}</tbody></table>'

    def _generate_connectivity_table(self, results: Results) -> str:
        if results.connectivity is None: return ""
        rows = "".join([f"<tr><td>{idx}</td><td>{', '.join(map(str, sorted(neighbors)))}</td></tr>" for idx, neighbors in results.connectivity.items()])
        return f'<h4>Conectividade (Índices dos Átomos do Sistema π)</h4><table><thead><tr><th>Índice do Átomo</th><th>Vizinhos</th></tr></thead><tbody>{rows}</tbody></table>'

    def _generate_population_table(self, results: Results) -> str:
        if results.atom_populations is None: return ""
        rows = "".join([f"<tr><td>{i}</td><td>{pop:.4f}</td><td>{charge:.4f}</td></tr>" for i, (pop, charge) in enumerate(zip(results.atom_populations, results.pi_charges))])
        return f'<div class="section"><h3>Populações e Cargas de Elétrons π</h3><table><thead><tr><th>Índice do Átomo</th><th>População (qᵢ)</th><th>Carga π (zᵢ - qᵢ)</th></tr></thead><tbody>{rows}</tbody></table></div>'

    def _generate_bond_order_table(self, results: Results, mol: Chem.Mol) -> str:
        if results.bond_orders is None: return ""
        rows = "".join([f"<tr><td>p<sub>{i}-{j}</sub></td><td>{order:.4f}</td></tr>" for (i, j), order in sorted(results.bond_orders.items(), key=lambda item: item[1], reverse=True)])
        return f'<div class="section"><h3>Ordens de Ligação π (pᵢⱼ)</h3><table><thead><tr><th>Ligação</th><th>Ordem</th></tr></thead><tbody>{rows}</tbody></table></div>'

    def export_results_json(self, results: Results, filename: str = None) -> str:
        """Exporta os resultados para o formato JSON."""
        try:
            if filename is None:
                method_slug = results.calculation_method.replace(' ', '_').replace('/', '-')
                filename = self.output_dir / f"{results.molecule_name}_{method_slug}_resultados.json"
            else:
                filename = self.output_dir / filename
            
            results_dict = asdict(results)
            
            for key, value in results_dict.items():
                if isinstance(value, np.ndarray):
                    results_dict[key] = value.tolist()

            if results_dict.get('bond_orders'):
                results_dict['bond_orders'] = {f"{k[0]}-{k[1]}": v for k, v in results_dict['bond_orders'].items()}

            results_dict['export_timestamp'] = pd.Timestamp.now().isoformat()
            results_dict['pipeline_version'] = "1.2.0 (Completo)"
            
            with open(filename, 'w') as f:
                json.dump(results_dict, f, indent=2)
            
            self.logger.info(f"Resultados exportados para JSON: {filename}")
            return str(filename)
            
        except Exception as e:
            self.logger.error(f"Erro ao exportar para JSON: {e}")
            return ""

# ============================================================================
# Classe Principal do Pipeline
# ============================================================================

class QuantumChemistryPipeline:
    """Pipeline principal que orquestra todos os componentes."""
    
    def __init__(self, output_dir: str = "output"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.logger = logging.getLogger(f"{__name__}.Pipeline")
        self.input_handler = InputHandler() if RDKIT_AVAILABLE else None
        self.visualizer = Visualizer(output_dir)
        self.report_generator = ReportGenerator(output_dir)
        self.huckel_calc = HuckelCalculator()
        self.logger.info("Pipeline de Química Quântica inicializado")
    
    def run_single_calculation(self, mol_config: MoleculeConfig, huckel_config: HuckelConfig = None) -> Results:
        self.logger.info(f"Iniciando cálculo para {mol_config.name}")
        try:
            if self.input_handler is None: raise RuntimeError("RDKit não está disponível")
            mol = self.input_handler.smiles_to_mol(mol_config.smiles, mol_config.name)
            if huckel_config: self.huckel_calc.config = huckel_config
            results = self.huckel_calc.solve_huckel(mol, mol_config)
            
            if not results.success:
                self.logger.warning(f"Cálculo falhou para {mol_config.name}: {results.error_message}")
                return results

            images = {}
            
            if results.pi_atoms:
                planarity_data = self.visualizer.compute_pi_planarity(mol, results.pi_atoms)
                results.planarity_status = planarity_data['is_planar']
                results.planarity_rms = planarity_data['rms_dev']
                results.planarity_max_dev = planarity_data['max_dev']
                results.planarity_method = planarity_data['method']
                
                bond_order_file = self.visualizer.plot_bond_order_overlay(results, mol, mol_config.name)
                if bond_order_file: images['bond_orders'] = bond_order_file

            structure_file = self.visualizer.plot_molecular_structure(mol, mol_config.name)
            if structure_file: images['structure'] = structure_file

            # Add Hückel matrix heatmap only if valid
            hmat = getattr(results, "huckel_matrix", None)
            if isinstance(hmat, np.ndarray) and hmat.ndim == 2 and hmat.size > 0:
                hmat_file = self.visualizer.plot_huckel_matrix_heatmap(hmat, results.pi_atoms, mol_config.name)
                if hmat_file:
                    images['huckel_matrix'] = hmat_file
            else:
                self.logger.info("Skipping Hückel matrix heatmap: no valid 2-D matrix present.")
                    
            energy_diagram_file = self.visualizer.plot_orbital_energy_diagram(results, f"{mol_config.name}_diagrama_energia")
            if energy_diagram_file: images['energy_diagram'] = energy_diagram_file
            
            if results.pi_atoms:
                num_mol_file = self.visualizer.plot_numbered_molecule(mol, results.pi_atoms, mol_config.name)
                if num_mol_file: images['numbered_molecule'] = num_mol_file
                if results.homo_index is not None:
                    homo_map_file = self.visualizer.plot_huckel_mo_on_molecule(results, mol, results.homo_index, results.pi_atoms, f"{mol_config.name}_HOMO")
                    if homo_map_file: images['homo_map'] = homo_map_file
                if results.lumo_index is not None:
                    lumo_map_file = self.visualizer.plot_huckel_mo_on_molecule(results, mol, results.lumo_index, results.pi_atoms, f"{mol_config.name}_LUMO")
                    if lumo_map_file: images['lumo_map'] = lumo_map_file
            
            html_report = self.report_generator.generate_html_report(results, mol, images)
            json_export = self.report_generator.export_results_json(results)
            
            self.logger.info(f"Cálculo concluído com sucesso para {mol_config.name}")
            self.logger.info(f"Relatório HTML: {html_report}")
            self.logger.info(f"Exportação JSON: {json_export}")
            
            return results
            
        except Exception as e:
            error_msg = f"Falha no pipeline para {mol_config.name}: {str(e)}"
            self.logger.error(f"{error_msg}\n{traceback.format_exc()}")
            return Results(molecule_name=mol_config.name, input_smiles=mol_config.smiles, success=False, error_message=error_msg)
    
    def run_batch_calculations(self, molecules: List[MoleculeConfig], huckel_config: HuckelConfig = None) -> List[Results]:
        self.logger.info(f"Iniciando cálculos em lote para {len(molecules)} moléculas")
        results_list = [self.run_single_calculation(mol_config, huckel_config) for mol_config in molecules]
        self.logger.info(f"Cálculos em lote concluídos. {sum(1 for r in results_list if r.success)}/{len(results_list)} bem-sucedidos")
        return results_list

# ============================================================================
# Demo
# ============================================================================

def demo_pipeline():
    print("Demonstração do Pipeline do Método de Hückel (com Verificação de Planaridade e Sobreposição de Ordem de Ligação)")
    pipeline = QuantumChemistryPipeline("demo_output")
    test_molecules = [
        MoleculeConfig(smiles="C1=CC=CC=C1", name="benzeno"),
    ]
    results = pipeline.run_batch_calculations(test_molecules)
    print("\nDemonstração concluída com sucesso! Verifique o diretório 'demo_output'.")
    return results

if __name__ == "__main__":
    demo_pipeline()