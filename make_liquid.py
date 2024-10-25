import os
import subprocess
import parmed as pmd
from openbabel import pybel
import re
import time

def convert_mol2_to_pdb(input_mol2, output_pdb):
    """
    Open Babelを使用してmol2ファイルをpdbファイルに変換
    """
    mol = next(pybel.readfile("mol2", input_mol2))
    mol.write("pdb", output_pdb, overwrite=True)
    print(f"{input_mol2} を {output_pdb} に変換しました。")

def calculate_molecular_weight(input_mol2):
    """
    mol2ファイルから分子量を計算
    """
    mol = next(pybel.readfile("mol2", input_mol2))
    molecular_weight = mol.molwt
    print(f"分子量は {molecular_weight:.2f} g/mol です。")
    return molecular_weight

def run_packmol(input_pdb, output_pdb, num_molecules, box_size):
    """
    Packmolを使用して複数分子系を作成
    """
    packmol_input = f"""
    seed 1234567
    tolerance 2.0
    filetype pdb
    output {output_pdb}
    structure {input_pdb}
      number {num_molecules}
      inside box  -{box_size/2}  -{box_size/2}  -{box_size/2}   {box_size/2}   {box_size/2}   {box_size/2}
    end structure
    """
    
    with open("packmol_input.inp", "w") as f:
        f.write(packmol_input)

    result = subprocess.run("packmol < packmol_input.inp", shell=True, capture_output=True, text=False)
    if result.returncode != 0:
        print(f"Packmolでエラー: {result.stderr}")
        raise RuntimeError("Packmolの実行に失敗しました。")
    else:
        print(f"Packmol実行中..")

    print(f"Packmolを使用して複数分子系を作成しました: {output_pdb}")

def extract_charges_from_mol2(input_mol2):
    """
    mol2ファイルから電荷情報を抽出
    """
    mol = pmd.load_file(input_mol2)
    charges = [atom.charge for atom in mol.atoms]
    return charges

def assign_charges_to_packed_system(input_mol2, packed_pdb, output_mol2):
    """
    元のmol2ファイルから電荷情報を取得し、packmolで作成した複数分子系に適用
    """
    # 元のmol2ファイルから電荷を取得
    charges = extract_charges_from_mol2(input_mol2)

    # 複数分子系の構造をロード
    packed_structure = pmd.load_file(packed_pdb, structure=True)
    print(packed_structure)
    # 元の分子数を取得
    num_original_atoms = len(charges)

    # 各分子に電荷を割り当て
    for i, atom in enumerate(packed_structure.atoms):
        atom.charge = charges[i % num_original_atoms]

    # 既存のoutput_mol2ファイルがある場合は削除
    if os.path.exists(output_mol2):
        os.remove(output_mol2)

    # 電荷が割り当てられたmol2ファイルを出力
    packed_structure.save(output_mol2, format='mol2')
    print(f"電荷を反映した複数分子系のmol2ファイルを作成しました: {output_mol2}")
    
def extract_atom_types_and_bonds_from_mol2(mol2_file):
    """
    mol2ファイルのATOMセクションから2列目（原子名）と6列目（原子タイプ）、
    およびBONDセクションから結合情報を抽出。
    """
    atom_types = []
    bond_info = []
    atom_section = False
    bond_section = False

    with open(mol2_file, 'r') as file:
        lines = file.readlines()
        for line in lines:
            # ATOMセクションの検出
            if line.startswith("@<TRIPOS>ATOM"):
                atom_section = True
                bond_section = False
                continue
            # BONDセクションの検出
            if line.startswith("@<TRIPOS>BOND"):
                atom_section = False
                bond_section = True
                continue
            # セクションの終了
            if line.startswith("@<TRIPOS>"):
                atom_section = False
                bond_section = False
            
            # ATOMセクションの処理
            if atom_section:
                parts = line.split()
                atom_type = parts[5]  # 6列目が原子タイプ
                atom_types.append(atom_type)
            # BONDセクションの処理
            if bond_section:
                parts = line.split()
                bond_info.append((int(parts[1]), int(parts[2]), parts[3]))  # 結合先と結合次数
    
    return atom_types, bond_info


def update_mol2_with_atom_types_and_bonds(packed_mol2, atom_types, bond_info, num_molecules, output_mol2):
    """
    packed_molecules_with_charges.mol2 の6列目を元のmol2の原子タイプに基づいて全分子に適用し、
    BOND セクションを単分子のmol2ファイルに基づいて修正し、SUBSTRUCTURE情報は残す。
    """
    with open(packed_mol2, 'r') as file:
        lines = file.readlines()

    with open(output_mol2, 'w') as file:
        atom_section = False
        bond_section = False
        atom_index = 0
        molecule_count = 0
        bond_index = 1
        current_molecule = 0  # 現在処理している分子番号を追跡

        for line in lines:
            if line.startswith("@<TRIPOS>ATOM"):
                atom_section = True
                file.write(line)
                continue
            if line.startswith("@<TRIPOS>BOND"):
                atom_section = False
                bond_section = True
                # 元のBONDセクションをスキップ
                continue
            if line.startswith("@<TRIPOS>SUBSTRUCTURE"):
                # 新しいBONDセクションヘッダーを追加
                file.write("@<TRIPOS>BOND\n")
                # 新しいBOND情報を追加
                for mol in range(num_molecules):
                    for bond in bond_info:
                        atom1 = bond[0] + mol * len(atom_types)
                        atom2 = bond[1] + mol * len(atom_types)
                        bond_order = bond[2]
                        new_bond_line = f"{bond_index} {atom1} {atom2} {bond_order}\n"
                        file.write(new_bond_line)
                        bond_index += 1
                # SUBSTRUCTUREセクションの書き込み
                bond_section = False
                file.write(line)
                continue
            if line.startswith("@<TRIPOS>"):
                atom_section = False
                bond_section = False
                file.write(line)
                continue
            
            # ATOMセクションの処理
            if atom_section:
                parts = line.split()
                # 元のmol2ファイルの原子タイプを繰り返し適用
                parts[5] = atom_types[atom_index % len(atom_types)]
                new_line = " ".join(parts) + "\n"
                file.write(new_line)
                atom_index += 1
                # 次の分子に移るたびにカウント
                if atom_index == len(atom_types) * (molecule_count + 1):
                    molecule_count += 1
                    atom_index = 0  # リセットして次の分子に
            elif not bond_section:
                file.write(line)
    
def run_tleap(input_mol2, prmtop_file, inpcrd_file, clearance):
    """
    tleapを使用してAmber形式のトポロジーと座標ファイルを生成。クリアランスを指定して周期境界条件を設定。
    """
    tleap_input = f"""
    source leaprc.gaff2
    
    mol = loadmol2 {input_mol2}
    setBox mol vdw {clearance}
    saveamberparm mol {prmtop_file} {inpcrd_file}
    quit
    """

    with open('leap.in', 'w') as f:
        f.write(tleap_input)

    tleap_cmd = ['tleap', '-f', 'leap.in']
    result = subprocess.run(tleap_cmd, capture_output=True, text=True)
    
    # エラー時に詳細メッセージを表示
    if result.returncode != 0:
        print(f"tleapでエラー: {result.stderr}")
        print(f"tleapの標準出力: {result.stdout}")
        raise RuntimeError("tleapの実行に失敗しました。")

    print(f"tleapでトポロジーファイルと座標ファイルを生成しました。")

    

def convert_to_gromacs(prmtop_file, inpcrd_file):
    """
    ParmEdを使用して、Amberのトポロジーファイルと座標ファイルをGROMACSのデータファイルに変換。
    """
    structure = pmd.load_file(prmtop_file, inpcrd_file)

    # 既存のtemp.groファイルがある場合は削除
    if os.path.exists('temp.gro'):
        os.remove('temp.gro')
    # 既存のtemp.topファイルがある場合は削除
    if os.path.exists('temp.top'):
        os.remove('temp.top')
    if os.path.exists('liq.mol2'):
        os.remove('liq.mol2')

    # GROMACS形式のファイルを保存
    structure.save('temp.gro', format='gro')
    structure.save('temp.top', format='gromacs')
    structure.save('liq.mol2', format='mol2')
    
    print("GROMACSのデータファイルを生成しました")


def main():
    # カレントディレクトリ内のmol_で始まるすべてのディレクトリを探索
    for dir_name in os.listdir('.'):
        if os.path.isdir(dir_name) and dir_name.startswith("mol_"):
            os.chdir(dir_name)  # 各ディレクトリに移動
            
            # ディレクトリ名からmol2ファイル名を決定
            input_mol2 = f"{dir_name}_.mol2"  # ディレクトリ名と同じ名前のmol2ファイル
            amber_input_mol2 = f"{dir_name}.mol2"
            single_pdb = f"{dir_name}.pdb"     # 一時的なpdbファイル
            packed_pdb = f"packed_{dir_name}.pdb"    # packmolで作成された複数分子系のpdbファイル
            output_mol2 = f"packed_{dir_name}_with_charges.mol2"  # 最終的な複数分子系のmol2ファイル
            num_molecules = 200                    # 分子数
            density = 0.6                          # 密度（g/cm^3）
            acpype_output_prefix = f"system_{dir_name}_gaff2"  # Acpypeの出力ファイル用の接頭辞
            prmtop_file = f"{acpype_output_prefix}.prmtop"
            inpcrd_file = f"{acpype_output_prefix}.inpcrd"

            # mol2からpdbへの変換
            convert_mol2_to_pdb(input_mol2, single_pdb)

            # mol2ファイルから分子量を計算
            mol_weight = calculate_molecular_weight(input_mol2)

            # 密度からボックスサイズを計算
            avogadro = 6.022e23  # アボガドロ数
            volume_per_molecule = mol_weight / (density * avogadro) * 1e24  # Å^3での1分子当たりの体積
            box_size = (volume_per_molecule * num_molecules) ** (1/3)
            print(f"計算されたボックスサイズ: {box_size:.2f} Å")

            # Packmolで複数分子系を作成
            run_packmol(single_pdb, packed_pdb, num_molecules, box_size)

            # 電荷を反映してmol2ファイルを作成
            assign_charges_to_packed_system(amber_input_mol2, packed_pdb, output_mol2)
            
            # 元のmol2ファイルから2列目と6列目の情報を抽出
            atom_types, bond_info = extract_atom_types_and_bonds_from_mol2(amber_input_mol2)

            # packed_molecules_with_charges.mol2の6列目を更新
            update_mol2_with_atom_types_and_bonds(output_mol2, atom_types, bond_info, num_molecules, output_mol2)

            # tleapを使用してAmber形式のトポロジーと座標ファイルを生成
            run_tleap(output_mol2, prmtop_file, inpcrd_file, 0)

            # ParmEdを使ってGROMACSのファイルを生成
            convert_to_gromacs(prmtop_file, inpcrd_file)
            
            # ディレクトリを元に戻す
            os.chdir("..")


if __name__ == "__main__":
    main()
