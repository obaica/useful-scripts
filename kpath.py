import seekpath
import numpy as np

element_to_atomic_number = {
    'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10,
    'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18,
    'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26,
    'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34,
    'Br': 35, 'Kr': 36, 'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42,
    'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50,
    'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56, 'La': 57, 'Ce': 58,
    'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66,
    'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72, 'Ta': 73, 'W': 74,
    'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82,
    'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90,
    'Pa': 91, 'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98,
    'Es': 99, 'Fm': 100, 'Md': 101, 'No': 102, 'Lr': 103, 'Rf': 104, 'Db': 105,
    'Sg': 106, 'Bh': 107, 'Hs': 108, 'Mt': 109, 'Ds': 110, 'Rg': 111, 'Cn': 112,
    'Nh': 113, 'Fl': 114, 'Mc': 115, 'Lv': 116, 'Ts': 117, 'Og': 118
}

def read_poscar(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    # 读取晶胞参数
    cell = np.zeros((3, 3))
    for i in range(3):
        cell[i] = [float(x) for x in lines[2 + i].split()]

    # 读取原子类型和数量
    species = lines[5].split()
    num_atoms = [int(x) for x in lines[6].split()]
    atomic_numbers = []
    for element, count in zip(species, num_atoms):
        atomic_numbers.extend([element_to_atomic_number[element]] * count)

    positions = []
    start_line = 8
    for i, num in enumerate(num_atoms):
        for j in range(num):
            pos = [float(x) for x in lines[start_line + sum(num_atoms[:i]) + j].split()[:3]]  # 只读取前三个坐标值！！！
            positions.append(pos)
    return cell, positions, atomic_numbers
    
# 从 POSCAR 文件读取结构
cell, positions, atomic_numbers = read_poscar('POSCAR')

# 转换为 seekpath 可接受的格式
structure = (cell, positions, atomic_numbers)

# 获取高对称性路径
path_data = seekpath.get_path(structure)

# 提取高对称点和建议路径
kpoints = path_data['point_coords']
band_path = path_data['path']

# 写入 KPOINTS 文件
with open('KPOINTS', 'w') as file:
    file.write('K-Path\n')
    file.write('100\n')  # 高对称点之间的点数（可调整）
    file.write('Reciprocal\n')
    for i, segment in enumerate(band_path):
        start, end = segment
        start_coords = kpoints[start]
        end_coords = kpoints[end]
        file.write(f"{np.array(start_coords)[0]:.4f} {np.array(start_coords)[1]:.4f} {np.array(start_coords)[2]:.4f} {start}\n")
        file.write(f"{np.array(end_coords)[0]:.4f} {np.array(end_coords)[1]:.4f} {np.array(end_coords)[2]:.4f} {end}\n")
        if i < len(band_path) - 1:
            file.write("\n")  # 在每两个高对称点之间添加一个空行

print("KPOINTS file generated for VASP energy band calculation.")
