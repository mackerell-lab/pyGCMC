# examples/basic/python_example.py

import pyGCMC

def main():
    # 创建一个系统实例，使用默认的 epsilon=1.0 和 sigma=1.0
    system = pyGCMC.System()

    # 从文件加载
    pdb_file = "examples/data/example.pdb"          # 确保存在该文件
    psf_file = "examples/data/example.psf"          # 确保存在该文件
    top_file = "examples/data/example.top"          # 确保存在该文件
    ff_file = "examples/data/charmm36.ff"           # 确保存在该文件

    system.load_pdb(pdb_file)
    system.load_psf(psf_file)
    system.load_top(top_file)
    system.load_forcefield(ff_file)

    # 输出粒子数量
    print(f"Particle count after loading files: {system.get_particle_count()}")

    # 计算并输出总能量
    energy = system.compute_total_energy()
    print(f"Total energy: {energy}")

if __name__ == "__main__":
    main()
