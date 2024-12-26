# examples/basic/example.py
import pyGCMC

def main():
    # 创建一个系统实例，使用默认的 epsilon=1.0 和 sigma=1.0
    system = pyGCMC.System()

    # 创建两个粒子
    p1 = pyGCMC.Particle(x=0.0, y=0.0, z=0.0, type=1, charge=0.5)
    p2 = pyGCMC.Particle(x=3.355, y=0.0, z=0.0, type=2, charge=-0.5)

    # 添加粒子到系统中
    system.add_particle(p1)
    system.add_particle(p2)

    # 输出粒子数量
    print(f"Particle count: {system.get_particle_count()}")

    # 计算并输出总能量
    energy = system.compute_total_energy()
    print(f"Total energy: {energy}")

if __name__ == "__main__":
    main()
