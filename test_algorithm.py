from Solver import Solver

dementia_ids = ["C0524851", "C0002736", "C0030567",
                "C0002395", "C0020179", "C0497327", 
                "C0338656", "C0011265", "C0038454",
                "C0242422", "C0233794", "C0026769",
                "C0338451", "C0752347"]

osa_ids = [ "C0520679", "C0037315", "C0520680"]


A = Solver.get_neighbors(dementia_ids)
B = Solver.get_neighbors(osa_ids)

Solver.solve_iterative_greedy(A, B)