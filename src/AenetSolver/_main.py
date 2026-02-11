import odatse

from .aenet import Solver


def main():
    """Entry point for the AENET solver with ODAT-SE.

    Initializes ODAT-SE, creates an :class:`AenetSolver.aenet.Solver`,
    wraps it with :class:`odatse.Runner`, and runs the selected algorithm.
    """
    info, run_mode = odatse.initialize()

    alg_module = odatse.algorithm.choose_algorithm(info.algorithm["name"])

    solver = Solver(info)
    runner = odatse.Runner(solver, info)
    alg = alg_module.Algorithm(info, runner, run_mode=run_mode)

    alg.main()
