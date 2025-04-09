import numpy as np

class MosfetSimulation:
    def __init__(self, simulator, netlist_gen, sweeps):
        """
        sweeps: dict mapping transistor name to TransistorSweep objects.
        """
        self.simulator = simulator
        self.netlist_gen = netlist_gen
        self.sweeps = sweeps
        self.lookup_table = {}

        # Initialize lookup table for each transistor model using each sweep's parameters.
        for mosfet_model, sweep in self.sweeps.items():
            n_vgs = int(round((sweep.vgs[1] - sweep.vgs[0]) / sweep.vgs[2])) + 1
            n_vds = int(round((sweep.vds[1] - sweep.vds[0]) / sweep.vds[2])) + 1
            n_vbs = int(round((sweep.vbs[1] - sweep.vbs[0]) / sweep.vbs[2])) + 1
            self.lookup_table[mosfet_model] = {
                p: np.zeros((len(sweep.length), n_vbs, n_vgs, n_vds), dtype=np.float32)
                for p in simulator.parameters_to_save
            }

    def linspace(self, start, stop, step):
        count = int(round((stop - start) / step)) + 1
        return np.linspace(start, stop, count)

    def get_single_netlist(self, sim_func):
        # Select the first transistor sweep for a sample netlist.
        transistor_name, sweep = next(iter(self.sweeps.items()))
        length = sweep.length[0]
        vbs_val = next(iter(self.linspace(*sweep.vbs)))
        netlist = self.netlist_gen.generate_netlist(transistor_name, length, vbs_val)
        sim_setup = sim_func(sweep.vgs, sweep.vds)
        return netlist + sim_setup

    def print_netlist(self):
        netlist = self.get_single_netlist(self.simulator.setup_dc_simulation)
        print("---------------------------------------------------")
        print("----- This is the netlist that gets simulated -----")
        print("---------------------------------------------------")
        print("\n".join(netlist))
        print("---------------------------------------------------")
        print("")

    def op_simulation(self):
        netlist = self.get_single_netlist(self.simulator.setup_op_simulation)
        self.simulator.run_simulation(netlist, verbose=True)

    def simulate(self):
        self.print_netlist()
        for transistor_name, sweep in self.sweeps.items():
            print(f"Simulating {transistor_name}")
            for l_idx, length in enumerate(sweep.length):
                print("Length:", length)
                for vbs_idx, vbs_val in enumerate(self.linspace(*sweep.vbs)):
                    netlist = self.netlist_gen.generate_netlist(transistor_name, length, vbs_val)
                    sim_setup = self.simulator.setup_dc_simulation(sweep.vgs, sweep.vds)
                    full_netlist = netlist + sim_setup
                    self.simulator.run_simulation(full_netlist)
                    analysis = self.simulator.parse_output()
                    self.simulator.save_parameters(
                        analysis, transistor_name, l_idx, vbs_idx,
                        self.lookup_table,
                        int(round((sweep.vgs[1]-sweep.vgs[0])/sweep.vgs[2]))+1,
                        int(round((sweep.vds[1]-sweep.vds[0])/sweep.vds[2]))+1
                    )
            print("")
