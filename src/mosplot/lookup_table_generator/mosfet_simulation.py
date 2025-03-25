# imports <<<
import numpy as np
# >>>

class MosfetSimulation:
    def __init__(self, simulator, netlist_gen, vgs, vds, vsb, lengths, n_vgs, n_vds, n_vsb, model_names):
        self.simulator = simulator
        self.netlist_gen = netlist_gen
        self.vgs = vgs
        self.vds = vds
        self.vsb = vsb
        self.lengths = lengths
        self.n_vgs = n_vgs
        self.n_vds = n_vds
        self.n_vsb = n_vsb
        self.model_names = model_names
        self.lookup_table = {}
        for mosfet_model in self.model_names.keys():
            self.lookup_table[mosfet_model] = {
                p: np.zeros((len(self.lengths), n_vsb, n_vgs, n_vds))
                for p in simulator.parameters_to_save
            }

    def linspace(self, start, stop, step):
        count = int(round((stop - start) / step)) + 1
        return np.linspace(start, stop, count)

    def get_single_netlist(self, sim_func):
        transistor_name = next(iter(self.model_names.keys()))
        length = next(iter(self.lengths))
        vsb_val = next(iter(self.linspace(*self.vsb)))
        netlist = self.netlist_gen.generate_netlist(transistor_name, length, vsb_val)
        sim_setup = sim_func(self.vgs, self.vds)
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
        for transistor_name, type_type in self.model_names.items():
            r = -1 if type_type == "pmos" else 1

            print(f"Simulating {transistor_name}")

            for l_idx, length in enumerate(self.lengths):

                print("Length:", length)

                for vsb_idx, vsb_val in enumerate(self.linspace(*self.vsb)):
                    netlist = self.netlist_gen.generate_netlist(transistor_name, length, vsb_val)
                    sim_setup = self.simulator.setup_dc_simulation(self.vgs * r, self.vds * r)
                    full_netlist = netlist + sim_setup
                    self.simulator.run_simulation(full_netlist)
                    analysis = self.simulator.parse_output()
                    self.simulator.save_parameters(
                        analysis, transistor_name, l_idx, vsb_idx,
                        self.lookup_table, self.n_vgs, self.n_vds
                    )

            print("")
