#!/usr/bin/env python3
import numpy as np
from mosplot.plot.mosplot import Mosfet
from mosplot.plot.helpers import load_lookup_table

def db_to_gain(db):
    return 10**(db/20)

def get_largest_closest_match(target, arr):
    return arr[np.argmin(np.abs(arr - target))]

class OTA5TDesigner:
    def __init__(self, *, table, nmos_model, pmos_model, VDD, DC_GAIN, GBW, CL, IBIAS, ICMR_LOW, ICMR_HIGH, CMRR):
        self.VDD = VDD                        # DC supply voltage
        self.DC_GAIN = db_to_gain(DC_GAIN)    # DC gain specification
        self.CMRR = db_to_gain(CMRR)          # Common-Mode Rejection Ratio
        self.GBW = GBW                        # Gain-bandwidth specification (Hz)
        self.CL = CL                          # Load capacitance (F)
        self.IBIAS = IBIAS                    # Bias current (A)
        self.ICMR_LOW = ICMR_LOW              # ICMR Low (V)
        self.ICMR_HIGH = ICMR_HIGH            # ICMR High (V)

        self.nmos_model = nmos_model
        self.pmos_model = pmos_model
        self.lookup_table = load_lookup_table(table)

        # Designed transistor parameters will be stored as dictionaries
        self.input_pair = None
        self.current_mirror_load = None
        self.tail_current_source = None

        # Range of voltages to use
        offset = 5
        self.vgs_range_pmos = (self.lookup_table[self.pmos_model]["vgs"][-offset], self.lookup_table[self.pmos_model]["vgs"][offset])
        self.vgs_range_nmos = (self.lookup_table[self.nmos_model]["vgs"][offset], self.lookup_table[self.nmos_model]["vgs"][-offset])

    def design_input_pair(self):
        pmos = Mosfet(
            lookup_table=self.lookup_table,
            mos=self.pmos_model,
            vsb=0.0,
            vds=-self.VDD/3,
            vgs=self.vgs_range_pmos,
        )

        # The current of the branch is from the specs.
        current = self.IBIAS

        # Get the required gm of the input pair from the GBW and CL specs.
        gm = 2 * np.pi * self.CL * self.GBW

        # Get the required gmid from the gm and IBIAS specs.
        gmid = np.ceil(gm / self.IBIAS)

        # DC_GAIN is given by gm1 * rout.
        # rout = rds2 || rds1 => gds_out = gds1 + gds2
        # Thus: gds1 + gds2 = gm1 / DC_GAIN

        # gds1,2 is given by:
        gds = gm / (2 * self.DC_GAIN)

        # The intrinsic gain of M1 and M2
        gain = 2 * self.DC_GAIN * 1.1

        # From the desired gmid and intrinsic gain, we can get the length.
        # Since gain depends on the output resistance, which in turn depends on the
        # length, the bigger the required gain, the longer the transistor.
        design_length = pmos.interpolate(
            x_expression=pmos.gmid_expression,
            x_value=gmid,
            y_expression=pmos.gain_expression,
            y_value=gain,
            z_expression=pmos.lengths_expression,
        )
        if not np.isnan(design_length):
            design_length = get_largest_closest_match(design_length, pmos.lengths)
        else:
            raise ValueError("The transistors do not have enough intrinsic gain to meet gain requirements.")

        # We now have everything. We simply commpute the rest of the required variables.
        current_density = pmos.interpolate(
            x_expression=pmos.gmid_expression,
            x_value=gmid,
            y_expression=pmos.lengths_expression,
            y_value=design_length,
            z_expression=pmos.current_density_expression,
        )
        design_width = current / current_density

        vgs = pmos.interpolate(
            x_expression=pmos.gmid_expression,
            x_value=gmid,
            y_expression=pmos.lengths_expression,
            y_value=design_length,
            z_expression=pmos.vgs_expression,
        )

        vdsat = pmos.interpolate(
            x_expression=pmos.gmid_expression,
            x_value=gmid,
            y_expression=pmos.lengths_expression,
            y_value=design_length,
            z_expression=pmos.vdsat_expression,
        )

        self.input_pair = {
            "current": current,
            "gmid": gmid,
            "gain": gain,
            "length": design_length,
            "width": design_width,
            "vsg": np.abs(vgs),
            "vdsat": vdsat,
            "gm": gm,
            "gds": gds
        }

    def design_current_mirror_load(self):
        nmos = Mosfet(
            lookup_table=self.lookup_table,
            mos=self.nmos_model,
            vsb=0.0,
            vds=self.VDD/3,
            vgs=self.vgs_range_nmos,
        )

        # The current of the branch is from the specs.
        current = self.IBIAS

        # The intrinsic gain of M1 and M2:
        gain = 2 * self.DC_GAIN

        # Since we don't know the gmid of this transistor, we can
        # sweep through various lengths until we get one that meets
        # the gain requirement.
        for design_length in nmos.lengths:
            gmid_for_gain = nmos.interpolate(
                x_expression=nmos.lengths_expression,
                x_value=design_length,
                y_expression=nmos.gain_expression,
                y_value=gain,
                z_expression=nmos.gmid_expression,
            )
            if not np.isnan(gmid_for_gain) and gmid_for_gain > 1 and gmid_for_gain < 20:
                break
        else:
            raise ValueError("Current mirror load transistors cannot meet gain requirement")

        # Use estimated length and ICMR_LOW spec to estimate required gmid.
        vgs3max = self.ICMR_LOW + self.input_pair["vsg"] - self.input_pair["vdsat"]
        gmid_for_icmr = nmos.interpolate(
            x_expression=nmos.lengths_expression,
            x_value=design_length,
            y_expression=nmos.vgs_expression,
            y_value=vgs3max[0],
            z_expression=nmos.gmid_expression,
        )

        if gmid_for_icmr < gmid_for_gain:
            # To bias deeper into saturation.
            gmid = np.ceil(gmid_for_icmr) + 1
        else:
            raise ValueError("Current mirror load transistors cannot satisfy gain and ICMR_LOW requirements")

        # We now have everything. We simply commpute the rest of the required variables.
        current_density = nmos.interpolate(
            x_expression=nmos.gmid_expression,
            x_value=gmid,
            y_expression=nmos.lengths_expression,
            y_value=design_length,
            z_expression=nmos.current_density_expression,
        )
        design_width = current / current_density

        vdsat = nmos.interpolate(
            x_expression=nmos.gmid_expression,
            x_value=gmid,
            y_expression=nmos.lengths_expression,
            y_value=design_length,
            z_expression=nmos.vdsat_expression,
        )
        design_width = current / current_density

        # The gm of this transistor is required in the CMRR calculation
        # for the tail current source.
        gm = gmid * current

        self.current_mirror_load = {
            "current": current,
            "gmid": gmid,
            "vdsat": vdsat,
            "length": design_length,
            "width": design_width,
            "vgs": vgs3max,
            "gm": gm,
        }

    def design_tail_current_source(self):
        pmos = Mosfet(
            lookup_table=self.lookup_table,
            mos=self.pmos_model,
            vsb=0.0,
            vds=-self.VDD/3,
            vgs=self.vgs_range_pmos,
        )

        # The tail current source supplies both branches.
        current = 2 * self.IBIAS

        # Since we don't know the gmid of this transistor, we can
        # sweep through the various lengths, looking for the one that
        # satisfies the ICMR_HIGH requirement.
        is_icmr_satisfied = False
        vdsat_max = self.VDD - self.input_pair["vsg"][0] - self.ICMR_HIGH
        for design_length in pmos.lengths:
            gmid = pmos.interpolate(
                x_expression=pmos.lengths_expression,
                x_value=design_length,
                y_expression=pmos.vdsat_expression,
                y_value=vdsat_max,
                z_expression=pmos.gmid_expression,
            )
            if not np.isnan(gmid) and gmid > 1 and gmid < 15:
               is_icmr_satisfied = True
               # We also check if the computed gmid meets the CMRR requirement.
               gm = current * gmid
               gds = 2 * self.DC_GAIN * self.current_mirror_load["gm"]  / self.CMRR
               desired_gain =  (gm / gds) * 1.1

               gain = pmos.interpolate(
                   x_expression=pmos.gmid_expression,
                   x_value=gmid[0],
                   y_expression=pmos.lengths_expression,
                   y_value=design_length,
                   z_expression=pmos.gain_expression,
               )

               if not np.isnan(gain) and gain > desired_gain:
                   break
        else:
            if is_icmr_satisfied:
                raise ValueError("Tail current transistor cannot CMRR requirements.")
            else:
                raise ValueError("Tail current transistor cannot meet ICMR_HIGH requirements.")

        # Bias transistor deeper into saturation
        gmid = np.ceil(gmid) + 1

        # We now have everything. We simply commpute the rest of the required variables.
        current_density = pmos.interpolate(
            x_expression=pmos.gmid_expression,
            x_value=gmid,
            y_expression=pmos.lengths_expression,
            y_value=design_length,
            z_expression=pmos.current_density_expression,
        )
        design_width = current / current_density

        vgs = pmos.interpolate(
            x_expression=pmos.gmid_expression,
            x_value=gmid,
            y_expression=pmos.lengths_expression,
            y_value=design_length,
            z_expression=pmos.vgs_expression,
        )

        vdsat = pmos.interpolate(
            x_expression=pmos.gmid_expression,
            x_value=gmid,
            y_expression=pmos.lengths_expression,
            y_value=design_length,
            z_expression=pmos.vdsat_expression,
        )

        self.tail_current_source = {
            "current": current,
            "gmid": gmid,
            "length": design_length,
            "width": design_width,
            "vgs": vgs,
            "vdsat": vdsat
        }

    def run_design(self):
        self.design_input_pair()
        self.design_current_mirror_load()
        self.design_tail_current_source()

    def print_design(self):
        print("5T OTA Design Parameters:")
        print("\nInput Pair (PMOS M1/M2):")
        print(f"  gm/ID       : {self.input_pair['gmid']:.2f} S/A")
        print(f"  I           : {self.input_pair['current']*1e6:.2f} μA")
        print(f"  L           : {self.input_pair['length']:.3e} m")
        print(f"  W           : {self.input_pair['width'][0]:.3e} m")
        print(f"  V_SG        : {self.input_pair['vsg'][0]:.3e} V")
        print(f"  V_DS_SAT    : {self.input_pair['vdsat'][0]:.3e} V")

        print("\nCurrent Mirror Load (NMOS M3/M4):")
        print(f"  gm/ID       : {self.current_mirror_load['gmid'][0]:.2f} S/A")
        print(f"  I           : {self.current_mirror_load['current']*1e6:.2f} μA")
        print(f"  L           : {self.current_mirror_load['length']:.3e} m")
        print(f"  W           : {self.current_mirror_load['width'][0]:.3e} m")
        print(f"  V_GS        : {self.current_mirror_load['vgs'][0]:.3e} V")
        print(f"  V_DS_SAT    : {self.current_mirror_load['vdsat'][0]:.3e} V")

        print("\nTail Current Source (PMOS M5):")
        print(f"  gm/ID         : {self.tail_current_source['gmid'][0]:.2f} S/A")
        print(f"  I             : {self.tail_current_source['current']*1e6:.2f} μA")
        print(f"  L             : {self.tail_current_source['length']:.3e} m")
        print(f"  W             : {self.tail_current_source['width'][0]:.3e} m")
        print(f"  V_SG          : {abs(self.tail_current_source['vgs'][0]):.3e} V")
        print(f"  V_DS_SAT      : {self.tail_current_source['vdsat'][0]:.3f} V")

        print("\nBias Current Mirror Transistor (PMOS M6):")
        print(f"  I             : {0.5*self.tail_current_source['current']*1e6:.2f} μA")
        print(f"  L             : {self.tail_current_source['length']:.3e} m")
        print(f"  W             : {0.5*self.tail_current_source['width'][0]:.3e} m")

if __name__ == "__main__":
    path = "/home/medwatt/git/gmid/tests/lookup_table_generator/tsmc65/tsmc_65nm.npy"
    designer = OTA5TDesigner(
        table=path,
        nmos_model="nch_lvt",
        pmos_model="pch_lvt",
        VDD=1.2,
        DC_GAIN=25,
        GBW=5e6,
        CL=5e-12,
        IBIAS=10e-6,
        ICMR_LOW=0.15,
        ICMR_HIGH=0.65,
        CMRR=50,
    )
    designer.run_design()
    designer.print_design()

