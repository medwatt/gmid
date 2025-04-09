class TransistorSweep:
    def __init__(self, vgs, vds, vbs, length):
        """
        vgs        : tuple (start, stop, step) for gate-source voltage
        vds        : tuple (start, stop, step) for drain-source voltage
        vbs        : tuple (start, stop, step) for bulk-source voltage
        length     : list of channel lengths
        """
        self.vgs = vgs
        self.vds = vds
        self.vbs = vbs
        self.length = length
