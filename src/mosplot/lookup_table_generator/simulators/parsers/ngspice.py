import numpy as np

class NgspiceRawFileReader:
    def __init__(self):
        self.bsize_sp = 512
        self.mdata_list = [
            b"title",
            b"date",
            b"plotname",
            b"flags",
            b"no. variables",
            b"no. points",
            b"dimensions",
            b"command",
            b"option",
        ]

    def read_file(self, fname):
        """Read ngspice binary raw files. Return tuple of the data, and the
        plot metadata. The dtype of the data contains field names. This is
        not very robust yet, and only supports ngspice.
        >>> darr, mdata = rawread('test.py')
        >>> darr.dtype.names
        >>> plot(np.real(darr['frequency']), np.abs(darr['v(out)']))
        """
        with open(fname, "rb") as fp:
            plot = {}
            arrs = []
            plots = []
            while True:
                mdata = fp.readline(self.bsize_sp).split(b":", maxsplit=1)
                if len(mdata) == 2:
                    key = mdata[0].lower()
                    if key in self.mdata_list:
                        plot[key] = mdata[1].strip()
                    elif key == b"variables":
                        nvars = int(plot[b"no. variables"])
                        npoints = int(plot[b"no. points"])
                        varspecs = [
                            fp.readline(self.bsize_sp).strip().decode("ascii").split()
                            for _ in range(nvars)
                        ]
                        plot["varnames"], plot["varunits"] = zip(
                            *[(spec[1], spec[2]) for spec in varspecs]
                        )
                    elif key == b"binary":
                        rowdtype = np.dtype(
                            {
                                "names": plot["varnames"],
                                "formats": [
                                    np.complex_
                                    if b"complex" in plot[b"flags"]
                                    else np.float64
                                ]
                                * nvars,
                            }
                        )
                        arrs.append(np.fromfile(fp, dtype=rowdtype, count=npoints))
                        plots.append(plot.copy())
                        fp.readline()  # Read to the end of line
                else:
                    break
        return (arrs, plots)
