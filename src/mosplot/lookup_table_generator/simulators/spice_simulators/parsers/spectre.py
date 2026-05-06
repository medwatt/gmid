import numpy as np

bsize_sp = 512
mdata_list = [
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


def parse_file(fname):
    """
    Parse a Spectre nutbin (binary nutmeg) raw output file.

    The file is a sequence of plot blocks. Each block has ASCII header lines
    (Title, Date, Plotname, Flags, No. Variables, No. Points), a Variables
    listing (one variable per line: index, name, unit), and a Binary section
    of `npoints` rows, each holding `nvars` IEEE-754 double-precision values
    in big-endian byte order (or pairs of doubles when Flags contains
    "complex").

    A parametric sweep produces multiple plot blocks back-to-back.

    Returns
    -------
    (arrs, plots) : tuple
        arrs  : list of structured numpy arrays, one per plot block, with
                field names matching the variable names from the header.
        plots : list of metadata dicts, one per plot block.
    """
    with open(fname, "rb") as fp:
        plots = []
        arrs = []
        plot = {}
        nvars = 0
        npoints = 0
        while True:
            raw = fp.readline(bsize_sp)
            if not raw:
                break
            mdata = raw.split(b":", maxsplit=1)
            if len(mdata) != 2:
                continue
            key = mdata[0].lower().strip()
            if key in mdata_list:
                # A new Plotname before we finished the previous block means
                # the previous block had no data — drop it and start fresh.
                if key == b"plotname" and plot:
                    plot = {}
                    nvars = 0
                    npoints = 0
                plot[key] = mdata[1].strip()
            elif key == b"variables":
                nvars = int(plot[b"no. variables"])
                npoints = int(plot[b"no. points"])
                varspecs = []
                # Spectre puts the first variable spec on the same line as
                # "Variables:" itself; subsequent specs are on their own lines.
                first_parts = mdata[1].strip().decode("ascii", errors="replace").split()
                if len(first_parts) >= 2:
                    varspecs.append(first_parts)
                while len(varspecs) < nvars:
                    line = fp.readline(bsize_sp)
                    if not line:
                        break
                    parts = line.strip().decode("ascii", errors="replace").split()
                    if len(parts) >= 2:
                        varspecs.append(parts)
                plot["varnames"], plot["varunits"] = zip(
                    *[(spec[1], spec[2] if len(spec) > 2 else "") for spec in varspecs]
                )
            elif key == b"binary":
                is_complex = b"complex" in plot.get(b"flags", b"")
                fmt = ">c16" if is_complex else ">f8"
                rowdtype = np.dtype({
                    "names": list(plot["varnames"]),
                    "formats": [fmt] * nvars,
                })
                arr = np.fromfile(fp, dtype=rowdtype, count=npoints)
                arrs.append(arr)
                plots.append(plot.copy())
                plot = {}
                nvars = 0
                npoints = 0
        return (arrs, plots)
