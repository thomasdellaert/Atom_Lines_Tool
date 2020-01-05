from pandas import read_csv

# TODO: Potentially make it possible to also parse the ADS lines form to get observed transitions. In this case we'd also want to categorize them as E1, M2, etc
# TODO: Make a parser for the DREAM table that reconciles with NIST data. This will be tedious.
def fraction_to_float(frac):
    if "/" in frac:
        return float(frac.split("/")[0]) / float(frac.split("/")[1])
    else:
        return float(frac)

def term_frac(term):
    num = int(term / 0.5)
    if num % 2 != 0:
        return str(num) + "/2"
    else:
        return str(num / 2)

def read_term(term):
    """
    Take a term in LS, JK, or jj coupling and returns the good quantum numbers as floats.
    LS coupling: L, S, parity
    JK coupling: K, S, parity
    jj coupling: j1, j2, parity
    """

    am = ["S", "P", "D", "F", "G", "H"]
    parity = "e"
    L = None
    S = None
    K = None
    j1 = None
    j2 = None
    if "*" in term:
        parity = "o"
        term = term.strip("*")
    if "[" in term:
        S = (float(term[0]) - 1) / 2
        K = fraction_to_float(term[1:].strip("[]"))
    elif "(" in term:
        jtup = term.strip("()").split(",")
        j1 = fraction_to_float(jtup[0])
        j2 = fraction_to_float(jtup[1])
    elif term == "":
        pass
    else:
        L = am.index(term[1])
        S = (float(term[0]) - 1) / 2
    return L, S, K, j1, j2, parity

def parse_table(path):
    # Parse the input file that contains all the atomic data
    df = read_csv(path)
    df = df.iloc[:250]

    for col in df.keys():
        df[col] = df[col].str.strip("=\" ?")

    # Drop empty/useless columns
    try:
        df = df.drop(["Uncertainty (cm-1)", "Leading percentages", "Reference"], axis=1)
    except KeyError:
        pass

    # Convert columns that should be floats to floats
    for i, j in enumerate(df["Level (cm-1)"]):
        try:
            df["Level (cm-1)"][i] = float(df["Level (cm-1)"][i])
        except ValueError:
            df["Level (cm-1)"][i] = None
        try:
            df["Lande"][i] = float(df["Lande"][i])
        except ValueError:
            df["Lande"][i] = None
        try:
            df["J"][i] = fraction_to_float(df["J"][i])
        except ValueError:
            df["J"][i] = None

    # convert level in wavenumbers to level in terahertz
    df["Level (cm-1)"] *= 2.99792458e-2
    df.rename(columns={"Level (cm-1)": "Level (THz)"}, inplace=True)

    # based ion the terms, populate the dataset with numerical columns that specify the state
    for i in ["j2", "j1", "K", "L", "S"]:
        df.insert(2, i, None)
    df.insert(8, "Parity", None)
    for i, term in enumerate(df["Term"]):
        df["L"][i], df["S"][i], df["K"][i], df["j1"][i], df["j2"][i], df["Parity"][i] = read_term(term)
    return df
