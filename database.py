class Atom():
    def __init__(self, label, atnum, mass):
        self.label = label
        self.atnum = atnum
        self.mass = mass

def get_name(obj):
    if isinstance(obj, str):
        for atom in periodic_table:
            if periodic_table[atom].label == obj:
                return atom
    elif isinstance(obj, int):
        for atom in periodic_table:
            if periodic_table[atom].atnum == obj:
                return atom
    else:
        raise Exception("Argument must be the atomic label or the atomic number.")

periodic_table = {
        "hydrogen": Atom("H", 1, 1.008),
        "helium": Atom("He", 2, 4.0026),
        "lithium": Atom("Li", 3, 6.94),
        "beryllium": Atom("Be", 4, 9.0122),
        "boron": Atom("B", 5, 10.81),
        "carbon": Atom("C", 6, 12.011),
        "nitrogen": Atom("N", 7, 14.007),
        "oxygen": Atom("O", 8, 15.999),
        "fluorine": Atom("F", 9, 18.998),
        "neon": Atom("Ne", 10, 20.180),
        "sodium": Atom("Na", 11, 22.990),
        "magnesium": Atom("Mg", 12, 24.305),
        "aluminium": Atom("Al", 13, 26.982),
        "silicon": Atom("Si", 14, 28.085),
        "phoshporus": Atom("P", 15, 30.974),
        "sulfur": Atom("S", 16, 32.06),
        "chlorine": Atom("Cl", 17, 35.45),
        "caesium": Atom("Cs", 55, 132.91),
        }
