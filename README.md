# Geometry Analysis

## How to install

```bash
git clone https://github.com/acaruso0/geometry_analysis.git
```

Append this line to your .bashrc:

```bash
export PYTHONPATH="${PYTHONPATH}:/path/to/geometry_analysis"
```

## How to use

### Import and file loading
The package is imported with:
```python
import geometry_analysis as geom
```

An .xyz file is imported with:
```python
xyzfile = geom.XYZfile('filename.xyz', (pbc_x, pbc_y, pbc_z))
```
where pbc_x, pbc_y and pbc_z are the lengths of the box (in Angstrom).

### Calculating rdfs
A rdf is calculated with:
```python
geom.rdf(xyzfile, 'atom1', 'atom2', initial_value(Angstrom), final_value(Angstrom), bins, 'output_file.dat')
```

### Examples
#### rdf calculation
A code to calculate a rdf would look like:
```python
import geometry_analysis as geom

xyzfile = geom.XYZfile("traj.xyz", (19.37, 19.37, 22.35))
geom.rdf(xyzfile, 'Cs', 'O', 1, 10, 300, "rdf_out.dat")
```
